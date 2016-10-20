#include "GraphObject.h"
#include <ncltech\Scene.h>
#include <ncltech\NCLDebug.h>
#include <sstream>

GraphObject::GraphObject()
{
	this->m_BoundingRadius = FLT_MAX;

	glGenVertexArrays(1, &m_glArrayObj);
	glBindVertexArray(m_glArrayObj);
	glGenBuffers(1, &m_glBufferVerts);
	glBindVertexArray(0);

	m_GraphPosX = 10;
	m_GraphPosY = 10;

	m_ElapsedTime = 0.0f;

	m_Visible = true;
}

GraphObject::~GraphObject()
{
	if (m_glArrayObj)
	{
		glDeleteVertexArrays(1, &m_glArrayObj);
		glDeleteBuffers(1, &m_glBufferVerts);
	}
}

void GraphObject::OnAttachedToScene()
{
	GameObject::OnAttachedToScene();

	m_GraphPosX = m_Scene->GetWidth() - m_Params.graph_width - 10;
	m_GraphPosY = 10;
}

void GraphObject::OnRenderObject()
{
	GameObject::OnRenderObject();

	if (m_Visible)
	{
		Matrix4 invProjView = Matrix4::Inverse(m_Scene->GetProjViewMatrix());
		Matrix4 ortho = Matrix4::Orthographic(-1.0f, 1.0f, m_Scene->GetWidth(), 0, 0, m_Scene->GetHeight());
		Matrix4 transGraph = Matrix4::Translation(Vector3((float)m_GraphPosX, (float)m_GraphPosY, 0.0f));

		projTransform = ortho;
		invprojTransform = invProjView*projTransform;

		projTransformGraph = ortho*transGraph;
		invprojTransformGraph = invProjView*projTransformGraph;




		DrawBackground();
		DrawSeriesLabels();
		float min_time_value = m_ElapsedTime - m_Params.axis_x_span_time;
		for (uint idx = 0; idx < m_Series.size(); ++idx)
		{
			float depth = m_Params.series_depth - idx * 0.01f;
			m_Series[idx]->Prune(min_time_value);
			DrawSeries(*m_Series[idx], depth);
		}

		DrawXAxis();
		DrawYAxis();
		DrawGraphTitle();
	}
}

void GraphObject::DrawBackground()
{
	Vector3 pTL = invprojTransformGraph * Vector3(0, 0, m_Params.background_depth);
	Vector3 pTR = invprojTransformGraph * Vector3((float)m_Params.graph_width, 0, m_Params.background_depth);
	Vector3 pBL = invprojTransformGraph * Vector3(0, (float)m_Params.graph_height, m_Params.background_depth);
	Vector3 pBR = invprojTransformGraph * Vector3((float)m_Params.graph_width, (float)m_Params.graph_height, m_Params.background_depth);

	NCLDebug::DrawTriangle(pTL, pBR, pTR, m_Params.background_colour);
	NCLDebug::DrawTriangle(pTL, pBL, pBR, m_Params.background_colour);

	Vector3 peTL = invprojTransformGraph * Vector3(0, 0, m_Params.background_edge_depth);
	Vector3 peTR = invprojTransformGraph * Vector3((float)m_Params.graph_width, 0, m_Params.background_edge_depth);
	Vector3 peBL = invprojTransformGraph * Vector3(0, (float)m_Params.graph_height, m_Params.background_edge_depth);
	Vector3 peBR = invprojTransformGraph * Vector3((float)m_Params.graph_width, (float)m_Params.graph_height, m_Params.background_edge_depth);

	NCLDebug::DrawHairLine(peTL, peTR, m_Params.background_edge_colour);
	NCLDebug::DrawHairLine(peBL, peBR, m_Params.background_edge_colour);
	NCLDebug::DrawHairLine(peTL, peBL, m_Params.background_edge_colour);
	NCLDebug::DrawHairLine(peTR, peBR, m_Params.background_edge_colour);
}

void GraphObject::DrawXAxis()
{
	float time_end = m_ElapsedTime;
	float time_start = time_end - m_Params.axis_x_span_time;

	float time_step = (float)m_Params.graph_width / m_Params.axis_x_span_time * m_Params.axis_x_span_interval;

	//	float time_start_nearest = max(floor(time_start / axis_x_span_interval) * axis_x_span_interval, 0.0f);// -(time_start - floor(time_start));
	float lambda = fmod(time_end, m_Params.axis_x_span_interval);
	float time_start_nearest = max(time_start - lambda, 0.0f);
	float offset = lambda / m_Params.axis_x_span_interval;

	uint max_itr = (uint)floor((time_end - time_start_nearest) / m_Params.axis_x_span_interval);
	for (uint i = 0; i < max_itr; ++i)
	{
		float x_from_right = (offset + i) * time_step;
		float x_offset = m_Params.graph_width - x_from_right;


		Vector4 colour = m_Params.axis_x_label_colour;
		if (x_from_right < 10.0f)
			colour.w = min(x_from_right / 10.0f, 1.0f);
		else
			colour.w = min(x_offset / 20.0f, 1.0f);

		float time_multiple = ((max_itr - i) * m_Params.axis_x_span_interval) + time_start_nearest;
		std::stringstream ss; ss << time_multiple << "s";
		NCLDebug::DrawTextClipSpace(projTransformGraph * Vector4(x_offset, m_Params.graph_height + m_Params.axis_x_label_offset, 0, 1), m_Params.axis_x_label_fontsize, ss.str(), TEXTALIGN_CENTRE, colour);

		NCLDebug::DrawHairLine(invprojTransformGraph * Vector3(x_offset, 0.0f, m_Params.axis_x_depth), invprojTransformGraph * Vector3(x_offset, (float)m_Params.graph_height, m_Params.axis_x_depth), m_Params.axis_x_span_colour);
	}
}

void GraphObject::DrawYAxis()
{
	float step_size;
	uint max_itr, i;

	step_size = m_Params.graph_height / m_Params.axis_y_span_max * m_Params.axis_y_span_sub_interval;
	max_itr = ceil(m_Params.axis_y_span_max / m_Params.axis_y_span_sub_interval);
	for (i = 0; i <= max_itr; ++i)
	{
		float value = (i * m_Params.axis_y_span_label_interval);
		float y_offset = m_Params.graph_height - (i * step_size);

		if (y_offset > 5 && y_offset < m_Params.graph_height - 5)
			NCLDebug::DrawHairLine(invprojTransformGraph * Vector3(0.0f, y_offset, m_Params.axis_y_depth), invprojTransformGraph * Vector3((float)m_Params.graph_width, y_offset, m_Params.axis_y_depth), m_Params.axis_y_span_sub_colour);
	}


	step_size = m_Params.graph_height / m_Params.axis_y_span_max * m_Params.axis_y_span_label_interval;
	max_itr = ceil(m_Params.axis_y_span_max / m_Params.axis_y_span_label_interval);
	for (i = 0; i <= max_itr; ++i)
	{
		float value = (i * m_Params.axis_y_span_label_interval);
		float y_offset = m_Params.graph_height - (i * step_size);

		if (y_offset > 5 && y_offset < m_Params.graph_height - 5)
			NCLDebug::DrawHairLine(invprojTransformGraph * Vector3(0.0f, y_offset, m_Params.axis_y_depth), invprojTransformGraph * Vector3((float)m_Params.graph_width, y_offset, m_Params.axis_y_depth), m_Params.axis_y_span_main_colour);
	
		std::stringstream ss; ss << value << m_Params.axis_y_label_subfix;
		
		Vector4 label_pos = Vector4(m_Params.axis_y_label_offset, y_offset, 0, 1);

		if (i == 0)
			label_pos.y -= m_Params.axis_y_label_fontsize * 0.5f;
		else if (i == max_itr)
			label_pos.y += m_Params.axis_y_label_fontsize * 0.5f;

		NCLDebug::DrawTextClipSpace(projTransformGraph * label_pos, m_Params.axis_y_label_fontsize, ss.str(), TEXTALIGN_RIGHT, m_Params.axis_y_label_colour);
	}
}

void GraphObject::DrawSeriesLabels()
{
	float offset_y = m_Params.series_labels_fontsize * 0.5f;

	const float series_label_key_width  = 10.0f;
	const float series_label_key_height = 10.0f;

	uint column_idx = 0;
	uint row_idx = 0;

	float incr_x = m_Params.graph_width / float(m_Params.series_labels_numcolumns);


	float start_x = m_Scene->GetWidth() - m_Params.graph_width - 10;
	float start_y = m_GraphPosY + m_Params.axis_x_label_fontsize + 5.0f + m_Params.graph_height;

	for (uint idx = 0; idx < m_Series.size(); ++idx)
	{
		GraphSeries* series = m_Series[idx];
		uint label_len = series->Label().length();


		float pos_y = start_y + offset_y + row_idx * (m_Params.series_labels_fontsize + m_Params.series_labels_spacing);
		float pos_x = start_x + column_idx * incr_x;

		float box_x   = pos_x;
		float label_x = pos_x;
		float offset_box_y = m_Params.series_labels_fontsize / 2.0f;


		label_x += series_label_key_width + 10.0f;

		Vector3 tl = invprojTransform * Vector3(box_x, pos_y - offset_box_y, m_Params.series_depth);
		Vector3 tr = invprojTransform * Vector3(box_x + series_label_key_width, pos_y - offset_box_y, m_Params.series_depth);
		Vector3 bl = invprojTransform * Vector3(box_x, pos_y + series_label_key_height - offset_box_y, m_Params.series_depth);
		Vector3 br = invprojTransform * Vector3(box_x + series_label_key_width, pos_y + series_label_key_height - offset_box_y, m_Params.series_depth);

		NCLDebug::DrawTriangle(tl, br, tr, series->Colour());
		NCLDebug::DrawTriangle(tl, bl, br, series->Colour());

		NCLDebug::DrawTextClipSpace(projTransform * Vector4(label_x, pos_y, m_Params.series_depth, 1.0f), m_Params.series_labels_fontsize, series->Label(), TEXTALIGN_LEFT, m_Params.series_labels_colour);

		column_idx++;
		if (column_idx == m_Params.series_labels_numcolumns)
		{
			column_idx = 0;
			row_idx++;
		}
	}
}

void GraphObject::DrawGraphTitle()
{
	float pos_x = (float)m_GraphPosX + (float)m_Params.graph_width * 0.5f;
	float pos_y = (float)m_GraphPosY - m_Params.title_text_fontsize * 0.5f;

	NCLDebug::DrawTextClipSpace(projTransform * Vector4(pos_x, pos_y, m_Params.title_depth, 1.0f), m_Params.title_text_fontsize, m_Params.title_text, TEXTALIGN_CENTRE, m_Params.title_text_colour);
}

void GraphObject::DrawSeries(const GraphSeries& series, float depth)
{
	GraphSeriesIteratorConst itr = series.Begin();
	GraphSeriesIteratorConst end = series.End();

	if (itr != end)
	{
		Vector3 posLast, posNew;
		while (!GraphEntryToPosition(*itr, &posLast))
		{ 
			itr++;
			if (itr == end)
				return;	//Nothing in the series is within the graph timescale!
		}
		posLast.z = depth;
		itr++;
		for (; itr != end; itr++)
		{
			if (m_Params.axis_y_span_max_isdynamic && itr->data > m_Params.axis_y_span_max)
			{
				m_Params.axis_y_span_max = itr->data;
			}

			if (GraphEntryToPosition(*itr, &posNew))
			{
				posNew.z = depth;

				Vector3 posLBTrans = invprojTransformGraph * Vector3(posLast.x, m_Params.graph_height, depth);
				Vector3 posNBTrans = invprojTransformGraph * Vector3(posNew.x, m_Params.graph_height, depth);

				Vector3 posLTTrans = invprojTransformGraph * posLast;
				Vector3 posNTTrans = invprojTransformGraph * posNew;

				Vector4 fill_colour = series.Colour();
				NCLDebug::DrawTriangle(posLTTrans, posNBTrans, posNTTrans, fill_colour);
				NCLDebug::DrawTriangle(posLTTrans, posLBTrans, posNBTrans, fill_colour);

				posLast = posNew;
			}			
		}
	}
}

bool GraphObject::GraphEntryToPosition(const GraphEntry& entry, Vector3* out_pos)
{
	float start_time = m_ElapsedTime - m_Params.axis_x_span_time;
	float time_pos = entry.time - start_time;

	if (time_pos > m_Params.axis_x_span_time)
		return false;

	float x_pos = time_pos * m_Params.graph_width / m_Params.axis_x_span_time;
	float y_pos = m_Params.graph_height - entry.data * m_Params.graph_height / m_Params.axis_y_span_max;

	if (out_pos)
		*out_pos = Vector3(x_pos, y_pos, m_Params.series_depth);

	return true;
}

void GraphObject::OnUpdateObject(float dt)
{
	GameObject::OnUpdateObject(dt);
}

void GraphObject::BufferGLArrays()
{
/*	glBindVertexArray(m_glArrayObj);

	unsigned int numVertices = m_Cloth->m_NumTotal;
	unsigned int numIndices = m_Cloth->m_NumTriangles * 3;

	//Buffer vertex data	
	glBindBuffer(GL_ARRAY_BUFFER, m_glBufferVerts);
	glBufferData(GL_ARRAY_BUFFER, numVertices*sizeof(Vector3), &m_Cloth->m_PhyxelsPos[0], GL_STATIC_DRAW);
	glVertexAttribPointer(VERTEX_BUFFER, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(VERTEX_BUFFER);

	glBindVertexArray(0);*/
}


