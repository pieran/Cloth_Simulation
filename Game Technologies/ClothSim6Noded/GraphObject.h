#pragma once

#include <nclgl\Mesh.h>
#include <ncltech\GameObject.h>
#include "EigenDefines.h"
#include <map>

class GraphObject;

struct GraphEntry
{
	GraphEntry() : time(0.0f), data(0.0f) {}
	GraphEntry(float entry_time, float entry_data) : time(entry_time), data(entry_data) {}

	float time;
	float data;
};


typedef std::vector<GraphEntry>::iterator GraphSeriesIterator;
typedef std::vector<GraphEntry>::const_iterator GraphSeriesIteratorConst;
class GraphSeries
{
public:
	GraphSeries(GraphObject* parent_graph, const std::string& label, const Vector4& colour) : m_Graph(parent_graph), m_Colour(colour), m_Label(label) {}
	virtual ~GraphSeries() {}

	inline std::string& Label() { return m_Label; }
	inline Vector4& Colour() { return m_Colour; }

	inline const std::string& Label() const { return m_Label; }
	inline const Vector4& Colour() const { return m_Colour; }

	inline void AddDataEntry(float time, float value) { m_Values.push_back(GraphEntry(time, value)); }
	inline void AddDataEntry(const GraphEntry& entry) { m_Values.push_back(entry); }

	inline GraphSeriesIterator Begin() { return m_Values.begin(); }
	inline GraphSeriesIterator End() { return m_Values.end(); }

	inline GraphSeriesIteratorConst Begin()	const { return m_Values.begin(); }
	inline GraphSeriesIteratorConst End() const { return m_Values.end(); }

	void Prune(float min_time_value)
	{
		GraphSeriesIterator itr = m_Values.begin();
		while (itr != m_Values.end() && itr->time < min_time_value)
		{
			itr = m_Values.erase(itr);
		}
	}
protected:
	GraphObject* m_Graph;

	Vector4		 m_Colour;
	std::string  m_Label;
	std::vector<GraphEntry> m_Values;
};

enum SeriesAlignment
{
	SERIES_ALIGN_TOP = 1,
	SERIES_ALIGN_LEFT = 2,


	SERIES_ALIGN_BOTTOMRIGHT = 0,
	SERIES_ALIGN_TOPRIGHT = SERIES_ALIGN_TOP,
	SERIES_ALIGN_BOTTOMLEFT = SERIES_ALIGN_LEFT,
	SERIES_ALIGN_TOPLEFT = SERIES_ALIGN_TOP | SERIES_ALIGN_LEFT
};

enum GraphSectionPositioning
{
	GRAPH_POSITION_HIDDEN = 0,
	GRAPH_POSITION_TOP = 1,
	GRAPH_POSITION_BOTTOM = 2,
	GRAPH_POSITION_LEFT = 3,
	GRAPH_POSITION_RIGHT = 4
};


struct GraphParameters
{
	GraphParameters() {}
	float	  graph_width = 300.0f;
	float	  graph_height = 150.0f;

	Vector4   background_colour = Vector4(0.7f, 0.7f, 0.7f, 0.9f);
	Vector4   background_edge_colour = Vector4(0.2f, 0.2f, 0.2f, 1.0f);
	float	  background_depth = 0.8f;
	float	  background_edge_depth = 0.85f;

	Vector4   axis_x_label_colour = Vector4(1.0f, 1.0f, 1.0f, 1.0f);
	float	  axis_x_label_fontsize = 12.0f;
	float	  axis_x_label_offset = axis_x_label_fontsize * 0.5f + 2.0f;

	Vector4   axis_x_span_colour = Vector4(0.4f, 0.4f, 0.4f, 0.5f);
	float	  axis_x_span_time = 4.0f;
	float	  axis_x_span_interval = 1.0f;
	float	  axis_x_depth = 0.95f;

	Vector4   axis_y_label_colour = Vector4(1.0f, 1.0f, 1.0f, 1.0f);
	float	  axis_y_label_fontsize = 12.0f;
	float	  axis_y_label_offset = -2.0f;
	std::string axis_y_label_subfix = "%";

	Vector4   axis_y_span_main_colour = Vector4(0.2f, 0.2f, 0.2f, 0.5f);
	Vector4   axis_y_span_sub_colour = Vector4(0.4f, 0.4f, 0.4f, 0.5f);
	float	  axis_y_span_max = 100.0f;
	bool	  axis_y_span_max_isdynamic = false;
	float	  axis_y_span_label_interval = 50.0f;
	float	  axis_y_span_sub_interval = 10.0f;
	float	  axis_y_depth = 1.0f;

	Vector4   series_labels_colour = Vector4(1.0f, 1.0f, 1.0f, 1.0f);
	float	  series_labels_fontsize = 12.0f;
	float     series_labels_spacing = 5.0f;
	uint	  series_labels_numcolumns = 2;
	float     series_depth = 0.9f;


	std::string	title_text = "";
	Vector4		title_text_colour = Vector4(1.0f, 1.0f, 1.0f, 1.0f);
	float		title_text_fontsize = 16.0f;
	float	    title_depth = 1.0f;

};

class GraphObject : public GameObject
{
public:
	GraphObject();
	virtual ~GraphObject();

	void UpdateYScale(float dt) { m_ElapsedTime += dt; }
	void AddSeries(const std::string& label, const Vector4& colour) { m_Series.push_back(new GraphSeries(this, label, colour)); }

	GraphSeries& GetSeries(uint idx) { return *m_Series[idx]; }
	const GraphSeries& GetSeries(uint idx) const { return *m_Series[idx]; }

	inline void AddDataEntry(uint seriesIdx, float value) {
		m_Series[seriesIdx]->AddDataEntry(m_ElapsedTime, value);
		if (m_Params.axis_y_span_max_isdynamic && value > m_Params.axis_y_span_max)
		{
			m_Params.axis_y_span_max = value;
		}
	}
	inline void SetVisibility(bool visible) { m_Visible = visible; }

	inline GraphParameters& Parameters() { return m_Params; }
	inline const GraphParameters& Parameters() const { return m_Params; }

	void SetPosition(uint x, uint y)
	{
		m_GraphPosX = x - m_Params.graph_width;
		m_GraphPosY = y;
	}

protected:
	void OnAttachedToScene() override;			//Called when attached to a scene object
	void OnRenderObject() override;				//Handles OpenGL calls to Render the object
	void OnUpdateObject(float dt) override;		//Override to handle things like AI etc on update loop

	void BufferGLArrays();

	void DrawBackground();			//Background
	void DrawXAxis();				//X Axis intervals + Labels
	void DrawYAxis();				//Y Axis intervals + Labels
	void DrawSeriesLabels();		//Series Labels & Colour-Key
	void DrawGraphTitle();

	void DrawSeries(const GraphSeries& series, float depth);
	bool GraphEntryToPosition(const GraphEntry& entry, Vector3* out_pos);	//returns false if outside of area

protected:
	GraphParameters m_Params;

	bool m_Visible;
	uint m_GraphPosX, m_GraphPosY;


	float m_ElapsedTime;

	Matrix4 projTransform;
	Matrix4 invprojTransform;
	Matrix4 projTransformGraph;
	Matrix4 invprojTransformGraph;


	std::vector<GraphSeries*> m_Series;

	GLuint m_glArrayObj;
	GLuint m_glBufferVerts;
};

