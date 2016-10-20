#include "FluidRenderer.h"
#include <ncltech\Scene.h>
#include <assert.h>
#include <algorithm>
#include <ncltech\NCLDebug.h>

FluidRenderer::FluidRenderer(MacGrid2D* grid, uint view_x, uint view_y)
	: GameObject("FluidRenderer")
	, m_glArray(NULL)
	, m_PressureShader(NULL)
	, m_FluidGrid(grid)
	, m_RenderVelocityField(true)
{
	this->SetBoundingRadius(FLT_MAX);

	//Load Particle Texture
	m_glParticleTexture = SOIL_load_OGL_texture(TEXTUREDIR"smoke_particle.png", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);
	glBindTexture(GL_TEXTURE_2D, m_glParticleTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST); //No linear interpolation to get crisp checkerboard no matter the scalling
	glBindTexture(GL_TEXTURE_2D, 0);


	//Load Shader
	m_PressureShader = new Shader(SHADERDIR"fluid/MACParticleVertex.glsl", SHADERDIR"fluid/MACParticleFragment.glsl", SHADERDIR"fluid/MACParticleGeometry.glsl");
	assert(m_PressureShader->LinkProgram());

	//Build Buffers
	int sx = view_x;
	int sy = view_y;

	float stepx = m_FluidGrid->m_DivX / (float)view_x;
	float stepy = m_FluidGrid->m_DivY / (float)view_y;

	m_NumVertices = sx * sy;
	m_Pixels.resize(m_NumVertices);

#pragma omp parallel for
	for (int y = 0; y < sy; ++y)
	{
		for (int x = 0; x < sx; ++x)
		{
			float fx = ((float)x + 0.5f) * stepx;
			float fy = ((float)y + 0.5f) * stepy;

			Vector3 ws_pos;
			ws_pos.x = m_FluidGrid->m_CellDimensions.x * fx + m_FluidGrid->m_BoundingRegion.offset_min.x;
			ws_pos.y = m_FluidGrid->m_CellDimensions.y * fy + m_FluidGrid->m_BoundingRegion.offset_min.y;
			ws_pos.z = 0.0f;

			int cellIdx = y * sx + x;
			m_Pixels[cellIdx].wsPos = ws_pos;
			m_Pixels[cellIdx].pressure = 0.0f;
			m_Pixels[cellIdx].velocity = Vector3(0.0f, 0.0f, 0.0f);
		}
	}
	

	glGenVertexArrays(1, &m_glArray);
	glBindVertexArray(m_glArray);

	glGenBuffers(1, &m_glPixelBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, m_glPixelBuffer);
	glBufferData(GL_ARRAY_BUFFER, m_NumVertices*sizeof(Pixel), &m_Pixels[0], GL_STREAM_DRAW);
	glVertexAttribPointer(VERTEX_BUFFER, 3, GL_FLOAT, GL_FALSE, sizeof(Pixel), 0);
	glVertexAttribPointer(COLOUR_BUFFER, 1, GL_FLOAT, GL_FALSE, sizeof(Pixel), (GLvoid*)12);
	glVertexAttribPointer(NORMAL_BUFFER, 3, GL_FLOAT, GL_FALSE, sizeof(Pixel), (GLvoid*)16);
	glEnableVertexAttribArray(VERTEX_BUFFER);
	glEnableVertexAttribArray(COLOUR_BUFFER);
	glEnableVertexAttribArray(NORMAL_BUFFER);

	glBindVertexArray(0);
}

FluidRenderer::~FluidRenderer()
{
	m_FluidGrid = NULL;

	if (m_glArray)
	{
		glDeleteBuffers(1, &m_glPixelBuffer);
		glDeleteVertexArrays(1, &m_glArray);
		m_glArray = NULL;
	}

	if (m_PressureShader)
	{
		delete m_PressureShader;
		m_PressureShader = NULL;
	}
}

void FluidRenderer::OnRenderObject()
{
#pragma omp parallel for
	for (int i = 0; i < m_NumVertices; ++i)
	{
		Vector2 pos = Vector2(m_Pixels[i].wsPos.x, m_Pixels[i].wsPos.y);
		Vector2 vel = m_FluidGrid->GetVelocityTrilinear(pos);
		m_Pixels[i].pressure = m_FluidGrid->GetPressure(pos);
		m_Pixels[i].velocity = Vector3(vel.x, vel.y, 0.0f);
	}

	Vector3 camPos = m_Scene->m_Camera->GetPosition();

	/*auto compareSort = [&camPos](const Pixel& a, const Pixel& b)
	{
		float disA = (a.wsPos - camPos).LengthSquared();
		float disB = (b.wsPos - camPos).LengthSquared();
		return disA > disB;
	};
	std::sort(m_Pixels.begin(), m_Pixels.end(), compareSort);*/
	const float e_offset = -1E-1;

	if (m_RenderVelocityField)
	{		
		for (int i = 0; i < m_NumVertices; ++i)
		{
			/*Vector4 col = Vector4(m_Pixels[i].pressure, m_Pixels[i].pressure * 0.5f, 0.0f, 1.0f);
			Vector3 vel = m_Pixels[i].velocity;
			vel.Normalise();
			Vector3 end = m_Pixels[i].wsPos + vel * 0.5f;
			NCLDebug::DrawHairLine(m_Pixels[i].wsPos, end, col);*/
			//NCLDebug::DrawPoint(m_Pixels[i].wsPos, 0.05f, col);


			Vector3 vel = m_Pixels[i].velocity;
			Vector4 col = Vector4(fabs(vel.x), fabs(vel.y), fabs(vel.z), 1.0f);
			//Vector4 col = Vector4(vel.x, vel.y, vel.z, 2.0f) * 0.5f + 0.5f;
			NCLDebug::DrawPoint(m_Pixels[i].wsPos + Vector3(0.0f, 0.0f, e_offset*2.0f), 0.05f, col);
		}
	}
	else
	{
		for (int i = 0; i < m_NumVertices; ++i)
		{
			float pressure = m_Pixels[i].pressure;
			Vector4 col = Vector4(pressure, pressure, pressure, 1.0f);
			NCLDebug::DrawPoint(m_Pixels[i].wsPos + Vector3(0.0f, 0.0f, e_offset*2.0f), 0.05f, col);
		}
	}

	//Draw Center
	Vector3 center;
	center.x = m_FluidGrid->m_BoundingRegion.offset_min.x + m_FluidGrid->m_BoundingRegion.scale.x * 0.5f;
	center.y = m_FluidGrid->m_BoundingRegion.offset_min.y + m_FluidGrid->m_BoundingRegion.scale.y * 0.5f;
	center.z = 0.0f;
	NCLDebug::DrawPoint(center, 0.1f, Vector4(1.0f, 0.0f, 1.0f, 1.0f));


	//Draw Border
	const Vector4 border_col = Vector4(0.0f, 0.0f, 0.0f, 1.0f);
	const float border_width = 0.05f;

	BoundingBox2D& region = m_FluidGrid->m_BoundingRegion;
	Vector3 a = Vector3(region.offset_min.x, region.offset_min.y, e_offset);
	Vector3 b = a + Vector3(region.scale.x, region.scale.y, 0.0f);

	NCLDebug::DrawThickLine(Vector3(a.x, a.y, a.z), Vector3(b.x, a.y, a.z), border_width, border_col);
	NCLDebug::DrawThickLine(Vector3(a.x, b.y, a.z), Vector3(b.x, b.y, a.z), border_width, border_col);

	NCLDebug::DrawThickLine(Vector3(a.x, a.y, a.z), Vector3(a.x, b.y, a.z), border_width, border_col);
	NCLDebug::DrawThickLine(Vector3(b.x, a.y, a.z), Vector3(b.x, b.y, a.z), border_width, border_col);

	NCLDebug::DrawThickLine(Vector3(a.x, a.y, a.z), Vector3(b.x, b.y, a.z), border_width, border_col);
	NCLDebug::DrawThickLine(Vector3(b.x, a.y, a.z), Vector3(a.x, b.y, a.z), border_width, border_col);

	/*Matrix4 mvp = m_Scene->GetProjViewMatrix(); //No model matrix as GRID is already in world space (probably)

	Matrix3 invview = Matrix3::Inverse(Matrix3(m_Scene->GetViewMatrix()));
	float size = 3.0f / m_SubDivisions;
	Vector3 camLeft = invview * Vector3(size, 0.0f, 0.0f);
	Vector3 camUp = invview * Vector3(0.0f, size, 0.0f);

	glUseProgram(m_PressureShader->GetProgram());
	glUniformMatrix4fv(glGetUniformLocation(m_PressureShader->GetProgram(), "mvpMatrix"), 1, false, (float*)&mvp);
	glUniform3fv(glGetUniformLocation(m_PressureShader->GetProgram(), "camLeft"), 1, (float*)&camLeft.x);
	glUniform3fv(glGetUniformLocation(m_PressureShader->GetProgram(), "camUp"), 1, (float*)&camUp.x);
	glUniform1i(glGetUniformLocation(m_PressureShader->GetProgram(), "particleTex"), 4);

	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D, m_glParticleTexture);

	//Draw Cells
	glDisable(GL_DEPTH_TEST);
	glBindVertexArray(m_glArray);
	glBindBuffer(GL_ARRAY_BUFFER, m_glPixelBuffer);
	glBufferSubData(GL_ARRAY_BUFFER, 0, m_NumVertices*sizeof(Pixel), &m_Pixels[0]);
	glDrawArrays(GL_POINTS, 0, m_NumVertices);
	glBindVertexArray(0);
	glEnable(GL_DEPTH_TEST);*/
}

void FluidRenderer::OnUpdateObject(float dt)
{
	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_F))
		m_RenderVelocityField = !m_RenderVelocityField;
}

void FluidRenderer::OnAttachedToScene()
{

}