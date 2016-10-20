#include "ClothRenderObject.h"
#include <ncltech\Scene.h>
#include <ncltech\NCLDebug.h>

ClothRenderObject::ClothRenderObject(CLOTH_SIM_TYPE* cloth) : m_Cloth(cloth)
{
	m_Paused = true;

	this->m_BoundingRadius = FLT_MAX;
	m_DrawElements = false;
	m_DrawStrain   = false;

	glGenVertexArrays(1, &m_glArrayObj);
	glBindVertexArray(m_glArrayObj);
	glGenBuffers(1, &m_glBufferVerts); 
	glGenBuffers(1, &m_glBufferColours);
	glGenBuffers(1, &m_glBufferNormals);
	glGenBuffers(1, &m_glBufferIndices);
	glGenBuffers(1, &m_glBufferTexCoords);
	glBindVertexArray(0);

	//m_glClothTexture = SOIL_load_OGL_texture(TEXTUREDIR"checkerboard.tga", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);

	m_glClothTexture = SOIL_load_OGL_texture(TEXTUREDIR"cloth_seamless_grey.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);

	glBindTexture(GL_TEXTURE_2D, m_glClothTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST); //No linear interpolation to get crisp checkerboard no matter the scalling
	glBindTexture(GL_TEXTURE_2D, 0);

	BuildGradientMap();

	m_RenderShader = new Shader(SHADERDIR"cloth/cloth_vertex.glsl", SHADERDIR"cloth/cloth_fragment_spec.glsl");
	if (!m_RenderShader->LinkProgram())
	{
		assert(false);
	}

	glUseProgram(m_RenderShader->GetProgram());
	glUniform1i(glGetUniformLocation(m_RenderShader->GetProgram(), "colorGradient"), 4);
}

ClothRenderObject::~ClothRenderObject()
{
	if (m_glArrayObj)
	{
		glDeleteVertexArrays(1, &m_glArrayObj);
		glDeleteBuffers(1, &m_glBufferVerts);
		glDeleteBuffers(1, &m_glBufferColours);
		glDeleteBuffers(1, &m_glBufferNormals);
		glDeleteBuffers(1, &m_glBufferIndices);
		glDeleteBuffers(1, &m_glBufferTexCoords);

		glDeleteTextures(1, &m_glClothTexture);
		glDeleteTextures(1, &m_glClothOutGradientMap);

		delete m_RenderShader;
	}
}

void ClothRenderObject::OnRenderObject()
{
	GameObject::OnRenderObject();

	if (m_Cloth)
	{
		GLboolean cullface;
		glGetBooleanv(GL_CULL_FACE, &cullface);

		GLint pid;
		glGetIntegerv(GL_CURRENT_PROGRAM, &pid);

		glDisable(GL_CULL_FACE);
		glUseProgram(m_RenderShader->GetProgram());

		Vector3 camPos = m_Scene->m_Camera->GetPosition();
		glUniformMatrix4fv(glGetUniformLocation(m_RenderShader->GetProgram(), "modelMatrix"), 1, false, (float*)&m_Scene->modelMatrix);
		glUniformMatrix4fv(glGetUniformLocation(m_RenderShader->GetProgram(), "viewMatrix"), 1, false, (float*)&m_Scene->viewMatrix);
		glUniformMatrix4fv(glGetUniformLocation(m_RenderShader->GetProgram(), "projMatrix"), 1, false, (float*)&m_Scene->projMatrix);
		glUniformMatrix4fv(glGetUniformLocation(m_RenderShader->GetProgram(), "textureMatrix"), 1, false, (float*)&m_Scene->textureMatrix);
		glUniform1i(glGetUniformLocation(m_RenderShader->GetProgram(), "diffuseTex"), 0);
		glUniform3fv(glGetUniformLocation(m_RenderShader->GetProgram(), "ambientColour"), 1, &m_Scene->m_AmbientColour.x);
		glUniform3fv(glGetUniformLocation(m_RenderShader->GetProgram(), "invLightDir"), 1, &m_Scene->m_InvLightDirection.x);
		glUniform3fv(glGetUniformLocation(m_RenderShader->GetProgram(), "cameraPos"), 1, &camPos.x);
		glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "specularIntensity"), m_Scene->m_SpecularIntensity);

		glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "isTextured"), m_DrawStrain ? 0.0f : 1.0f);
		glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "minVal"), m_Cloth->m_RenderValueMinMax.x);
		glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "maxVal"), m_Cloth->m_RenderValueMinMax.y);

		//NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Min/Max: %5.2f / %5.2f", m_Cloth->m_RenderValueMinMax.x, m_Cloth->m_RenderValueMinMax.y);
		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_1D, m_glClothOutGradientMap);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, m_glClothTexture);

		glBindVertexArray(m_glArrayObj);
		unsigned int numIndices = m_Cloth->m_NumTriangles * 3;
		glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		if (cullface)
			glEnable(GL_CULL_FACE);

		glUseProgram(pid);
	}
}

void ClothRenderObject::OnUpdateObject(float dt)
{
	GameObject::OnUpdateObject(dt);
	if (m_Cloth)
	{
		if (!m_Paused)
		{
			m_Cloth->Simulation_StepSimulation(1.0f / 60.0f);
		}

		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_R))
			m_Cloth->Simulation_StepSimulation(1.0f / 60.0f);

		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_X))
			m_DrawStrain = !m_DrawStrain;

		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_F))
			m_DrawElements = !m_DrawElements;

		if (m_DrawElements)
			m_Cloth->Render_DrawingToVisualiser();



		BufferGLArrays();
	}
}

void ClothRenderObject::BufferGLArrays()
{
	if (m_Cloth)
	{
		glBindVertexArray(m_glArrayObj);

		unsigned int numVertices = m_Cloth->m_NumTotal;
		unsigned int numIndices = m_Cloth->m_NumTriangles * 3;

		//Buffer vertex data	
		glBindBuffer(GL_ARRAY_BUFFER, m_glBufferVerts);
		glBufferData(GL_ARRAY_BUFFER, numVertices*sizeof(Vector3), &m_Cloth->m_PhyxelsPos[0], GL_STATIC_DRAW);
		glVertexAttribPointer(VERTEX_BUFFER, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(VERTEX_BUFFER);


		//Buffer vertex colour data	
		glBindBuffer(GL_ARRAY_BUFFER, m_glBufferColours);
		glBufferData(GL_ARRAY_BUFFER, numVertices*sizeof(float), &m_Cloth->m_RenderValue[0], GL_STATIC_DRAW);
		glVertexAttribPointer(COLOUR_BUFFER, 1, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(COLOUR_BUFFER);
		

		//Buffer texture coord data
		glBindBuffer(GL_ARRAY_BUFFER, m_glBufferTexCoords);
		glBufferData(GL_ARRAY_BUFFER, numVertices*sizeof(Vector2), &m_Cloth->m_PhyxelTexCoords[0], GL_STATIC_DRAW);
		glVertexAttribPointer(TEXTURE_BUFFER, 2, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(TEXTURE_BUFFER);

		//Buffer normal data
		glBindBuffer(GL_ARRAY_BUFFER, m_glBufferNormals);
		glBufferData(GL_ARRAY_BUFFER, numVertices*sizeof(Vector3), &m_Cloth->m_PhyxelNormals[0], GL_STATIC_DRAW);
		glVertexAttribPointer(NORMAL_BUFFER, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(NORMAL_BUFFER);

		//buffer index data
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_glBufferIndices);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, numIndices*sizeof(GLuint), &m_Cloth->m_RenderIndices[0], GL_STATIC_DRAW);

		glBindVertexArray(0);
	}
}

void ClothRenderObject::BuildGradientMap()
{	
	const uint num_vals = 360;



	glGenTextures(1, &m_glClothOutGradientMap);
	glBindTexture(GL_TEXTURE_1D, m_glClothOutGradientMap);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	Vector3* gradient = new Vector3[num_vals];
	float step = 1.0f / num_vals * 0.25f;
	for (uint i = 0; i < num_vals; ++i)
	{
		hsv2rgb(gradient[i], Vector3(i * 0.25f + 270.0f, 1.0f, 1.0f));
	}

	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB32F, num_vals, 0, GL_RGB, GL_FLOAT, gradient);

	delete[] gradient;

	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_1D, 0);
}