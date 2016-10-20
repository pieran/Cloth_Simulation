#include "ClothRenderObject.h"

ClothRenderObject::ClothRenderObject(PositionBasedDynamicsMS* cloth) : m_Cloth(cloth)
{
	m_Paused = true;

	this->m_BoundingRadius = FLT_MAX;

	glGenVertexArrays(1, &m_glArrayObj);
	glBindVertexArray(m_glArrayObj);
	glGenBuffers(1, &m_glBufferVerts);
	glGenBuffers(1, &m_glBufferNormals);
	glGenBuffers(1, &m_glBufferIndices); 
	glGenBuffers(1, &m_glBufferTexCoords);
	glBindVertexArray(0);

	m_glClothTexture = SOIL_load_OGL_texture(TEXTUREDIR"checkerboard.tga", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);
	glBindTexture(GL_TEXTURE_2D, m_glClothTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST); //No linear interpolation to get crisp checkerboard no matter the scalling
	glBindTexture(GL_TEXTURE_2D, 0);
}

ClothRenderObject::~ClothRenderObject()
{
	if (m_glArrayObj)
	{
		glDeleteVertexArrays(1, &m_glArrayObj);
		glDeleteBuffers(1, &m_glBufferVerts);
		glDeleteBuffers(1, &m_glBufferNormals);
		glDeleteBuffers(1, &m_glBufferIndices); 
		glDeleteBuffers(1, &m_glBufferTexCoords);

		glDeleteTextures(1, &m_glClothTexture);
	}
}

void ClothRenderObject::OnRenderObject()
{
	GameObject::OnRenderObject();

	if (m_Cloth)
	{
		GLboolean cullface;
		glGetBooleanv(GL_CULL_FACE, &cullface);

		glDisable(GL_CULL_FACE);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, m_glClothTexture);

		glBindVertexArray(m_glArrayObj);
		unsigned int numIndices = m_Cloth->m_NumTriangles * 3;
		glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		if (cullface)
			glEnable(GL_CULL_FACE);
	}
}

void ClothRenderObject::OnUpdateObject(float dt)
{
	GameObject::OnUpdateObject(dt);
	if (m_Cloth)
	{
		if (!m_Paused)		
			m_Cloth->Simulation_StepSimulation(CLOTH_UPDATETIMESTEP);
		
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
		glBufferData(GL_ARRAY_BUFFER, numVertices*sizeof(Vector3), &m_Cloth->m_PhyxelsPosNew[0], GL_STATIC_DRAW);
		glVertexAttribPointer(VERTEX_BUFFER, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(VERTEX_BUFFER);


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
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, numIndices*sizeof(GLuint), &m_Cloth->m_Triangles[0], GL_STATIC_DRAW);

		glBindVertexArray(0);
	}
}


