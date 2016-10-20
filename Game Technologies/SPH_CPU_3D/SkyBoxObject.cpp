#include "SkyBoxObject.h"
#include <ncltech\Scene.h>


SkyBoxObject::SkyBoxObject(const std::string& name)
	: GameObject(name)
	, m_cubeTex(0)
	, m_shdrCubeMap(NULL)
	, m_cubeMesh(NULL)
{
	/*m_cubeTex = SOIL_load_OGL_cubemap(
		TEXTUREDIR"cubemap/plain_sky/plain_sky_back.jpg", TEXTUREDIR"cubemap/plain_sky/plain_sky_front.jpg",
		TEXTUREDIR"cubemap/plain_sky/plain_sky_top.jpg", TEXTUREDIR"cubemap/plain_sky/plain_sky_top.jpg",
		TEXTUREDIR"cubemap/plain_sky/plain_sky_right.jpg", TEXTUREDIR"cubemap/plain_sky/plain_sky_left.jpg",
		SOIL_LOAD_RGB,
		SOIL_CREATE_NEW_ID, 0);*/
	
	m_cubeTex = SOIL_load_OGL_cubemap(
		TEXTUREDIR"cubemap/fjaderholmarna/posx.jpg", TEXTUREDIR"cubemap/fjaderholmarna/negx.jpg",
		TEXTUREDIR"cubemap/fjaderholmarna/posy.jpg", TEXTUREDIR"cubemap/fjaderholmarna/negy.jpg",
		TEXTUREDIR"cubemap/fjaderholmarna/posz.jpg", TEXTUREDIR"cubemap/fjaderholmarna/negz.jpg",
		SOIL_LOAD_RGB,
		SOIL_CREATE_NEW_ID, 0);

	glEnable(GL_TEXTURE_CUBE_MAP);
	glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
	// Set parameters
	glBindTexture(GL_TEXTURE_CUBE_MAP, m_cubeTex);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_REPEAT);



	m_shdrCubeMap = new Shader(SHADERDIR"skybox.vert", SHADERDIR"skybox.frag");// , SHADERDIR"skybox.geom");
	m_shdrCubeMap->LinkProgram();


	m_cubeMesh = new OBJMesh(MESHDIR"cube.obj");
}


SkyBoxObject::~SkyBoxObject()
{
	if (m_shdrCubeMap)
	{
		delete m_shdrCubeMap;
		m_shdrCubeMap = NULL;
	}

	if (m_cubeTex)
	{
		glDeleteTextures(1, &m_cubeTex);
		m_cubeTex = 0;
	}
	
	if (m_cubeMesh)
	{
		delete m_cubeMesh;
		m_cubeMesh = NULL;
	}
}


void SkyBoxObject::OnRenderObject()
{
	//Matrix3 invViewMatrix = Matrix3::Inverse(Matrix3(m_Scene->GetViewMatrix()));

	Matrix4 viewMatrix = m_Scene->GetViewMatrix();
	viewMatrix.SetPositionVector(Vector3(0.0f, 0.0f, 0.0f));
	Matrix4 projViewMatrix = m_Scene->GetProjMatrix() * viewMatrix;

	glUseProgram(m_shdrCubeMap->GetProgram());
	glUniformMatrix4fv(glGetUniformLocation(m_shdrCubeMap->GetProgram(), "projViewMatrix"), 1, false, (float*)&projViewMatrix);
	glUniform1i(glGetUniformLocation(m_shdrCubeMap->GetProgram(), "texReflect"), 0);
	glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_CUBE_MAP, m_cubeTex);

	glDisable(GL_CULL_FACE);
	glDepthMask(GL_FALSE);
	m_cubeMesh->Draw();
	glDepthMask(GL_TRUE);
}
