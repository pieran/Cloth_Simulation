#pragma once
#include <ncltech\GameObject.h>
#include <nclgl\Shader.h>
#include <nclgl\ObjMesh.h>
#include "SPHCommon.h"

class SkyBoxObject : public GameObject
{
public:
	SkyBoxObject(const std::string& name = "");
	virtual ~SkyBoxObject();

	inline GLuint GetTexture() { return m_cubeTex; }

protected:
	virtual void OnRenderObject();

protected:
	Shader* m_shdrCubeMap;
	GLuint m_cubeTex;

	OBJMesh* m_cubeMesh;
};