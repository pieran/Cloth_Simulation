#pragma once
#include <nclgl\Mesh.h>
#include <ncltech\GameObject.h>
#include "PositionBasedDynamics.h"

#define CLOTH_UPDATETIMESTEP (1.0f / 240.0f)

class ClothRenderObject : public GameObject
{
public:
	ClothRenderObject(PositionBasedDynamicsMS* cloth);
	virtual ~ClothRenderObject();

	inline void SetPaused(bool paused) { m_Paused = paused; }
	inline void TogglePaused() { m_Paused = !m_Paused; }

protected:
	void OnRenderObject() override;				//Handles OpenGL calls to Render the object
	void OnUpdateObject(float dt) override;		//Override to handle things like AI etc on update loop

	void BufferGLArrays();

protected:
	bool m_Paused;

	PositionBasedDynamicsMS* m_Cloth;

	GLuint m_glClothTexture;

	GLuint m_glArrayObj;
	GLuint m_glBufferVerts;
	GLuint m_glBufferTexCoords;
	GLuint m_glBufferNormals;
	GLuint m_glBufferIndices;
};

