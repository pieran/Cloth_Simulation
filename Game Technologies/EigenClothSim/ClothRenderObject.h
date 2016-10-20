#pragma once

#define MODE 2
#include <nclgl\Shader.h>
#include <nclgl\Mesh.h>
#include <ncltech\GameObject.h>

#if (MODE==0)
#include "FEM3Noded.h"
#define CLOTH_SIM_TYPE FEM3Noded
#else
	#if (MODE==1)
	#include "FEM3NodedLMA.h"
	#define CLOTH_SIM_TYPE FEM3NodedLMA
	#else 
		#if (MODE==2)
		#include "FEM3NodedFAST.h"
		#define CLOTH_SIM_TYPE FEM3NodedFast
		#else
		#include "FEM6Noded.h"
		#define CLOTH_SIM_TYPE FEM6Noded
		#endif
	#endif
#endif


class ClothRenderObject : public GameObject
{
public:
	ClothRenderObject(CLOTH_SIM_TYPE* cloth);

	virtual ~ClothRenderObject();

	inline void SetPaused(bool paused) { m_Paused = paused; }
	inline void TogglePaused() { m_Paused = !m_Paused; }
	inline bool IsPaused() const { return m_Paused; }

	CLOTH_SIM_TYPE* m_Cloth;

protected:
	void OnRenderObject() override;				//Handles OpenGL calls to Render the object
	void OnUpdateObject(float dt) override;		//Override to handle things like AI etc on update loop

	void BufferGLArrays();
	void BuildGradientMap();

protected:
	bool m_Paused;
	bool m_DrawElements;
	bool m_DrawStrain;

	Shader* m_RenderShader;
	GLuint  m_glClothOutGradientMap;

	GLuint m_glClothTexture;

	GLuint m_glArrayObj;
	GLuint m_glBufferVerts;
	GLuint m_glBufferColours;
	GLuint m_glBufferTexCoords;
	GLuint m_glBufferNormals;
	GLuint m_glBufferIndices;
};

