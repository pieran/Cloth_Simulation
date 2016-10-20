#pragma once

#include <nclgl\OGLRenderer.h>
#include <ncltech\GameObject.h>
#include <nclgl\Mesh.h>
#include "MacGrid2D.h"


struct Pixel
{
	Vector3 wsPos;
	float pressure;
	Vector3 velocity;
};

class FluidRenderer : public GameObject
{
public:
	FluidRenderer(MacGrid2D* grid, uint view_x, uint view_y);
	~FluidRenderer();

protected:
	void OnRenderObject() override;				//Handles OpenGL calls to Render the object
	void OnUpdateObject(float dt) override;				//Override to handle things like AI etc on update loop
	void OnAttachedToScene() override;				//Called when attached to a scene object

protected:
	uint m_NumVertices;
	std::vector<Pixel> m_Pixels;
	bool m_RenderVelocityField;

	GLuint m_glArray;
	GLuint m_glPixelBuffer;

	GLuint m_glParticleTexture;

	Shader* m_PressureShader;

	MacGrid2D* m_FluidGrid;
};