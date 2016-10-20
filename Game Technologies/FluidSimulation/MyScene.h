
#pragma once

#include <nclgl\Mesh.h>
#include <ncltech\Scene.h>
#include <ncltech\SimpleMeshObject.h>

#include "FluidSim2D.h"
#include "FluidRenderer.h"

class MyScene : public Scene
{
public:
	MyScene(Window& window);
	~MyScene();

	bool InitialiseGL()	override;
	void RenderScene() override;
	void UpdateScene(float dt)  override;

protected:
	FluidSim2D* m_FluidSimulation;
	FluidRenderer* m_FluidRenderer;
};