
#pragma once

#include <nclgl\Mesh.h>
#include <ncltech\Scene.h>
#include <ncltech\SimpleMeshObject.h>
#include "sph_system.h"
#include "SPH_GPU\sph3d_system.h"
#include "FluidScreenSpaceRenderer.h"

//#define SPH_SYSTEM SPHSystem
#define USE_GPU TRUE


class MyScene : public Scene
{
public:
	MyScene(Window& window);
	~MyScene();

	bool InitialiseGL()	override;
	void RenderScene() override;
	void UpdateScene(float dt)  override;

protected:
	FluidScreenSpaceRenderer* m_fluidRenderer;

	Vector3 real_world_origin;
	Vector3 real_world_side;
	Vector3 sim_ratio;

#if USE_GPU == TRUE
	SPH3DSystem* sph;
#else
	SPHSystem* sph;
#endif
};