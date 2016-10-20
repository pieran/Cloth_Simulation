
#include "MyScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\NCLDebug.h>
#include "SkyBoxObject.h"

#include <cuda_gl_interop.h>

bool init_cuda(void)
{
	int count = 0;
	int i = 0;

	cudaGetDeviceCount(&count);
	if (count == 0)
	{
		fprintf(stderr, "There is no device.\n");
		return false;
	}

	for (i = 0; i < count; i++)
	{
		cudaDeviceProp prop;
		if (cudaGetDeviceProperties(&prop, i) == cudaSuccess)
		{
			if (prop.major >= 1)
			{
				break;
			}
		}
	}

	if (i == count)
	{
		fprintf(stderr, "There is no device supporting CUDA.\n");
		return false;
	}

	cudaSetDevice(i);

	printf("CUDA initialized.\n");
	return true;
}

MyScene::MyScene(Window& window) : Scene(window)
{
	if (init == true)
		init = InitialiseGL();

	projMatrix = Matrix4::Perspective(1.0f, 1000.0f, width / (float)height, 45.0f);
	UpdateWorldMatrices(m_RootGameObject, Matrix4());


	
}

MyScene::~MyScene()
{
	delete sph;
}

bool MyScene::InitialiseGL()
{
	init_cuda();

	m_Camera->SetPosition(Vector3(-6.25f, 2.0f, 10.0f));

	PhysicsEngine::Instance()->SetGravity(Vector3(0.0f, 0.0f, 0.0f));		//No Gravity
	PhysicsEngine::Instance()->SetDampingFactor(1.0f);						//No Damping

	//Create Skybox
	SkyBoxObject* sky = new SkyBoxObject("Skybox");
	sky->SetBoundingRadius(FLT_MAX);
	this->AddGameObject(sky);

	//Create Ground
	/*SimpleMeshObject* ground = new SimpleMeshObject("Ground");
	ground->SetMesh(CommonMeshes::Cube(), false);
	ground->SetLocalTransform(Matrix4::Scale(Vector3(80.0f, 0.1f, 2.f)));
	ground->SetColour(Vector4(0.2f, 1.0f, 0.5f, 1.0f));
	ground->SetBoundingRadius(80.0f * 80.f);

	ground->Physics()->SetPosition(Vector3(-6.25f, -0.2f, 0.0f));

	this->AddGameObject(ground);*/



	real_world_origin.x = -10.0f;
	real_world_origin.y = -1.0f;
	real_world_origin.z = -10.0f;

	real_world_side.x = 20.0f;
	real_world_side.y = 20.0f;
	real_world_side.z = 20.0f;

#if USE_GPU == TRUE
	sph = new SPH3DSystem();
	sph->init_system();

	sim_ratio.x = real_world_side.x / sph->hParam->world_size.x;
	sim_ratio.y = real_world_side.y / sph->hParam->world_size.y;
	sim_ratio.z = real_world_side.z / sph->hParam->world_size.z;
#else


	sph = new SPH_SYSTEM();
	sph->init_system();

	sim_ratio.x = real_world_side.x / sph->world_size.x;
	sim_ratio.y = real_world_side.y / sph->world_size.y;
	sim_ratio.z = real_world_side.z / sph->world_size.z;
#endif

	m_fluidRenderer = new FluidScreenSpaceRenderer();
	m_fluidRenderer->SetWorldSpace(real_world_origin, sim_ratio);
#if USE_GPU == TRUE
	m_fluidRenderer->BuildBuffers(sph->hParam->max_particle);
#else
	m_fluidRenderer->BuildBuffers(sph->max_particle);
#endif
	m_fluidRenderer->SetBackgroundSceneCubemap(sky->GetTexture());
	this->AddGameObject(m_fluidRenderer);

	return true;
}

void MyScene::UpdateScene(float msec)
{
	Scene::UpdateScene(msec);
	sph->animation();

	sph->sys_running = Window::GetKeyboard()->KeyDown(KEYBOARD_R);
#if USE_GPU == TRUE
	m_fluidRenderer->UpdateBuffers(sph->hParam->num_particle, sph->hPoints);
#else
	m_fluidRenderer->UpdateBuffers(sph->num_particle, sph->mem);
#endif
}

void MyScene::RenderScene()
{
	/*for (uint i = 0; i < sph->num_particle; i++)
	{
		Vector3 pos = Vector3(sph->mem[i].pos.x*sim_ratio.x + real_world_origin.x,
			sph->mem[i].pos.y*sim_ratio.y + real_world_origin.y,
			sph->mem[i].pos.z*sim_ratio.z + real_world_origin.z);
		NCLDebug::DrawPoint(pos, 0.1f, Vector4(0.2f, 0.2f, 1.0f, 1.0f));
	}*/

	Scene::RenderScene();
}

