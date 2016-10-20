
#include "MyScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\NCLDebug.h>


MyScene::MyScene(Window& window) : Scene(window), m_FluidSimulation(NULL)
{
	if (init == true)
		init = InitialiseGL();

	UpdateWorldMatrices(m_RootGameObject, Matrix4());

}

MyScene::~MyScene()
{
	if (m_FluidSimulation)
	{
		delete m_FluidSimulation;
		m_FluidSimulation = NULL;
	}
}

bool MyScene::InitialiseGL()
{
	m_Camera->SetPosition(Vector3(-6.25f, 2.0f, 10.0f));

	PhysicsEngine::Instance()->SetGravity(Vector3(0.0f, 0.0f, 0.0f));		//No Gravity
	PhysicsEngine::Instance()->SetDampingFactor(1.0f);						//No Damping

	//Create Ground
	SimpleMeshObject* ground = new SimpleMeshObject("Ground");
	ground->SetMesh(CommonMeshes::Cube(), false);
	ground->SetLocalTransform(Matrix4::Scale(Vector3(20.0f, 0.1f, 20.f)));
	ground->SetColour(Vector4(0.2f, 1.0f, 0.5f, 1.0f));
	ground->SetBoundingRadius(80.0f * 80.f);

	ground->Physics()->SetPosition(Vector3(-6.25f, -0.2f, 0.0f));

	this->AddGameObject(ground);

	BoundingBox2D region;
	region.offset_min = Vector2(0.0f, 5.0f);
	region.scale = Vector2(10.0f, 10.0f);

	const uint div = 12;
	m_FluidSimulation = new FluidSim2D(region, div, div);

	//Initialize Fluid
	uint mid = div / 2;
	uint quater = div / 4;
	mid -= quater / 2;
	m_FluidSimulation->m_MACGrid.SetCellPressure(1.0f, mid, mid, mid + quater, mid + quater);


	m_FluidRenderer = new FluidRenderer(&m_FluidSimulation->m_MACGrid, 128, 128);
	this->AddGameObject(m_FluidRenderer);

	return true;
}

void MyScene::UpdateScene(float msec)
{
	Scene::UpdateScene(msec);

	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_T))
		m_FluidSimulation->StepSimulation(1.0f / 60.0f);

	if (Window::GetKeyboard()->KeyDown(KEYBOARD_R))
		m_FluidSimulation->StepSimulation(1.0f / 60.0f);
}

void MyScene::RenderScene()
{
	Scene::RenderScene();


}

