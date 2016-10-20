
#include "ClothScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\DistanceConstraint.h>

#include "PositionBasedDynamics.h"

ClothScene::ClothScene(Window& window) : Scene(window), m_ClothRenderer(NULL)
{
	if (init == true)
		init = InitialiseGL();

	UpdateWorldMatrices(m_RootGameObject, Matrix4());
}

ClothScene::~ClothScene()
{
}
const uint n_divisions = 65;
void ClothScene::InitialiseCloth()
{

	for (int i = 0; i < 4; ++i)
	{
		m_Corners[i] = BuildSphereObject("", Vector3(0.0f, 0.0f, 0.0f), 0.02f, 0.0f, Vector4(1.0f, 0.0f, 0.0f, 1.0f));
		this->AddGameObject(m_Corners[i]);
	}

	Ellipsoid e;
	e.Transform = Matrix4::Translation(Vector3(0.0f, 0.5f, 0.0f));
	e.radius = 0.2f;
	e.invTransform = Matrix4::Inverse(e.Transform);


	GameObject* e_obj = BuildSphereObject("", Vector3(0.0f, 0.5f, 0.0f), 0.2f, 2.0f);
	this->AddGameObject(e_obj);

	Matrix4 transform = Matrix4::Translation(Vector3(0, 1.0f, 0.0f)) * Matrix4::Rotation(90.0f, Vector3(1, 0, 0));



	ClothDesignEntity cloth_design;
	cloth_design.GenerateDefaultClothVertices(n_divisions, 1.0f, 1.0f, transform); //1m x 1m (subdivided 33 times along x/y
	cloth_design.GenerateDefaultClothTriangles(n_divisions);

	//Set Static Vertices
	/*for (uint i = 0; i < n_divisions - 1; ++i)
	{
		cloth_design.Vertices()[i].flags = VERTEXFLAGS_IS_STATIC;
	}*/

	//cloth_design.Vertices()[(n_divisions - 1) * n_divisions + 0].flags = VERTEXFLAGS_IS_STATIC;
	//cloth_design.Vertices()[(n_divisions - 1) * n_divisions + n_divisions - 1].flags = VERTEXFLAGS_IS_STATIC;

	cloth_design.Vertices()[0].force												 = Vector3(0, -1.0f, 0);
	cloth_design.Vertices()[n_divisions - 1].force									 = Vector3(0, -1.0f, 0);
	cloth_design.Vertices()[(n_divisions - 1) * n_divisions + 0].force				 = Vector3(0, -1.0f, 0);
	cloth_design.Vertices()[(n_divisions - 1) * n_divisions + n_divisions - 1].force = Vector3(0, -1.0f, 0);

	m_Cloth = new PositionBasedDynamicsMS();
	m_Cloth->simulation_OnClothDesignChanged(&cloth_design);
	m_Cloth->AddCollisionEllipsoid(e);

	if (!m_ClothRenderer)
	{
		m_ClothRenderer = new ClothRenderObject(m_Cloth);
		this->AddGameObject(m_ClothRenderer);
	}
}

bool ClothScene::InitialiseGL()
{
	InitialiseCloth();

	m_Camera->SetPosition(Vector3(-3.0f, 10.0f, 10.0f));
	m_Camera->SetPitch(-20.f);


	//Create Ground
	GameObject* ground = BuildCuboidObject("Ground", Vector3(0.0f, -1.001f, 0.0f), Vector3(20.0f, 1.0f, 20.0f), 0.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f));
	this->AddGameObject(ground);

	return true;
}

void ClothScene::RenderScene()
{
	Scene::RenderScene();
}

void ClothScene::UpdateScene(float dt)
{
	Scene::UpdateScene(dt);

	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_P))
		m_ClothRenderer->TogglePaused();

	m_Corners[0]->Physics()->SetPosition(m_Cloth->m_PhyxelsPosNew[0]);
	m_Corners[1]->Physics()->SetPosition(m_Cloth->m_PhyxelsPosNew[n_divisions - 1]);
	m_Corners[2]->Physics()->SetPosition(m_Cloth->m_PhyxelsPosNew[(n_divisions - 1) * n_divisions + 0]);
	m_Corners[3]->Physics()->SetPosition(m_Cloth->m_PhyxelsPosNew[(n_divisions - 1) * n_divisions + n_divisions - 1]);
}

GameObject* ClothScene::BuildSphereObject(const std::string& name, const Vector3& pos, float radius, float invmass, const Vector4& color)
{
	//Create Ground
	SimpleMeshObject* sphere = new SimpleMeshObject(name);
	sphere->SetMesh(CommonMeshes::Sphere(), false);
	sphere->SetLocalTransform(Matrix4::Scale(Vector3(radius, radius, radius)));
	sphere->SetColour(color);
	sphere->SetBoundingRadius(radius);

	sphere->Physics()->SetPosition(pos);
	sphere->Physics()->SetInverseMass(invmass);

	//Build Inverse Inertia Matrix of a sphere
	Matrix3 inertia;
	{
		float i = 2.5f * invmass * radius * radius;

		inertia._11 = i;
		inertia._22 = i;
		inertia._33 = i;
	}
	sphere->Physics()->SetInverseInertia(inertia);

	return sphere;
}

GameObject* ClothScene::BuildCuboidObject(const std::string& name, const Vector3& pos, const Vector3& halfdims, float invmass, const Vector4& color)
{
	SimpleMeshObject* cuboid = new SimpleMeshObject(name);
	cuboid->SetMesh(CommonMeshes::Cube(), false);
	cuboid->SetLocalTransform(Matrix4::Scale(halfdims));
	cuboid->SetColour(color);
	cuboid->SetBoundingRadius(80.0f * 80.f);

	cuboid->Physics()->SetPosition(pos);

	cuboid->Physics()->SetInverseMass(invmass);
	//Build Inverse Inertia Matrix of a cuboid
	Matrix3 inertia;
	{
		Vector3 dimsSq = (halfdims + halfdims);
		dimsSq = dimsSq * dimsSq;

		inertia._11 = 12.f * invmass * 1.f / (dimsSq.y + dimsSq.z);
		inertia._22 = 12.f * invmass * 1.f / (dimsSq.x + dimsSq.z);
		inertia._33 = 12.f * invmass * 1.f / (dimsSq.x + dimsSq.y);
	}
	cuboid->Physics()->SetInverseInertia(inertia);

	return cuboid;
}
