
#include "MyScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\DistanceConstraint.h>

MyScene::MyScene(Window& window) : Scene(window)
{
	if (init == true)
		init = InitialiseGL();

	UpdateWorldMatrices(m_RootGameObject, Matrix4());
}

MyScene::~MyScene()
{

}

bool MyScene::InitialiseGL()
{
	m_Camera->SetPosition(Vector3(-3.0f, 10.0f, 10.0f));
	m_Camera->SetPitch(-20.f);


	//Create Ground
	GameObject* ground = BuildCuboidObject("Ground", Vector3(0.0f, 0.0f, 0.0f), Vector3(20.0f, 1.0f, 20.0f), 0.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f));
	this->AddGameObject(ground);


	//Create Rope
	GameObject* prevRopeNode = NULL;
	for (int i = 0; i < 5; i++)
	{
		GameObject* ropeNode = BuildSphereObject("", 
			Vector3(i * 1.1f, 7.f, 0.0f), 
			0.5f, 
			(prevRopeNode == NULL) ? 0.0f : 1.0f,
			Vector4(1.0f, 0.2f, 0.5f, 1.0f));

		if (prevRopeNode != NULL)
		{
			DistanceConstraint* constraint = new DistanceConstraint(prevRopeNode->Physics(), ropeNode->Physics(), Vector3((i - 1) * 1.1f , 7.f, 0.0f), Vector3(i * 1.1f, 7.f, 0.0f));
			PhysicsEngine::Instance()->AddConstraint(constraint);
		}

		this->AddGameObject(ropeNode);
		prevRopeNode = ropeNode;
	}


	//Create Bridge
	GameObject* prevBridgeNode = NULL;
	Vector3 bridgeStart = Vector3(-10.f, 5.f, -5.f);
	for (int i = 0; i < 10; i++)
	{
		GameObject* bridgeNode = BuildCuboidObject("",
			Vector3(i * 1.1f, 0.0f, 0.0f) + bridgeStart,
			Vector3(0.5f, 0.1f, 1.0f),
			(i == 0 || i == 9) ? 0.0f : 1.0f,
			Vector4(1.0f, 0.2f, 0.5f, 1.0f));

		if (prevBridgeNode != NULL)
		{
			DistanceConstraint* constraint = new DistanceConstraint(prevBridgeNode->Physics(), bridgeNode->Physics(), Vector3((i - 1) * 1.1f + 0.5f, 0.0f, 0.0f) + bridgeStart, Vector3(i * 1.1f - 0.5f, 0.0f, 0.0f) + bridgeStart);
			PhysicsEngine::Instance()->AddConstraint(constraint);
		}

		this->AddGameObject(bridgeNode);
		prevBridgeNode = bridgeNode;
	}

	return true;
}


GameObject* MyScene::BuildSphereObject(const std::string& name, const Vector3& pos, float radius, float invmass, const Vector4& color)
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

GameObject* MyScene::BuildCuboidObject(const std::string& name, const Vector3& pos, const Vector3& halfdims, float invmass, const Vector4& color)
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
