
#include "MyScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\CuboidCollisionShape.h>
#include <ncltech\SphereCollisionShape.h>
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
	m_Camera->SetPosition(Vector3(-3.0f, 10.0f, 15.0f));
	m_Camera->SetYaw(-10.f);
	m_Camera->SetPitch(-30.f);

	//Create Ground
	this->AddGameObject(BuildCuboidObject("Ground", Vector3(0.0f, -1.0f, 0.0f), Vector3(20.0f, 1.0f, 20.0f), 0.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f)));

	//Create Bouncing Boxes
	this->AddGameObject(BuildCuboidObject("bb_lower", Vector3(5.0f, 2.5f, 4.0f), Vector3(1.0f, 0.5f, 0.5f), 1.0f, Vector4(1.0f, 0.7f, 0.5f, 1.0f)));
	this->AddGameObject(BuildCuboidObject("bb_upper", Vector3(5.0f, 5.5f, 4.0f), Vector3(2.0f, 0.5f, 0.5f), 1.0f, Vector4(1.0f, 0.7f, 0.5f, 1.0f)));

	//Create Colliding Spheres
	this->AddGameObject(BuildSphereObject("", Vector3(-1.0f, 0.5f, 6.0f), 0.5f, 1.0f, Vector4(0.7f, 1.0f, 0.7f, 1.0f)));
	this->AddGameObject(BuildSphereObject("", Vector3(0.0f, 0.5f, 6.0f), 0.5f, 1.0f, Vector4(0.7f, 1.0f, 0.7f, 1.0f)));
	this->AddGameObject(BuildSphereObject("", Vector3(1.0f, 0.5f, 6.0f), 0.5f, 1.0f, Vector4(0.7f, 1.0f, 0.7f, 1.0f)));

	GameObject* movingSphere = BuildSphereObject("", Vector3(5.0f, 0.5f, 6.0f), 0.5f, 1.0f, Vector4(1.0f, 0.7f, 0.7f, 1.0f));
	movingSphere->Physics()->SetLinearVelocity(Vector3(-5.f, 0.0f, 0.0f));
	movingSphere->Physics()->SetAngularVelocity(Vector3(0.0f, 0.0f, 0.5f));
	movingSphere->Physics()->SetFriction(0.5f);
	this->AddGameObject(movingSphere);

	//Create Ramp
	GameObject* ramp = BuildCuboidObject("Ramp", Vector3(4.0f, 3.5f, -5.0f), Vector3(5.0f, 0.5f, 4.0f), 0.0f, Vector4(1.0f, 0.7f, 1.0f, 1.0f));
	ramp->Physics()->SetOrientation(Quaternion::AxisAngleToQuaterion(Vector3(0.0f, 0.0f, 1.0f), 20.0f));
	ramp->Physics()->SetFriction(1.0f);
	this->AddGameObject(ramp);

	//Create Cubes to roll on ramp
	for (int i = 0; i < 5; ++i)
	{
		Vector4 colour = Vector4(i * 0.25f, 0.7f, (2 - i) * 0.25f, 1.0f);
		GameObject* cube = BuildCuboidObject("", Vector3(8.0f, 6.0f, -7.0f + i * 1.1f), Vector3(0.5f, 0.5f, 0.5f), 1.f, colour);
		cube->Physics()->SetFriction(i * 0.05f);
		cube->Physics()->SetOrientation(Quaternion::AxisAngleToQuaterion(Vector3(0.0f, 0.0f, 1.0f), 200.0f));
		this->AddGameObject(cube);
	}

	//Create Cubes Tower
	for (int i = 0; i < 5; ++i)
	{
		Vector4 colour = Vector4(i * 0.25f, 0.7f, (2 - i) * 0.25f, 1.0f);
		GameObject* cube = BuildCuboidObject("", Vector3(8.0f, 0.5f + i * 1.0f, -12.0f), Vector3(0.5f, 0.5f, 0.5f), 1.f, colour);
		cube->Physics()->SetFriction(1.0f);
		this->AddGameObject(cube);
	}



	//Create Top of Cradle
	GameObject* cradleTop = BuildCuboidObject("", Vector3(0.0f, 12.0f, -12.0f), Vector3(5.f, 0.5f, 0.5f), 0.0f, Vector4(0.2f, 0.2f, 0.2f, 1.0f));
	this->AddGameObject(cradleTop);

	//Create the balls hanging from the cradle
	GameObject* cradleBall;
	for (int i = 0; i < 4; ++i)
	{
		if (i == 0)
			cradleBall = BuildSphereObject("", Vector3(-4.0f, 7.0f, -12.0f), 1.0f, 1.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f));
		else
		//	cradleBall = BuildSphereObject("", Vector3(-4.0f + i * 2.001f, 7.0f, -12.0f), 1.0f, 1.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f));
			cradleBall = BuildCuboidObject("", Vector3(-4.0f + i * 2.001f, 7.0f, -12.0f), Vector3(1.0f, 1.0f, 1.0f), 1.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f));

		PhysicsEngine::Instance()->AddConstraint(new DistanceConstraint(cradleBall->Physics(), cradleTop->Physics(), cradleBall->Physics()->GetPosition(), Vector3(-4.0f + i * 2.001f, 11.5f, -12.f)));
		this->AddGameObject(cradleBall);
	}

	cradleBall = BuildSphereObject("", Vector3(8.503f, 11.5f, -12.0f), 1.0f, 1.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f));
	PhysicsEngine::Instance()->AddConstraint(new DistanceConstraint(cradleBall->Physics(), cradleTop->Physics(), cradleBall->Physics()->GetPosition(), Vector3(4.003, 11.5f, -12.f)));
	this->AddGameObject(cradleBall);


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

	sphere->Physics()->SetCollisionShape(new SphereCollisionShape(radius));

	sphere->Physics()->SetInverseMass(invmass);
	sphere->Physics()->SetInverseInertia(sphere->Physics()->GetCollisionShape()->BuildInverseInertia(sphere->Physics()->GetInverseMass()));

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

	cuboid->Physics()->SetCollisionShape(new CuboidCollisionShape(halfdims));

	cuboid->Physics()->SetInverseMass(invmass);
	cuboid->Physics()->SetInverseInertia(cuboid->Physics()->GetCollisionShape()->BuildInverseInertia(cuboid->Physics()->GetInverseMass()));

	return cuboid;
}
