
#include "MyScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\CuboidCollisionShape.h>
#include <ncltech\SphereCollisionShape.h>
#include <ncltech\NCLDebug.h>

#include "Player.h"

MyScene::MyScene(Window& window) : Scene(window)
{
	glGenTextures(1, &m_whiteTexture);
	glBindTexture(GL_TEXTURE_2D, m_whiteTexture);
	int white_pixel = 0xFFFFFFFF;
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, &white_pixel);


	if (init == true)
		init = InitialiseGL();

	UpdateWorldMatrices(m_RootGameObject, Matrix4());


	PhysicsEngine::Instance()->SetPaused(false);
}

MyScene::~MyScene()
{
	if (m_whiteTexture)
	{
		glDeleteTextures(1, &m_whiteTexture);
		m_whiteTexture = NULL;
	}
}

bool MyScene::InitialiseGL()
{
	m_Camera->SetPosition(Vector3(-3.0f, 10.0f, 15.0f));
	m_Camera->SetYaw(-10.f);
	m_Camera->SetPitch(-30.f);

	//Create Ground
	this->AddGameObject(BuildCuboidObject("Ground", Vector3(0.0f, -1.001f, 0.0f), Vector3(20.0f, 1.0f, 20.0f), 0.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f)));

	//Create Player
	this->AddGameObject(new Player("Player1"));


	//Create Some Objects
	{
		SimpleMeshObject* obj = new SimpleMeshObject("House");
		obj->SetLocalTransform(Matrix4::Scale(Vector3(2.0f, 2.0f, 2.f)));
		obj->SetMesh(new OBJMesh(MESHDIR"house.obj"), true);
		obj->SetTexture(m_whiteTexture, false);
		obj->SetColour(Vector4(0.8f, 0.3f, 0.1f, 1.0f));
		obj->Physics()->SetPosition(Vector3(-5.0f, 2.f, -5.0f));
		obj->Physics()->SetCollisionShape(new CuboidCollisionShape(Vector3(2.0f, 2.f, 2.f)));
		obj->Physics()->SetOnCollisionCallback([](PhysicsObject* collidingObject){
			NCLDebug::Log(Vector3(0.6f, 0.3f, 0.1f), "You are inside the house!");
			return false;
		});

		this->AddGameObject(obj);
	}

	{
		SimpleMeshObject* obj = new SimpleMeshObject("Garden");
		obj->SetLocalTransform(Matrix4::Scale(Vector3(2.0f, 1.0f, 2.f)));
		obj->SetMesh(new OBJMesh(MESHDIR"garden.obj"), true);
		obj->SetTexture(m_whiteTexture, false);
		obj->SetColour(Vector4(0.5f, 1.0f, 0.5f, 1.0f));
		obj->Physics()->SetPosition(Vector3(5.0f, 0.5f, -5.0f));
		obj->Physics()->SetCollisionShape(new CuboidCollisionShape(Vector3(2.0f, 0.5f, 2.f)));
		obj->Physics()->SetOnCollisionCallback([](PhysicsObject* collidingObject){
			NCLDebug::Log(Vector3(0.0f, 1.0f, 0.0f), "You are inside the garden!");
			return false;
		});
		
		this->AddGameObject(obj);
	}

	{
		PhysicsObject* obj = new PhysicsObject();
		obj->SetPosition(Vector3(5.0f, 1.0f, 0.0f));
		obj->SetCollisionShape(new SphereCollisionShape(1.0f));
		obj->SetOnCollisionCallback(
			[obj](PhysicsObject* collidingObject) {

			NCLDebug::Log(Vector3(1.0f, 0.0f, 0.0f), "You found the secret!");

			float r_x = 5.f * ((rand() % 200) / 100.f - 1.0f);
			float r_z = 3.f * ((rand() % 200) / 100.f - 1.0f);
			obj->SetPosition(Vector3(r_x, 1.0f, r_z + 3.0f));
			return false;
		});
		PhysicsEngine::Instance()->AddPhysicsObject(obj);
	}

	return true;
}

GameObject* MyScene::BuildCuboidObject(const std::string& name, const Vector3& pos, const Vector3& halfdims, float invmass, const Vector4& color)
{
	SimpleMeshObject* cuboid = new SimpleMeshObject(name);
	cuboid->SetMesh(CommonMeshes::Cube(), false);
	cuboid->SetLocalTransform(Matrix4::Scale(halfdims));
	cuboid->SetColour(color);
	cuboid->SetBoundingRadius(Vector3::Dot(halfdims, halfdims));

	cuboid->Physics()->SetPosition(pos);

	return cuboid;
}
