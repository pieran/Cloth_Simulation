
#include "MyScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\CuboidCollisionShape.h>
#include <ncltech\SphereCollisionShape.h>
#include <ncltech\NCLDebug.h>

MyScene::MyScene(Window& window) : Scene(window)
{
	if (init == true)
		init = InitialiseGL();

	UpdateWorldMatrices(m_RootGameObject, Matrix4());


	PhysicsEngine::Instance()->SetPaused(false); //Start physics immediately
}

MyScene::~MyScene()
{

}

bool MyScene::InitialiseGL()
{
	m_Camera->SetPosition(Vector3(-3.0f, 10.0f, 15.0f));
	m_Camera->SetYaw(-10.f);
	m_Camera->SetPitch(-30.f);

	//Create Some Objects
	{
		SimpleMeshObject* obj = new SimpleMeshObject("CubeA");
		obj->SetLocalTransform(Matrix4::Scale(Vector3(2.0f, 2.0f, 2.f)));
		obj->SetMesh(CommonMeshes::Cube(), false);
		obj->SetColour(Vector4(0.8f, 0.3f, 0.1f, 1.0f));
		obj->Physics()->SetPosition(Vector3(-5.0f, 2.f, -5.0f));
		obj->Physics()->SetCollisionShape(new CuboidCollisionShape(Vector3(2.0f, 2.f, 2.f)));
		obj->Physics()->SetOnCollisionCallback([](PhysicsObject* collidingObject){
			return false;
		});

		this->AddGameObject(obj);
	}

	{
		SimpleMeshObject* obj = new SimpleMeshObject("CubeB");
		obj->SetLocalTransform(Matrix4::Scale(Vector3(2.0f, 1.0f, 2.f)));
		obj->SetMesh(CommonMeshes::Cube(), false);
		obj->SetColour(Vector4(0.5f, 1.0f, 0.5f, 1.0f));
		obj->Physics()->SetPosition(Vector3(0.0f, 0.f, 0.0f));	
		obj->Physics()->SetCollisionShape(new CuboidCollisionShape(Vector3(2.0f, 0.5f, 2.f)));
		obj->Physics()->SetOnCollisionCallback([](PhysicsObject* collidingObject){
			return false;
		});

		_cube = obj;

		this->AddGameObject(obj);
	}


	{
		SimpleMeshObject* obj = new SimpleMeshObject("SphereA");
		obj->SetLocalTransform(Matrix4::Scale(Vector3(2.0f, 2.0f, 2.f)));
		obj->SetMesh(CommonMeshes::Sphere(), false);
		obj->SetColour(Vector4(0.8f, 0.3f, 0.1f, 1.0f));
		obj->Physics()->SetPosition(Vector3(5.0f, 0.5f, -5.0f));
		obj->Physics()->SetCollisionShape(new SphereCollisionShape(2.f));
		obj->Physics()->SetOnCollisionCallback([](PhysicsObject* collidingObject){
			return false;
		});

		this->AddGameObject(obj);
	}

	return true;
}


void MyScene::UpdateScene(float dt)
{
	Scene::UpdateScene(dt);

	const float speed = 1.0f * dt;
	Vector3 movement = Vector3(0, 0, 0);

	if (Window::GetKeyboard()->KeyDown(KEYBOARD_J)) { movement.x -= speed; }
	if (Window::GetKeyboard()->KeyDown(KEYBOARD_L)) { movement.x += speed; }
	if (Window::GetKeyboard()->KeyDown(KEYBOARD_I)) { movement.z -= speed; }
	if (Window::GetKeyboard()->KeyDown(KEYBOARD_K)) { movement.z += speed; }	
	if (Window::GetKeyboard()->KeyDown(KEYBOARD_H)) { movement.y -= speed; }
	if (Window::GetKeyboard()->KeyDown(KEYBOARD_Y)) { movement.y += speed; }

	_cube->Physics()->SetPosition(_cube->Physics()->GetPosition() + movement);
}