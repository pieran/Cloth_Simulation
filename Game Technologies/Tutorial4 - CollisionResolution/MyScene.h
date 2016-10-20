
#pragma once

#include <nclgl\Mesh.h>
#include <ncltech\Scene.h>
#include <ncltech\SimpleMeshObject.h>

class MyScene : public Scene
{
public:
	MyScene(Window& window);
	~MyScene();

	bool InitialiseGL()	override;

protected:
	GameObject* BuildSphereObject(const std::string& name, const Vector3& pos, float radius, float invmass, const Vector4& color = Vector4(1.0f, 1.0f, 1.0f, 1.0f));
	GameObject* BuildCuboidObject(const std::string& name, const Vector3& pos, const Vector3& scale, float invmass, const Vector4& color = Vector4(1.0f, 1.0f, 1.0f, 1.0f));

};