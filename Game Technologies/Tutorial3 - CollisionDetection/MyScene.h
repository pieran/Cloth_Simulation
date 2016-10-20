
#pragma once

#include <nclgl\Mesh.h>
#include <ncltech\Scene.h>

class MyScene : public Scene
{
public:
	MyScene(Window& window);
	~MyScene();

	bool InitialiseGL()	override;

protected:
	GameObject* BuildCuboidObject(const std::string& name, const Vector3& pos, const Vector3& scale, float invmass, const Vector4& color = Vector4(1.0f, 1.0f, 1.0f, 1.0f));
	
	GLuint m_whiteTexture;

};