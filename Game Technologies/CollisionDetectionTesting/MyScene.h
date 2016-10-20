
#pragma once

#include <nclgl\Mesh.h>
#include <ncltech\Scene.h>

class MyScene : public Scene
{
public:
	MyScene(Window& window);
	~MyScene();

	bool InitialiseGL()	override;
	void UpdateScene(float dt)  override;

protected:
	GameObject* _cube;
};