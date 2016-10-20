
#pragma once

#include <nclgl\Mesh.h>
#include <ncltech\Scene.h>
#include <ncltech\SimpleMeshObject.h>
#include "ClothDesignEntity.h"
#include "ClothRenderObject.h"

class ClothScene : public Scene
{
public:
	ClothScene(Window& window);
	~ClothScene();

	void InitialiseCloth();

	bool InitialiseGL()	override;
	void RenderScene() override;
	void UpdateScene(float dt)  override;

protected:
	GameObject* BuildSphereObject(const std::string& name, const Vector3& pos, float radius, float invmass, const Vector4& color = Vector4(1.0f, 1.0f, 1.0f, 1.0f));
	GameObject* BuildCuboidObject(const std::string& name, const Vector3& pos, const Vector3& scale, float invmass, const Vector4& color = Vector4(1.0f, 1.0f, 1.0f, 1.0f));

	PositionBasedDynamicsMS* m_Cloth;
	GameObject* m_Corners[4];

	ClothRenderObject* m_ClothRenderer;
};