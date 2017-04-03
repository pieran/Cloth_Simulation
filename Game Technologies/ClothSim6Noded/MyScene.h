
#pragma once

#define USE_MAN_TANGENTS 0

#include <nclgl\Mesh.h>
#include <ncltech\Scene.h>
#include <ncltech\SimpleMeshObject.h>

#include "GraphObject.h"
#include "VideoEncoder.h"

#if USE_MAN_TANGENTS
#include "Sim_6Noded_ManTangents.h"
typedef Sim_6Noded_ManTangents FESimulation;
#else
#include "FEM6Noded.h"
typedef FEM6Noded FESimulation ;
#endif

#define MAX_HUD_VISIBILITY_TYPES 4


class MyScene : public Scene
{
public:
	MyScene(Window& window);
	~MyScene();

	void InitialiseCloth(uint visual_subdivisions);

	bool InitialiseGL()			override;
	void RenderScene()			override;
	void UpdateScene(float dt)  override;

	void OnResize()				override;


	void SaveCameraData();
	void LoadCameraData();


	GameObject* BuildSphereObject(const std::string& name, const Vector3& pos, float radius, float invmass, const Vector4& color = Vector4(1.0f, 1.0f, 1.0f, 1.0f));
	GameObject* BuildCuboidObject(const std::string& name, const Vector3& pos, const Vector3& scale, float invmass, const Vector4& color = Vector4(1.0f, 1.0f, 1.0f, 1.0f));

	FESimulation*		m_Sim;
	float m_SimTimestep;
	bool m_SimPaused;
	bool m_Drag2D;

	int				m_EncoderVideoNum;
	VideoEncoder*	m_Encoder;

	float           m_ClothUpdateMs;
	GraphObject*	m_GraphObject;
	GraphObject*	m_GraphObjectSolver;

	uint m_HudVisibility;
};