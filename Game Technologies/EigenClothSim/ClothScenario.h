#pragma once

#include "ClothBase.h"
#include <ncltech\GameObject.h>
#include <nclgl/Plane.h>

class ClothScenario : public GameObject
{
	friend class ClothScenarioViewer;
public:
	ClothScenario(ClothBase* clothsim) : m_ClothSim(clothsim) {}

	virtual ~ClothScenario() {}

	virtual void OnDeleteScene() = 0;			//Cleanup any non-cloth data structures
	virtual void OnInitializeScene() = 0;		//Set All data in cloth class and scenario class NOTE: Called by the cloth class as a callback upon being connected to the scenario
	virtual void OnClothCut(int start_idx, int end_idx, const Ray& start_ray, const Ray& end_ray, const Plane& cut_plane) {}

	void OnResetScene()
	{
		OnDeleteScene();
		OnInitializeScene();
	}

	void OnUpdateScene(float full_timestep)
	{
		m_ClothSim->UpdateSimulation(full_timestep);
	}

	ClothBase* GetCloth() { return m_ClothSim; }

protected:
	virtual void OnRenderObject() { };				//Handles OpenGL calls to Render the object
	virtual void OnUpdateObject(float dt) {};		//NOT USED!!!!
	virtual void OnAttachedToScene() {};			//NOT USED!!!!

	/* !!!!!POSSIBLE CALLBACKS!!!!!!
	virtual void OnOuterUpdateScene(float cloth_timestep) = 0;	//Called once per cloth update
	virtual void OnInnerUpdateScene(float sub_timestep) = 0;	//Called once per timestep (possibly multiple cloth steps to each scene timestep) NOTE: Sum of sub_timesteps each frame will equal cloth_timestep so no need to update on both callbacks!

	virtual void OnOuterCollisions() = 0;				//Called once per timestep (possibly multiple cloth steps to each scene timestep)
	virtual void OnInnerCollisions() = 0;				//Called once per cloth update
	*/

protected:
	ClothBase* m_ClothSim;
};