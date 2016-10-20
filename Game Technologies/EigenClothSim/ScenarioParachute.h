#pragma once

#include "ClothScenario.h"
#include "ClothDesignEntity.h"
#include <nclgl\Vector3.h>
#include <nclgl\Vector2.h>

class ScenarioParachute : public ClothScenario
{
public:
	ScenarioParachute(ClothBase* clothsim, int numDivsX = 10, int numDivsY = 10, const Vector2& scale = Vector2(4.0f, 4.0f), const Vector3& pos = Vector3(0, 1.0f, -5.0f))
		: ClothScenario(clothsim), m_ClothDivsX(numDivsX), m_ClothDivsY(numDivsY), m_ClothScale(scale), m_ClothPosition(pos), m_AccumTime(0.0f)
	{
		using namespace std::placeholders;
		clothsim->SetCallback_OnInnerFullTimestep(std::bind(&ScenarioParachute::OnInnerFullTimestep, this, _1));
		clothsim->SetCallback_OnInnerApplySubTimestep(std::bind(&ScenarioParachute::OnInnerApplySubTimestep, this, _1, _2, _3, _4));
		m_AccumTime = 0.0f;
	}

	virtual ~ScenarioParachute() {}

	//Cleanup any non-cloth data structures
	virtual void OnDeleteScene() override
	{

	}

	//Set All data in cloth class and scenario class
	virtual void OnInitializeScene() override
	{
		m_AccumTime = 0.0f;


		Matrix4 transform = Matrix4::Translation(m_ClothPosition)*Matrix4::Rotation(90.0f, Vector3(1, 0, 0));
		Mat44 tmp;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				tmp(i, j) = transform[j * 4 + i];
			}
		}


		//Create vertical cloth with fixed points at the top and bottom rows
		ClothDesignEntity<3> design;
		design.GenerateDefaultClothVertices(m_ClothDivsX, m_ClothScale.x, m_ClothScale.y, &tmp); 
		design.GenerateDefaultClothTriangles(m_ClothDivsX);

		//Set Static Vertices
		for (Vertex& v : design.Vertices())
		{
			v.flags = 0;
		}
		design.Vertices()[0].flags = VERTEXFLAGS_IS_STATIC;
		design.Vertices()[m_ClothDivsX - 1].flags = VERTEXFLAGS_IS_STATIC;

		design.Vertices()[(m_ClothDivsY - 1) * m_ClothDivsX + 0].flags = VERTEXFLAGS_IS_STATIC;
		design.Vertices()[(m_ClothDivsY - 1) * m_ClothDivsX + m_ClothDivsX - 1].flags = VERTEXFLAGS_IS_STATIC;

		m_ClothSim->InitializeClothDesign(&design);
	}

	//Called once per cloth update
	void OnOuterUpdateScene(float cloth_timestep)
	{

	}

	//Called once per timestep (possibly multiple cloth steps to each scene timestep) NOTE: Sum of sub_timesteps each frame will equal cloth_timestep so no need to update on both callbacks!
	void OnInnerFullTimestep(float timestep)
	{
		m_AccumTime += timestep;
	}

	void OnInnerApplySubTimestep(float sub_timestep, uint numPhyxels, Vector3* clothPos, Vector3* clothVel)
	{

		for (uint i = 0; i < numPhyxels; ++i)
		{
			if (m_ClothSim->IsPhyxelStatic(i))
			{
				float spd = 3.85f;

				clothVel[i].y = sin((m_AccumTime+ sub_timestep) * 4.0f) * spd;
			}
		}
	}

	//Called once per timestep (possibly multiple cloth steps to each scene timestep)
	void OnOuterCollisions()
	{

	}

	//Called once per cloth update
	void OnInnerCollisions()
	{

	}

protected:
	int		m_ClothDivsX, m_ClothDivsY;
	Vector2 m_ClothScale;
	Vector3 m_ClothPosition;
	float	m_AccumTime;
};
