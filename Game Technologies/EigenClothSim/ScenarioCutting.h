#pragma once

#include "ClothScenario.h"
#include "ClothDesignEntity.h"
#include <ncltech\NCLDebug.h>
#include <nclgl\Vector3.h>
#include <nclgl\Vector2.h>

class ScenarioCutting : public ClothScenario
{
public:
	ScenarioCutting(ClothBase* clothsim, int numDivsX = 10, int numDivsY = 10, const Vector2& scale = Vector2(4.0f, 4.0f), const Vector3& pos = Vector3(0, 1.0f, -5.0f))
		: ClothScenario(clothsim), m_ClothDivsX(numDivsX), m_ClothDivsY(numDivsY), m_ClothScale(scale), m_ClothPosition(pos), m_AccumTime(0.0f)
	{
		using namespace std::placeholders;
		clothsim->SetCallback_OnInnerFullTimestep(std::bind(&ScenarioCutting::OnInnerFullTimestep, this, _1));
		clothsim->SetCallback_OnInnerApplySubTimestep(std::bind(&ScenarioCutting::OnInnerApplySubTimestep, this, _1, _2, _3, _4));
		m_AccumTime = 0.0f;
	}

	virtual ~ScenarioCutting() {}

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

				clothVel[i].y = sin((m_AccumTime + sub_timestep) * 4.0f) * spd;
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

	//Called each time cloth is 'Cut'. Giving the start and end particle indicies
	virtual void OnClothCut(int start_idx, int end_idx, const Ray& start_ray, const Ray& end_ray, const Plane& cut_plane) override
	{
		NCLDebug::Log(Vector3(1.0f, 0.0f, 0.0f), "%d: Cloth Cut!", time(NULL));

		//Find Triangle Idx
		uint numTris = m_ClothSim->m_NumTris * 3;
		for (int i = 0; i < numTris; i+=3)
		{
			uint indices[3] = { m_ClothSim->m_TriIndicies[i]
								, m_ClothSim->m_TriIndicies[i + 1]
								, m_ClothSim->m_TriIndicies[i + 2] };
			Vector3 positions[3] = { m_ClothSim->m_ClothStatesPos[indices[0]]
								, m_ClothSim->m_ClothStatesPos[indices[1]]
								, m_ClothSim->m_ClothStatesPos[indices[2]] };

			//Iterate over each edge
			bool intersects = false;
			Vector3 outP;
			for (int j = 0; intersects == false && j < 3; ++j)
			{
				int k = (j + 1) % 3; //A-B, B-C, C-A

				uint idxA = indices[j];
				uint idxB = indices[k];

				//Check if the edge intersects the plane
				intersects = cut_plane.GetSegmentPlaneIntersection(positions[j], positions[k], outP);
			}

			if (intersects)
			{
				NCLDebug::DrawThickLine(positions[0], positions[1], 0.1f);
				NCLDebug::DrawThickLine(positions[1], positions[2], 0.1f);
				NCLDebug::DrawThickLine(positions[2], positions[0], 0.1f);

				//Sub-Divide Triangle Quickly/Hackily
				

				//1. Add new Point at one of the edges
				uint idxD = m_ClothSim->AddPoint(i, Vector3(0.5f, 0.5f, 0.0f)); //Half way between triangle points A-B

				//2. Update Original Triangle to be A-D-C (not A-B-C)
				m_ClothSim->m_TriIndicies[i + 1] = idxD;
				m_ClothSim->UpdateTriangleMaterial(i);

				//3. Add new triangle to cover the other half (B-D-C)
				uint idx = m_ClothSim->AddTriangle(indices[1], idxD, indices[2]);
				m_ClothSim->UpdateTriangleMaterial(idx);
			}

		}

	}

protected:
	int		m_ClothDivsX, m_ClothDivsY;
	Vector2 m_ClothScale;
	Vector3 m_ClothPosition;
	float	m_AccumTime;
};
