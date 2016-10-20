#pragma once



#include "ClothDesignEntity.h"
#include "mpcg.h"

#include "EigenDefines.h"
#include "common.h"
#include "FECommon.h"

#include <nclgl\Matrix4.h>
#include <nclgl\Matrix3.h>
#include <nclgl\Vector3.h>
#include <nclgl\Vector2.h>
#include "ProfilingTimer.h"


#define USE_DYNAMIC_MINMAX FALSE

class FEM3Noded
{
	friend class ClothRenderObject;
	friend class MyScene;
public:
	FEM3Noded();
	~FEM3Noded();

	std::vector<Ellipsoid>		m_Ellipsoids;

	//Simulation
	void simulation_OnClothDesignChanged(ClothDesignEntity<3>* design);
	void Simulation_StepSimulation(float dt);

	void UpdateConstraints();
	bool ValidateVelocityTimestep();

	//Misc
	void Render_DrawingToVisualiser();



	//Profiling
	ProfilingTimer    m_ProfilingRotationNormals;
	ProfilingTimer    m_ProfilingBuildMatricies;
	ProfilingTimer	  m_ProfilingConjugateGradient;
	ProfilingTimer	  m_ProfilingExternalCollisions;
	ProfilingTimer	  m_ProfilingTotalTime;

protected:
	void GenFEMTriangle(unsigned int idx, FETriangleLite& tri);
	Matrix3 BuildRotationMatrix(const Vector3& a, const Vector3& b, const Vector3& c);


	void BuildRotationAndNormals();
	void SimpleCorotatedBuildAMatrix(float dt);

	void CollideEllipsoid(const Ellipsoid& e);
protected:



	uint m_NumWidth, m_NumHeight;
	uint m_NumTotal, m_NumTriangles;
	float m_TotalArea;

	float						m_TimeStep;
	uint						m_StepCounter;
	uint						m_StepsBeforeIncrease;
	//Phyxel Data

	std::vector<Vector3>		m_PhyxelsPos;
	std::vector<Vector3>		m_PhyxelsPosTemp;	//Temp Positions to validate timestep
	std::vector<Vector3>		m_PhyxelsPosInitial;
	std::vector<Vector3>		m_PhyxelsVel;
	std::vector<Vector3>		m_PhyxelsVelChange;
	std::vector<Vector3>		m_PhyxelForces;
	std::vector<Vector3>		m_PhyxelExtForces;
	std::vector<Vector3>		m_PhyxelNormals;
	std::vector<Vector2>		m_PhyxelTexCoords;
	std::vector<bool>			m_PhyxelIsStatic;
	std::vector<float>			m_PhyxelsMass;
	std::vector<float>			m_RenderValue;	//Stress/Strain etc for rendering
	Vector2						m_RenderValueMinMax;
	std::vector<uint>			m_RenderIndices;

	//Structural Data
	std::vector<FETriangleLite> m_Triangles;
	std::vector<Quad>			m_Quads;	//For generating smooth normals


	MPCG<SparseRowMatrix<Matrix3>> m_Solver;	//Solver

	float angle = 0.0f;

};
