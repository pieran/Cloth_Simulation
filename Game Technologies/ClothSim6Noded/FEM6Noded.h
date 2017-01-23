#pragma once


#include "mpcg.h"

#include "EigenDefines.h"
#include "SimulationDefines.h"

#include <nclgl\Matrix4.h>
#include <nclgl\Matrix3.h>
#include <nclgl\Vector3.h>
#include <nclgl\Vector2.h>
#include "ProfilingTimer.h"


#define USE_DYNAMIC_MINMAX FALSE



typedef Eigen::Matrix<float, 12, 12> StiffnessMatrix;

typedef Eigen::Matrix<float, 3, 18> BMatrix;
typedef Eigen::Matrix<float, 2, 2>  JaMatrix;
typedef Eigen::Matrix<float, 6, 18> GMatrix;
typedef Eigen::Matrix<float, 2, 6> Mat26;
typedef Eigen::Matrix<float, 18, 1> VDisplacement;


class FEM6Noded
{
	friend class ClothRenderObject;
	friend class MyScene;
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	FEM6Noded();
	~FEM6Noded();

	//Simulation
	void simulation_OnClothDesignChanged(uint nTris, FETriangle* tris, uint nVerts, FEVertDescriptor* verts);
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
	float GenFEMTriangle(unsigned int idx, FETriangle& tri);
	void BuildStiffnessMatrix(unsigned int idx, FETriangle& tri);

	Matrix3 BuildRotationMatrix(const Vector3& a, const Vector3& b, const Vector3& c);


	void BuildRotationAndNormals();
	void SimpleCorotatedBuildAMatrix(float dt);

	void InitGaussWeights();

	void BuildTransformationMatrix(const FETriangle& tri, const std::vector<Vector3>& pos, const Vec3& gaussPoint, const float warp_angle, Mat33& out_t, JaMatrix& out_ja, Mat26& out_dn);

	void CalcBMatrix(const FETriangle& tri, const std::vector<Vector3>& pos, const Vec3& gaussPoint, const float warp_angle, const VDisplacement& displacements, BMatrix& out_b, JaMatrix& out_ja, GMatrix& out_g);

protected:

	Eigen::Matrix<float, 12, 3> GaussPoint12;
	Eigen::Matrix<float, 12, 1> GaussWeight12;


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
	std::vector<FETriangle> m_Triangles;
	std::vector<StiffnessMatrix, Eigen::aligned_allocator<StiffnessMatrix>> m_TriangleStiffness;


	MPCG<SparseRowMatrix<Matrix3>> m_Solver;	//Solver

	float angle = 0.0f;



};
