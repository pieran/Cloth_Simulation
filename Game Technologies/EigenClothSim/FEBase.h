#pragma once
#include "ClothBase.h"
#include "mpcg.h"

class FEBase : public ClothBase
{
public:
	FEBase();
	virtual ~FEBase();

	virtual void InitializeClothDesign(const ClothDesignEntity<3>* design) override;
	virtual uint AddPoint(const uint& triIdx, const Vector3& barycentric_coords) override;
	virtual void UpdateTriangleMaterial(const uint& triIdx) override;

	void SetSolverConstraint(int idx, const Matrix3& constraint);

protected:
	virtual void GenRotations(const Vector3* in_clothpos) = 0;
	virtual void GenStressStrainForces(float timestep, const Vector3* in_clothpos, const Vector3* in_clothvel) = 0;
	
	void ComputeVelocity(float subtimestep, const Vector3* in_clothpos, const Vector3* in_clothvel, Vector3* out_clothvel) override;

	void ConstructGlobalMatrix(float subtimestep, const Vector3* in_clothpos, const Vector3* in_phyxelVel);
	void BuildPhyxelTriMappings();
	void UpdateStaticPhyxelConstraints();
	void ReAllocateSolverMemory();

public:
	//POINT-MASS STRUCTURE
	//------------------------
	PArray<Vector3> m_PhyxelPosInitial;				//Initial position

	//SOLVER
	//------------------------
	uint m_SolverAllocated;
	MPCG<SparseRowMatrix<Matrix3>> m_Solver;

	//ADDITIONAL LOCK-FREE INTERMEDIATE STRUCTURES
	//------------------------
	uint											m_NumTriMappings;
	PArray<Vector3>									m_TriB;
	PArray<Matrix3>									m_TriAii;
	PArray<Matrix3>									m_TriAij;

	std::vector<std::vector<uint>>					m_PhyxelTriLookup;
	std::vector<std::map<uint, std::vector<uint>>>	m_PhyxelTriLookupIJ;
};