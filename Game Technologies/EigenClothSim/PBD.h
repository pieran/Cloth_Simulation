#pragma once

#include "ClothBase.h"

struct PBDDistanceConstraint { uint p1, p2;	float rest_length, k, k_prime; };
struct PBDBendingConstraint { uint p1, p2, p3; float rest_length, w, k, k_prime; };

class PBD : public ClothBase
{
public:
	PBD();
	virtual ~PBD();

	virtual void InitializeClothDesign(const ClothDesignEntity<3>* design) override;
	void InitializePBDConstraints();

	inline uint& SolverIterations() { return m_SolverIterations; }
	inline const uint& SolverIterations() const { return m_SolverIterations; }

protected:
	void ComputeVelocity(float subtimestep, const Vector3* in_clothpos, const Vector3* in_clothvel, Vector3* out_clothvel) override;

	void BuildLookupArrays();

	bool BuildDistanceConstraint(uint idxA, uint idxB, float k, PBDDistanceConstraint& out_constraint);
	bool BuildBendingConstraint(uint idxA, uint idxB, uint idxC, float k, PBDBendingConstraint& out_constraint);

	void SolveDistanceConstraints(float weighting);
	void SolveBendingConstraints(float weighting);

protected:
	uint m_SolverIterations;
	float m_MaterialBend = 0.99999f;
	float m_MaterialStretchWeft = 0.9999f;
	float m_MaterialStretchWarp = 0.9999f;
	float m_MaterialShear = 0.99f;
	float m_MaterialDamp = 0.001f;

	PArray<Vector3> m_PhyxelPosNew;
	PArray<float> m_PhyxelInvMass; //1/mass (or 0.0 if PhyxelIsStatic==true)

	PArray<PBDDistanceConstraint> m_ConstraintsDistance;
	PArray<PBDBendingConstraint> m_ConstraintsBending;

	std::vector<std::vector<uint>> m_LookupsDistance;
	std::vector<std::vector<uint>> m_LookupsBending;

	PArray<Vector3> m_LargeArrayParallisation;
};