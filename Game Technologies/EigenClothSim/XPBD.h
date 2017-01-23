#pragma once

#include "ClothBase.h"
#define USE_TRIANGLE_BENDING_CONSTRAINT

struct XPBDDistanceConstraint { uint p1, p2;	float rest_length, k, k_prime, lamdaij; };
#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
struct XPBDBendingConstraint { uint p1, p2, p3; float rest_length, w, k, k_prime, lambdaij; };
#else
struct XPBDBendingConstraint { uint p1, p2, p3, p4;	float phi0, rest_length1, rest_length2, w1, w2, k; float k_prime, lambdaij; };
#endif
struct XPBDSphereConstraint { float radius; Vector3 centre; };

class XPBD : public ClothBase
{
public:
	XPBD();
	virtual ~XPBD();

	virtual void InitializeClothDesign(const ClothDesignEntity<3>* design) override;
	void InitializePBDConstraints(int num_width, int num_height);

	inline uint& SolverIterations() { return m_SolverIterations; }
	inline const uint& SolverIterations() const { return m_SolverIterations; }

protected:
	void ComputeVelocity(float subtimestep, const Vector3* in_clothpos, const Vector3* in_clothvel, Vector3* out_clothvel) override;

	void BuildLookupArrays();

	bool BuildDistanceConstraint(uint idxA, uint idxB, float k, XPBDDistanceConstraint& out_constraint);

#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
	bool BuildBendingConstraint(uint idxA, uint idxB, uint idxC, float k, XPBDBendingConstraint& out_constraint);
#else
	bool BuildBendingConstraint(uint idxA, uint idxB, uint idxC, uint idxD, float k, XPBDBendingConstraint& out_constraint);
#endif

	void SolveDistanceConstraint(float weighting, XPBDDistanceConstraint& c);

	void SolveDistanceConstraints(float weighting);
	void SolveBendingConstraints(float weighting);
	void SolveSphereConstraints(float weighting);

public:
	uint m_SolverIterations;
	float m_MaterialBend = 0.75f;
	float m_MaterialStretchWeft = 0.75f;
	float m_MaterialStretchWarp = 0.759f;
	float m_MaterialShear = 0.3f;
	float m_MaterialDamp = 0.1f;

	PArray<Vector3> m_PhyxelPosNew;
	PArray<float> m_PhyxelInvMass; //1/mass (or 0.0 if PhyxelIsStatic==true)

	PArray<XPBDDistanceConstraint> m_ConstraintsDistance;
	PArray<XPBDBendingConstraint> m_ConstraintsBending;
	std::vector<XPBDSphereConstraint> m_SphereConstraints; //sphere collision shapes

	std::vector<uint> m_DistanceBatches;

	std::vector<std::vector<uint>> m_LookupsDistance;
	std::vector<std::vector<uint>> m_LookupsBending;

	PArray<Vector3> m_LargeArrayParallisation;
};