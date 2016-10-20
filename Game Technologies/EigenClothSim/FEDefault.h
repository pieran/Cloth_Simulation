#pragma once

#include "FEBase.h"

class FEDefault : public FEBase
{
public:
	FEDefault();
	virtual ~FEDefault();

	virtual void UpdateTriangleMaterial(const uint& triIdx) override;

	

protected:
	virtual void InitializeClothDesign(const ClothDesignEntity<3>* design) override;

	virtual void GenRotations(const Vector3* in_clothpos) override;
	virtual void GenStressStrainForces(float timestep, const Vector3* in_clothpos, const Vector3* in_clothvel) override;

	void BuildRotationMatrix(Matrix3* out, const Vector3& a, const Vector3& b, const Vector3& c);
	void BuildStiffnessMatricies(uint triIdx);
protected:
	PArray<Matrix3> m_TriRotBase;
	PArray<Matrix3> m_TriRot;

	PArray<Matrix3> m_TriStiffness;
	PArray<Matrix3> m_TriBendStiffness;
};
