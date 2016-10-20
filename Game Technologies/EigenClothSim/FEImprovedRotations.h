#pragma once

#include "FEBase.h"

class FEImprovedRotations : public FEBase
{
public:
	FEImprovedRotations(bool useSmoothedNormals = false);
	virtual ~FEImprovedRotations();

protected:
	virtual void InitializeClothDesign(const ClothDesignEntity<3>* design) override;

	virtual void GenRotations(const Vector3* in_clothpos) override;
	virtual void GenStressStrainForces(float timestep, const Vector3* in_clothpos, const Vector3* in_clothvel) override;

	void BuildStiffnessMatricies(uint triIdx);

	void FEImprovedRotations::Diagonalize(const Matrix3&  A, Matrix3& Q, Matrix3& D);
protected:
	
	bool m_UseSmoothedNormals; //Use smoothed normals for calculating out of plane points on triangle

	PArray<Matrix3> m_TriRot;

	PArray<Matrix3> m_TriStiffness;
	PArray<Matrix3> m_TriBendStiffness;
};
