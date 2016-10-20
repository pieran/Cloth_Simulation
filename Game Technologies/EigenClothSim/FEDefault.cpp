#include "FEDefault.h"

FEDefault::FEDefault() : FEBase()
{

}

FEDefault::~FEDefault()
{

}

void FEDefault::UpdateTriangleMaterial(const uint& triIdx)
{
	BuildStiffnessMatricies(triIdx / 3);
}

void FEDefault::InitializeClothDesign(const ClothDesignEntity<3>* design)
{
	FEBase::InitializeClothDesign(design);
	
	m_TriRotBase.resize(m_NumAllocatedTris);
	m_TriRot.resize(m_NumAllocatedTris);

	m_TriStiffness.resize(m_NumAllocatedTris * 9);
	m_TriBendStiffness.resize(m_NumAllocatedTris * 3);

#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTris; ++i)
	{
		BuildStiffnessMatricies(i);
	}
}

void FEDefault::GenRotations(const Vector3* in_clothpos)
{
	Matrix3 mtx;
	Vector3 a, b, c;

//#pragma omp parallel for private(mtx, a, b, c)
	for (int i = 0; i < (int)m_NumTris; ++i)
	{
		int i3 = i * 3;
		uint idxA = m_TriIndicies[i3];
		uint idxB = m_TriIndicies[i3 + 1];
		uint idxC = m_TriIndicies[i3 + 2];

		a = in_clothpos[idxA];
		b = in_clothpos[idxB];
		c = in_clothpos[idxC];

		BuildRotationMatrix(&mtx, a, b, c);
		InplaceMatrix3MultMatrix3(&m_TriRot[i], mtx, m_TriRotBase[i]);
	}
}

void FEDefault::GenStressStrainForces(float timestep, const Vector3* in_clothpos, const Vector3* in_clothvel)
{
	const float dt2 = timestep * timestep;
	const float dt2_visos = (timestep + V_SCALAR) * timestep;

	Matrix3 RK, RKR, B_rot, Q, ReT;
	Vector3 posAB, posIAB;

//#pragma omp parallel for private(RK, RKR, ReT, B_rot, Q, posAB, posIAB)
	for (int i = 0; i < (int)m_NumTris; ++i)
	{
		int i3 = i * 3;
		uint triIndicies[3] = { m_TriIndicies[i3], m_TriIndicies[i3 + 1], m_TriIndicies[i3 + 2] };

		Vector3* B_i = &m_TriB[i3];
		Matrix3* A_ii = &m_TriAii[i3];
		Matrix3* A_ij = &m_TriAij[i3];
		memset(&A_ii[0], 0, 3 * sizeof(Matrix3));
		memset(&A_ij[0], 0, 3 * sizeof(Matrix3));


		Vector3 triNormals[3] = { m_PhyxelNormals[triIndicies[0]], m_PhyxelNormals[triIndicies[1]], m_PhyxelNormals[triIndicies[2]] };

		Matrix3 Pa = Matrix3::OuterProduct(triNormals[0], triNormals[0]);
		Matrix3 Pb = Matrix3::OuterProduct(triNormals[1], triNormals[1]);
		Matrix3 Pc = Matrix3::OuterProduct(triNormals[2], triNormals[2]);


		const Matrix3& Re = m_TriRot[i];
		Matrix3 ReT = Matrix3::Transpose(Re);

		Vector3 tmpMult;

		for (uint j = 0; j < 3; ++j)
		{
			for (uint k = 0; k < 3; ++k)
			{
				//const Vector3& dNa = dN[j];
				//const Vector3& dNb = dN[k];
				Matrix3& Ke = m_TriStiffness[i * 9 + j * 3 + k];

				const int idxA = triIndicies[j];
				const int idxB = triIndicies[k];

				posAB = in_clothpos[idxA] - in_clothpos[idxB];
				posIAB = m_PhyxelPosInitial[idxA] - m_PhyxelPosInitial[idxB];

				//RK = Re * ke;
				//RKR = Re * ke * ReT;
				InplaceMatrix3MultMatrix3(&RK, Re, Ke);
				InplaceMatrix3MultMatrix3(&RKR, RK, ReT);

				if (j != k)
				{
					InplaceMatrix3MultVector3(&tmpMult, RKR, posAB);
					B_i[j] += tmpMult;
					//B_i[k] -= tmpMult;

					InplaceMatrix3MultVector3(&tmpMult, RK, posIAB);
					B_i[j] -= tmpMult;
					//B_i[k] += tmpMult;

					RKR *= dt2_visos;
					if (j < k)
					A_ij[j] += RKR;
				}
				else if (j == k)
				{
					RKR *= dt2_visos;
					A_ii[j] -= RKR;
					//A_ii[k] -= RKR;
				}			
			}
		}

		for (uint j = 0; j < 3; ++j)
		{
			uint k = (j + 1) % 3;

			const int idxA = triIndicies[j];
			const int idxB = triIndicies[k];

			posAB = in_clothpos[idxA] - in_clothpos[idxB];
			posIAB = m_PhyxelPosInitial[idxA] - m_PhyxelPosInitial[idxB];

			//RK = Re * ke;
			//RKR = Re * ke * ReT;
			/*InplaceMatrix3MultMatrix3(&RK, Re, m_TriStiffness[i3 + j]);
			InplaceMatrix3MultMatrix3(&RKR, RK, ReT);

			InplaceMatrix3MultVector3(&tmpMult, RKR, posAB);
			B_i[j] += tmpMult;
			B_i[k] -= tmpMult;

			InplaceMatrix3MultVector3(&tmpMult, RK, posIAB);
			B_i[j] -= tmpMult;
			B_i[k] += tmpMult;

			RKR *= dt2_visos;
			A_ii[j] -= RKR;
			A_ii[k] -= RKR;
			A_ij[j] += RKR;*/


			//Bending			
			/*B_rot = Re * t.B[j] * ReT;
			Q = Pa * B_rot * Pa
			+ Pb * B_rot * Pb
			+ Pc * B_rot * Pc;
			Q /= 3.f;*/
			/*InplaceMatrix3MultMatrix3(&RK, Re, m_TriBendStiffness[i3 + j]);
			InplaceMatrix3MultMatrix3(&B_rot, RK, ReT);

			InplaceMatrix3MultMatrix3(&RK, Pa, B_rot);
			InplaceMatrix3MultMatrix3(&Q, RK, Pa);

			InplaceMatrix3MultMatrix3(&RK, Pb, B_rot);
			InplaceMatrix3MultMatrix3Additive(&Q, RK, Pb);

			InplaceMatrix3MultMatrix3(&RK, Pc, B_rot);
			InplaceMatrix3MultMatrix3Additive(&Q, RK, Pc);

			Q *= (1.0f / 3.0f);

			InplaceMatrix3MultVector3(&tmpMult, Q, posAB);
			B_i[j] -= tmpMult;
			B_i[k] += tmpMult;

			Q *= dt2;
			A_ii[j] -= Q;
			A_ii[k] -= Q;
			A_ij[j] += Q;*/
		}
	}
}

void FEDefault::BuildRotationMatrix(Matrix3* out, const Vector3& a, const Vector3& b, const Vector3& c)
{
	Vector3 t_rows[3];
	t_rows[0] = a - c;
	t_rows[1] = b - c;
	t_rows[2] = Vector3::Cross(t_rows[0], t_rows[1]);

	//Build Rotation Matrtix
/*	t_rows[0].Normalise();
	t_rows[1] -= t_rows[0] * Vector3::Dot(t_rows[0], t_rows[1]);
	t_rows[2] -= t_rows[0] * Vector3::Dot(t_rows[0], t_rows[2]);

	t_rows[1].Normalise();
	t_rows[2] -= t_rows[1] * Vector3::Dot(t_rows[1], t_rows[2]);

	t_rows[2].Normalise();*/

	t_rows[2].Normalise();
	t_rows[0] -= t_rows[2] * Vector3::Dot(t_rows[2], t_rows[0]);
	t_rows[1] -= t_rows[2] * Vector3::Dot(t_rows[2], t_rows[1]);

	t_rows[1].Normalise();
	t_rows[0] -= t_rows[1] * Vector3::Dot(t_rows[1], t_rows[0]);

	t_rows[0].Normalise();

	(*out) = Matrix3(t_rows[0], t_rows[1], t_rows[2]);
}

void FEDefault::BuildStiffnessMatricies(uint triIdx)
{
	int i, j;
	Vector3 pA = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3]];
	Vector3 pB = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3 + 1]];
	Vector3 pC = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3 + 2]];

	BuildRotationMatrix(&m_TriRotBase[triIdx], pA, pB, pC);
	m_TriRotBase[triIdx] = Matrix3::Transpose(m_TriRotBase[triIdx]);


	//pA = m_TriRotBase[triIdx] * pA;
	//pB = m_TriRotBase[triIdx] * pB;
	//pC = m_TriRotBase[triIdx] * pC;

	float triArea = 0.5f * abs((pA.x - pC.x) * (pB.y - pC.y)
								- (pA.y - pC.y) * (pB.x - pC.x));

	Vector3 dN[3];


	//Computes Shape Function and Derivitives (Assumes triangle only exists in the x/y planes)
	Vector2 r1 = Vector2(pA.x - pC.x, pB.x - pC.x);
	Vector2 r2 = Vector2(pA.y - pC.y, pB.y - pC.y);

	float detJ = (float)(1.0 / (r1.x * r2.y - r1.y * r2.x));

	dN[0] = Vector3(r2.y, -r1.y, 0.0f) * detJ;
	dN[1] = Vector3(-r2.x, r1.x, 0.0f) * detJ;
	dN[2] = -dN[0] - dN[1];

	/*float det = Vector3::Dot(e1, e1) * Vector3::Dot(e2, e2) - Vector3::Dot(e1, e2) * Vector3::Dot(e1, e2);
	dN[0] = (e1 * Vector3::Dot(e2, e2) - e2 * Vector3::Dot(e1, e2)) / det;
	dN[1] = (e2 * Vector3::Dot(e1, e1) - e1 * Vector3::Dot(e1, e2)) / det;
	dN[2] = -tri.dN[0] - tri.dN[1];*/





	//Stiffness Matrix for each vertex pair
	for (i = 0; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			const Vector3& dNa = dN[i];
			const Vector3& dNb = dN[j];

			Matrix3& Ke = m_TriStiffness[triIdx * 9 + i * 3 + j];

			float lambda = 1102.0f;// 3000.f;// 1366.f;
			float mu = 1000.f;//3000.f;//800.f;

			Matrix3 tmp2 = Matrix3::OuterProduct(dNa, dNb);// +Matrix3::OuterProduct(dNb, dNa)) * 0.5f;
			Ke = (Matrix3::Identity * (lambda * tmp2.Trace())) + tmp2 * (mu * 2.0f);
			Ke(2, 2) = 0.f;
			Ke *= triArea;
		}
	}

	//Bending Matrix for Each Vertex Pair
	for (i = 0; i < 3; ++i)
	{
		j = (i + 1) % 3;

		Matrix3& B = m_TriBendStiffness[triIdx * 3 + i];
		B = Matrix3::ZeroMatrix;

		const Vector3& dNa = dN[i];
		const Vector3& dNb = dN[j];

		B(0, 0) = dNa.x * dNb.x * B1 * triArea;
		B(1, 1) = dNa.y * dNb.y * B2 * triArea;
	}
}

