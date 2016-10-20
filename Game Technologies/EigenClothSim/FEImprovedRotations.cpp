#include "FEImprovedRotations.h"

FEImprovedRotations::FEImprovedRotations(bool useSmoothedNormals) : FEBase()
{
	m_UseSmoothedNormals = useSmoothedNormals;
}

FEImprovedRotations::~FEImprovedRotations()
{

}


void FEImprovedRotations::InitializeClothDesign(const ClothDesignEntity<3>* design)
{
	FEBase::InitializeClothDesign(design);

	m_TriRot.resize(m_NumTris);

	m_TriStiffness.resize(m_NumTris * 3);
	m_TriBendStiffness.resize(m_NumTris * 3);

#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTris; ++i)
	{
		BuildStiffnessMatricies(i);
	}
}

// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0
typedef double Real;
void FEImprovedRotations::Diagonalize(const Matrix3&  Amtx, Matrix3& Qmtx, Matrix3& Dmtx)
{
	// A must be a symmetric matrix.
	// returns Q and D such that 
	// Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
	const int maxsteps = 24;  // certainly wont need that many.
	int k0, k1, k2;
	Real o[3], m[3];
	Real q[4] = { 0.0,0.0,0.0,1.0 };
	Real jr[4];
	Real sqw, sqx, sqy, sqz;
	Real tmp1, tmp2, mq;
	Real AQ[3][3];
	Real thet, sgn, t, c;
	Real A[9], Q[9], D[9];
	
	for (int i = 0; i < 9; i++)
	{
		A[i] = (double)Amtx.mat_array[i];
		Q[i] = (double)Qmtx.mat_array[i];
		D[i] = (double)Dmtx.mat_array[i];
	}

	for (int i = 0; i < maxsteps; ++i)
	{
		// quat to matrix
		sqx = q[0] * q[0];
		sqy = q[1] * q[1];
		sqz = q[2] * q[2];
		sqw = q[3] * q[3];
		Q[0] = (sqx - sqy - sqz + sqw);
		Q[3 + 1] = (-sqx + sqy - sqz + sqw);
		Q[6 + 2] = (-sqx - sqy + sqz + sqw);
		tmp1 = q[0] * q[1];
		tmp2 = q[2] * q[3];
		Q[3 + 0] = 2.0 * (tmp1 + tmp2);
		Q[0 + 1] = 2.0 * (tmp1 - tmp2);
		tmp1 = q[0] * q[2];
		tmp2 = q[1] * q[3];
		Q[6 + 0] = 2.0 * (tmp1 - tmp2);
		Q[0 + 2] = 2.0 * (tmp1 + tmp2);
		tmp1 = q[1] * q[2];
		tmp2 = q[0] * q[3];
		Q[6 + 1] = 2.0 * (tmp1 + tmp2);
		Q[3 + 2] = 2.0 * (tmp1 - tmp2);

		// AQ = A * Q
		AQ[0][0] = Q[0 + 0] * A[0 + 0] + Q[3 + 0] * A[0 + 1] + Q[6 + 0] * A[0 + 2];
		AQ[0][1] = Q[0 + 1] * A[0 + 0] + Q[3 + 1] * A[0 + 1] + Q[6 + 1] * A[0 + 2];
		AQ[0][2] = Q[0 + 2] * A[0 + 0] + Q[3 + 2] * A[0 + 1] + Q[6 + 2] * A[0 + 2];
		AQ[1][0] = Q[0 + 0] * A[0 + 1] + Q[3 + 0] * A[3 + 1] + Q[6 + 0] * A[3 + 2];
		AQ[1][1] = Q[0 + 1] * A[0 + 1] + Q[3 + 1] * A[3 + 1] + Q[6 + 1] * A[3 + 2];
		AQ[1][2] = Q[0 + 2] * A[0 + 1] + Q[3 + 2] * A[3 + 1] + Q[6 + 2] * A[3 + 2];
		AQ[2][0] = Q[0 + 0] * A[0 + 2] + Q[3 + 0] * A[3 + 2] + Q[6 + 0] * A[6 + 2];
		AQ[2][1] = Q[0 + 1] * A[0 + 2] + Q[3 + 1] * A[3 + 2] + Q[6 + 1] * A[6 + 2];
		AQ[2][2] = Q[0 + 2] * A[0 + 2] + Q[3 + 2] * A[3 + 2] + Q[6 + 2] * A[6 + 2];
		// D = Qt * AQ
		D[0 + 0] = AQ[0][0] * Q[0 + 0] + AQ[1][0] * Q[3 + 0] + AQ[2][0] * Q[6 + 0];
		D[0 + 1] = AQ[0][0] * Q[0 + 1] + AQ[1][0] * Q[3 + 1] + AQ[2][0] * Q[6 + 1];
		D[0 + 2] = AQ[0][0] * Q[0 + 2] + AQ[1][0] * Q[3 + 2] + AQ[2][0] * Q[6 + 2];
		D[3 + 0] = AQ[0][1] * Q[0 + 0] + AQ[1][1] * Q[3 + 0] + AQ[2][1] * Q[6 + 0];
		D[3 + 1] = AQ[0][1] * Q[0 + 1] + AQ[1][1] * Q[3 + 1] + AQ[2][1] * Q[6 + 1];
		D[3 + 2] = AQ[0][1] * Q[0 + 2] + AQ[1][1] * Q[3 + 2] + AQ[2][1] * Q[6 + 2];
		D[6 + 0] = AQ[0][2] * Q[0 + 0] + AQ[1][2] * Q[3 + 0] + AQ[2][2] * Q[6 + 0];
		D[6 + 1] = AQ[0][2] * Q[0 + 1] + AQ[1][2] * Q[3 + 1] + AQ[2][2] * Q[6 + 1];
		D[6 + 2] = AQ[0][2] * Q[0 + 2] + AQ[1][2] * Q[3 + 2] + AQ[2][2] * Q[6 + 2];
		o[0] = D[3 + 2];
		o[1] = D[0 + 2];
		o[2] = D[0 + 1];
		m[0] = fabs(o[0]);
		m[1] = fabs(o[1]);
		m[2] = fabs(o[2]);

		k0 = (m[0] > m[1] && m[0] > m[2]) ? 0 : (m[1] > m[2]) ? 1 : 2; // index of largest element of offdiag
		k1 = (k0 + 1) % 3;
		k2 = (k0 + 2) % 3;
		if (o[k0] == 0.0)
		{
			break;  // diagonal already
		}
		thet = (D[k2 * 3 + k2] - D[k1 * 3 + k1]) / (2.0*o[k0]);
		sgn = (thet > 0.0) ? 1.0 : -1.0;
		thet *= sgn; // make it positive
		t = sgn / (thet + ((thet < 1.E6) ? sqrt(thet*thet + 1.0) : thet)); // sign(T)/(|T|+sqrt(T^2+1))
		c = 1.0 / sqrt(t*t + 1.0); //  c= 1/(t^2+1) , t=s/c 
		if (c == 1.0)
		{
			break;  // no room for improvement - reached machine precision.
		}
		jr[0] = jr[1] = jr[2] = jr[3] = 0.0;
		jr[k0] = sgn*sqrt((1.0 - c) / 2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)  
		jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
		jr[3] = sqrt(1.0f - jr[k0] * jr[k0]);
		if (jr[3] == 1.0)
		{
			break; // reached limits of floating point precision
		}
		q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
		q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
		q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
		q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
		mq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		q[0] /= mq;
		q[1] /= mq;
		q[2] /= mq;
		q[3] /= mq;
	}

	for (int i = 0; i < 9; i++)
	{
		Qmtx.mat_array[i] = (float)Q[i];
		Dmtx.mat_array[i] = (float)D[i];
	}
}

void FEImprovedRotations::GenRotations(const Vector3* in_clothpos)
{
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTris; ++i)
	{
		int i3 = i * 3;
		uint idxA = m_TriIndicies[i3];
		uint idxB = m_TriIndicies[i3 + 1];
		uint idxC = m_TriIndicies[i3 + 2];

		//SHAPE MATCHING!!!
		Vector3 pos[3]{
			in_clothpos[idxA],
			in_clothpos[idxB],
			in_clothpos[idxC]
		};

		Vector3 posi[3]{
			m_PhyxelPosInitial[idxA],
			m_PhyxelPosInitial[idxB],
			m_PhyxelPosInitial[idxC]
		};

		Vector3 center = (pos[0] + pos[1] + pos[2]) / 3.0f;
		Vector3 centeri = (posi[0] + posi[1] + posi[2]) / 3.0f;

		Matrix3 Apq = Matrix3::ZeroMatrix;
		for (int i = 0; i < 3; ++i)
		{
			Apq += Matrix3::OuterProduct(pos[i] - center, posi[i] - centeri);
		}

		//Compute Out of Plane Point (centre + normal)
		Vector3 normal;
		if (m_UseSmoothedNormals)
		{
			normal = m_PhyxelNormals[idxA] + m_PhyxelNormals[idxB] + m_PhyxelNormals[idxC];		
		}
		else
		{
			normal = Vector3::Cross(pos[0] - pos[2], pos[1] - pos[2]);
		}
		normal.Normalise();

		Apq += Matrix3::OuterProduct(normal, Vector3(0, 0, 1));

		Apq /= 4.0f;


		Matrix3 S2 = Matrix3::Transpose(Apq) * Apq;

		Matrix3 Q, D;
		Diagonalize(S2, Q, D);
		const float epsilon = 1E-6f;
		Matrix3 sD(sqrtf(D._11), 0.0f, 0.0f,
			0.0f, sqrtf(D._22), 0.0f,
			0.0f, 0.0f, sqrtf(D._33));

		Matrix3 S = Matrix3::Transpose(Q) * sD * Q;
		Matrix3 rot = Apq * Matrix3::Inverse(S);

		S2 = Matrix3::Transpose(rot) * rot;
		Diagonalize(S2, Q, D);
		sD = Matrix3(sqrtf(D._11), 0.0f, 0.0f,
			0.0f, sqrtf(D._22), 0.0f,
			0.0f, 0.0f, sqrtf(D._33));

		S = Matrix3::Transpose(Q) * sD * Q;
		rot = rot * Matrix3::Inverse(S);


		m_TriRot[i] = rot;
	}
}

void FEImprovedRotations::GenStressStrainForces(float timestep, const Vector3* in_clothpos, const Vector3* in_clothvel)
{
	const float dt2_visos = timestep*timestep + V_SCALAR * timestep;
	const float dt2 = timestep*timestep;

	Matrix3 RK, RKR, B_rot, Q, ReT;
	Vector3 posAB, posIAB;

#pragma omp parallel for private(RK, RKR, ReT, B_rot, Q, posAB, posIAB)
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

		Vector3 face_normal = (triNormals[0] + triNormals[1] + triNormals[2]);
		face_normal.Normalise();


		const Matrix3& Re = m_TriRot[i];
		Matrix3 ReT = Matrix3::Transpose(Re);

		Vector3 tmpMult;


		for (uint j = 0; j < 3; ++j)
		{
			uint k = (j + 1) % 3;

			const int idxA = triIndicies[j];
			const int idxB = triIndicies[k];

			posAB = in_clothpos[idxA] - in_clothpos[idxB];
			posIAB = m_PhyxelPosInitial[idxA] - m_PhyxelPosInitial[idxB];

			//RK = Re * ke;
			//RKR = Re * ke * ReT;
			InplaceMatrix3MultMatrix3(&RK, Re, m_TriStiffness[i3 + j]);
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
			A_ij[j] += RKR;

			//Bending			
			/*B_rot = Re * t.B[j] * ReT;
			Q = Pa * B_rot * Pa
			+ Pb * B_rot * Pb
			+ Pc * B_rot * Pc;
			Q /= 3.f;*/
			InplaceMatrix3MultMatrix3(&RK, Re, m_TriBendStiffness[i3 + j]);
			InplaceMatrix3MultMatrix3(&B_rot, RK, ReT);

			InplaceMatrix3MultMatrix3(&RK, Pa, B_rot);
			InplaceMatrix3MultMatrix3(&Q, RK, Pa);

			InplaceMatrix3MultMatrix3(&RK, Pb, B_rot);
			InplaceMatrix3MultMatrix3Additive(&Q, RK, Pb);
			InplaceMatrix3MultMatrix3(&RK, Pc, B_rot);
			InplaceMatrix3MultMatrix3Additive(&Q, RK, Pc);
			Q *= (1.0f / 3.0f);

			InplaceMatrix3MultVector3(&tmpMult, Q, posAB);
			B_i[j] += tmpMult;
			B_i[k] -= tmpMult;


			Q *= dt2;
			A_ii[j] -= Q;
			A_ii[k] -= Q;
			A_ij[j] += Q;
		}
	}
}

void FEImprovedRotations::BuildStiffnessMatricies(uint triIdx)
{
	int i, j;
	Vector3 pA = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3]];
	Vector3 pB = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3 + 1]];
	Vector3 pC = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3 + 2]];


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



	{
		/*Vector3 a = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3]];
		Vector3 b = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3 + 1]];
		Vector3 c = m_PhyxelPosInitial[m_TriIndicies[triIdx * 3 + 2]];

		Vector3 e1 = a - c;
		Vector3 e2 = b - c;
		Vector3 e3 = a - b;


		float detJ = e1.x * e2.y - e1.y - e2.x;
		Eigen::MatrixXf B(3, 6);
		B.setZero();
		B(0, 0) = e2.y;
		B(2, 0) = -e2.x;
		B(1, 1) = -e2.x;
		B(2, 1) = e2.y;
		B(0, 2) = -e1.y;
		B(2, 2) = e1.x;
		B(1, 3) = e1.x;
		B(2, 3) = -e1.y;
		B(0, 4) = e3.y;
		B(2, 4) = -e3.x;
		B(1, 5) = -e3.x;
		B(2, 5) = e3.y;
		B *= 1.0f / detJ;

		float area = 0.5f * detJ;


		const float E = 1000.0f;		//Youngs Modulus
		const float v = 0.3f;			//Poisson coefficient
		Eigen::Matrix3f D;
		D.setZero();
		D(0, 0) = 1.0f;
		D(1, 0) = v;
		D(0, 1) = v;
		D(1, 1) = 1.0f;
		D(2, 2) = (1.0f - v) * 0.5f;
		D *= E / (1.0f - v * v);

		Eigen::MatrixXf eke = B.transpose() * D * B * area;

		/*for (i = 0; i < 3; ++i)
		{
			j = (i + 1) % 3;

			const Vector3& dNa = dN[i];
			const Vector3& dNb = dN[j];
			Matrix3& Ke = m_TriStiffness[triIdx * 3 + i];
			Ke.ToZero();

			for (int k = 0; k < 2; ++k)
			{
				for (int l = 0; l < 2; ++l)
				{
					Ke(k, l) = eke(k + 2 * i, l + 2 * j);
				}
			}
		}*/

		/*auto copy_matrix22 = [](Matrix3& dst, const Eigen::MatrixXf& src, int r_offset, int c_offset)
		{
			dst.ToZero();
			for (int k = 0; k < 2; ++k)
			{
				for (int l = 0; l < 2; ++l)
				{
					dst(k, l) = src(k + r_offset, l + c_offset);
				}
			}
		};
		uint triIdx3 = triIdx * 3;
		copy_matrix22(m_TriStiffness[triIdx3], eke, 0, 2);
		copy_matrix22(m_TriStiffness[triIdx3 + 1], eke, 2, 4);
		copy_matrix22(m_TriStiffness[triIdx3 + 1], eke, 0, 4);*/
	}


	//{

	//	/*const float C1111 = 245;			// Weft Stretch
	//	const float C2222 = 366;			// Warp Stretch
	//	const float C1212 = 2.38;			// Shear Modulus
	//	const float C1122 = 61.1;			// Transverse Contraction*/

	//	const float C1111 = 3057;			// Weft Stretch
	//	const float C2222 = 1534;			// Warp Stretch
	//	const float C1212 = 1.22;			// Shear Modulus
	//	const float C1122 = 459.1;			// Transverse Contraction

	//	const float C2121 = C1212;			// Shear Modulus
	//	const float C2112 = C1212;
	//	const float C1221 = C1212;
	//	const float C2211 = C1122;			// Transverse Contraction
	//	const float C2122 = 0.f;
	//	const float C1211 = 0.f;
	//	const float C1112 = 0.f;
	//	const float C2111 = 0.f;
	//	const float C2221 = 0.f;
	//	const float C2212 = 0.f;
	//	const float C1121 = 0.f;
	//	const float C1222 = 0.f;

	//	uint triIdx3 = triIdx * 3;
	//	Matrix3& ke_ab = m_TriStiffness[triIdx3];
	//	Matrix3& ke_bc = m_TriStiffness[triIdx3 + 1];
	//	Matrix3& ke_ca = m_TriStiffness[triIdx3 + 2];


	//	ke_ab(0, 0) = dN[0].x * C1111 * dN[1].x
	//		+ dN[0].x * C1112 * dN[1].y
	//		+ dN[0].x * C1211 * dN[1].x
	//		+ dN[0].x * C1212 * dN[1].y;
	//	ke_ab(0, 1) = dN[0].x * C1121 * dN[1].x
	//		+ dN[0].x * C1122 * dN[1].y
	//		+ dN[0].x * C1221 * dN[1].x
	//		+ dN[0].x * C1222 * dN[1].y;
	//	ke_ab(1, 0) = dN[0].y * C2111 * dN[1].x
	//		+ dN[0].y * C2112 * dN[1].y
	//		+ dN[0].y * C2211 * dN[1].x
	//		+ dN[0].y * C2212 * dN[1].y;
	//	ke_ab(1, 1) = dN[0].y * C2121 * dN[1].x
	//		+ dN[0].y * C2122 * dN[1].y
	//		+ dN[0].y * C2221 * dN[1].x
	//		+ dN[0].y * C2222 * dN[1].y;

	//	ke_bc(0, 0) = dN[1].x * C1111 * dN[2].x
	//		+ dN[1].x * C1112 * dN[2].y
	//		+ dN[1].x * C1211 * dN[2].x
	//		+ dN[1].x * C1212 * dN[2].y;
	//	ke_bc(0, 1) = dN[1].x * C1121 * dN[2].x
	//		+ dN[1].x * C1122 * dN[2].y
	//		+ dN[1].x * C1221 * dN[2].x
	//		+ dN[1].x * C1222 * dN[2].y;
	//	ke_bc(1, 0) = dN[1].y * C2111 * dN[2].x
	//		+ dN[1].y * C2112 * dN[2].y
	//		+ dN[1].y * C2211 * dN[2].x
	//		+ dN[1].y * C2212 * dN[2].y;
	//	ke_bc(1, 1) = dN[1].y * C2121 * dN[2].x
	//		+ dN[1].y * C2122 * dN[2].y
	//		+ dN[1].y * C2221 * dN[2].x
	//		+ dN[1].y * C2222 * dN[2].y;


	//	/*ke_ca(0, 0) = dN[2].x * C1111 * dN[0].x
	//		+ dN[2].x * C1112 * dN[0].y
	//		+ dN[2].x * C1211 * dN[0].x
	//		+ dN[2].x * C1212 * dN[0].y;
	//	ke_ca(0, 1) = dN[2].x * C1121 * dN[0].x
	//		+ dN[2].x * C1122 * dN[0].y
	//		+ dN[2].x * C1221 * dN[0].x
	//		+ dN[2].x * C1222 * dN[0].y;
	//	ke_ca(1, 0) = dN[2].y * C2111 * dN[0].x
	//		+ dN[2].y * C2112 * dN[0].y
	//		+ dN[2].y * C2211 * dN[0].x
	//		+ dN[2].y * C2212 * dN[0].y;
	//	ke_ca(1, 1) = dN[2].y * C2121 * dN[0].x
	//		+ dN[2].y * C2122 * dN[0].y
	//		+ dN[2].y * C2221 * dN[0].x
	//		+ dN[2].y * C2222 * dN[0].y;*/

	//	ke_ca(0, 0) = dN[0].x * C1111 * dN[2].x
	//		+ dN[0].x * C1112 * dN[2].y
	//		+ dN[0].x * C1211 * dN[2].x
	//		+ dN[0].x * C1212 * dN[2].y;
	//	ke_ca(0, 1) = dN[0].x * C1121 * dN[2].x
	//		+ dN[0].x * C1122 * dN[2].y
	//		+ dN[0].x * C1221 * dN[2].x
	//		+ dN[0].x * C1222 * dN[2].y;
	//	ke_ca(1, 0) = dN[0].y * C2111 * dN[2].x
	//		+ dN[0].y * C2112 * dN[2].y
	//		+ dN[0].y * C2211 * dN[2].x
	//		+ dN[0].y * C2212 * dN[2].y;
	//	ke_ca(1, 1) = dN[0].y * C2121 * dN[2].x
	//		+ dN[0].y * C2122 * dN[2].y
	//		+ dN[0].y * C2221 * dN[2].x
	//		+ dN[0].y * C2222 * dN[2].y;

	//	ke_ab *= triArea;
	//	ke_bc *= triArea;
	//	ke_ca *= triArea;
	//}

	//Stiffness Matrix for each vertex pair
	for (i = 0; i < 3; ++i)
	{
		j = (i + 1) % 3;

		const Vector3& dNa = dN[i];
		const Vector3& dNb = dN[j];

		Matrix3& Ke = m_TriStiffness[triIdx * 3 + i];
		//Matrix3 comp = m_TriStiffness[triIdx * 3 + i];
		float lambda =  1366.f;
		float mu = 800.f;

		Matrix3 tmp2 = (Matrix3::OuterProduct(dNa, dNb) + Matrix3::OuterProduct(dNb, dNa)) * 0.5f;
		Ke = (Matrix3::Identity * (lambda * tmp2.Trace())) + tmp2 * (mu * 2.0f);
		Ke(2, 2) = 0.f;
		Ke *= triArea;

		//printf("moo");
	}


	//Bending Matrix for Each Vertex Pair
	for (i = 0; i < 3; ++i)
	{
		j = (i + 1) % 3;

		Matrix3& B = m_TriBendStiffness[triIdx * 3 + i];

		const Vector3& dNa = dN[i];
		const Vector3& dNb = dN[j];

		B(0, 0) = dNa.x * dNb.x * B1 *triArea;
		B(1, 1) = dNa.y * dNb.y * B2 *triArea;
		B(2, 2) = 0.0f;
	}
}
