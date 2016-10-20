#pragma once

#pragma once

#include "common.h"
#include <nclgl\Vector3.h>
#include <nclgl\Matrix3.h>
#include <vector>
#include <stdio.h>

//Sparse Matrix implementation exploiting the fact that each row will have some data

template<class T>
struct LMAMatrixItem
{
	uint    idxA, idxB, idxC;
	T		Kii[3];
	T		Kij[3];
};


template<class T>
class LMAMatrix
{
public:
	LMAMatrix<T>();
	~LMAMatrix<T>();

	void resize(uint num_elements, uint num_rows);

	void BuildInversePrecondition(std::vector<T>& precondition);

	inline LMAMatrixItem<T>& operator[](uint idx) { return m_Elements[idx]; }
	inline void SetMass(uint idx, float mass) { m_AiiMasses[idx] = mass; }


	void SolveAMultX(std::vector<Vector3>& out_residual, std::vector<Vector3>& out_previous, float& out_r0z0, float& out_beta,
		const std::vector<Matrix3>& constraints, const std::vector<Matrix3>& inv_precondition, const std::vector<Vector3>& bvec, const std::vector<Vector3>& xvec);

	float SolveAMultU(std::vector<Vector3>& out, const std::vector<Matrix3>& constraints, const std::vector<Vector3>& u);


protected:
	std::vector<float>		   m_AiiMasses;
	std::vector<LMAMatrixItem<Matrix3>> m_Elements;
	std::vector<Vector3> m_TmpVec1, m_TmpVec2;

};

template<class T>
LMAMatrix<T>::LMAMatrix()
{

}

template<class T>
LMAMatrix<T>::~LMAMatrix()
{

}

template<class T>
void LMAMatrix<T>::resize(uint num_elements, uint num_rows)
{
	m_Elements.resize(num_elements);
	m_AiiMasses.resize(num_rows);
	m_TmpVec1.resize(num_rows);
	m_TmpVec2.resize(num_rows);
}


template<class T>
void LMAMatrix<T>::SolveAMultX(std::vector<Vector3>& out_residual, std::vector<Vector3>& out_previous, float& out_r0z0, float& out_beta,
	const std::vector<Matrix3>& constraints, const std::vector<Matrix3>& inv_precondition, const std::vector<Vector3>& bvec, const std::vector<Vector3>& xvec)
{
	uint memcpy_size = bvec.size() * sizeof(Vector3);
	memcpy(&out_residual[0], &bvec[0], memcpy_size);
	memcpy(&m_TmpVec2[0], &bvec[0], memcpy_size);

	float beta = 0.0f, r0z0 = 0.0f;
	int elen = (int)m_Elements.size();
	int nlen = (int)out_residual.size();

#pragma omp parallel for
	for (int i = 0; i < nlen; ++i)
	{
		 m_TmpVec1[i].x = (1.0f - constraints[i]._11) * xvec[i].x;
		 m_TmpVec1[i].y = constraints[i]._21 * xvec[i].x;
		 m_TmpVec1[i].z = constraints[i]._31 * xvec[i].x;

		 m_TmpVec1[i].x += constraints[i]._12 * xvec[i].y;
		 m_TmpVec1[i].y += (1.0f - constraints[i]._22) * xvec[i].y;
		 m_TmpVec1[i].z += constraints[i]._32 * xvec[i].y;

		 m_TmpVec1[i].x += constraints[i]._13 * xvec[i].z;
		 m_TmpVec1[i].y += constraints[i]._23 * xvec[i].z;
		 m_TmpVec1[i].z += (1.0f - constraints[i]._33) * xvec[i].z;

		 out_residual[i] -= xvec[i] * m_AiiMasses[i];
		 m_TmpVec2[i] -= m_TmpVec1[i] * m_AiiMasses[i];

	}

#pragma omp parallel for
	for (int i = 0; i < elen; ++i)
	{
		const LMAMatrixItem<T>& ele = m_Elements[i];


		Vector3 out_res[3];
		Vector3 out_tmp[3];
		Vector3 rhs_res[3]{ xvec[ele.idxA], xvec[ele.idxB], xvec[ele.idxC] };
		Vector3 rhs_tmp[3]{ m_TmpVec1[ele.idxA], m_TmpVec1[ele.idxB], m_TmpVec1[ele.idxC] };

		for (int j = 0; j < 3; ++j)
		{
			int k = (j + 1) % 3;

			out_res[j] += ele.Kii[j] * rhs_res[j];
			out_res[j] += ele.Kij[j] * rhs_res[k];
			out_res[k] += ele.Kij[j] * rhs_res[j];

			out_tmp[j] += ele.Kii[j] * rhs_tmp[j];
			out_tmp[j] += ele.Kij[j] * rhs_tmp[k];
			out_tmp[k] += ele.Kij[j] * rhs_tmp[j];
		}

#pragma omp atomic
		out_residual[ele.idxA].x -= out_res[0].x;
#pragma omp atomic
		out_residual[ele.idxA].y -= out_res[0].y;
#pragma omp atomic
		out_residual[ele.idxA].z -= out_res[0].z;

#pragma omp atomic
		out_residual[ele.idxB].x -= out_res[1].x;
#pragma omp atomic
		out_residual[ele.idxB].y -= out_res[1].y;
#pragma omp atomic
		out_residual[ele.idxB].z -= out_res[1].z;

#pragma omp atomic
		out_residual[ele.idxC].x -= out_res[2].x;
#pragma omp atomic
		out_residual[ele.idxC].y -= out_res[2].y;
#pragma omp atomic
		out_residual[ele.idxC].z -= out_res[2].z;


#pragma omp atomic
		m_TmpVec2[ele.idxA].x -= out_tmp[0].x;
#pragma omp atomic
		m_TmpVec2[ele.idxA].y -= out_tmp[0].y;
#pragma omp atomic
		m_TmpVec2[ele.idxA].z -= out_tmp[0].z;

#pragma omp atomic
		m_TmpVec2[ele.idxB].x -= out_tmp[1].x;
#pragma omp atomic
		m_TmpVec2[ele.idxB].y -= out_tmp[1].y;
#pragma omp atomic
		m_TmpVec2[ele.idxB].z -= out_tmp[1].z;

#pragma omp atomic
		m_TmpVec2[ele.idxC].x -= out_tmp[2].x;
#pragma omp atomic
		m_TmpVec2[ele.idxC].y -= out_tmp[2].y;
#pragma omp atomic
		m_TmpVec2[ele.idxC].z -= out_tmp[2].z;
/*#pragma omp atomic
		out_residual[ele.idxB] -= out_res[1];

#pragma omp atomic
		out_residual[ele.idxC] -= out_res[2];

#pragma omp atomic
		m_TmpVec2[ele.idxA] -= out_tmp[0];

#pragma omp atomic
		m_TmpVec2[ele.idxB] -= out_tmp[1];

#pragma omp atomic
		m_TmpVec2[ele.idxC] -= out_tmp[2];*/
	}


#pragma omp parallel for reduction(+:r0z0) reduction(+:beta)
	for (int row = 0; row < nlen; ++row)
	{
		out_residual[row] = constraints[row] * out_residual[row];
		out_previous[row] = constraints[row] * inv_precondition[row] * out_residual[row];
		r0z0 += Vector3::Dot(out_residual[row], out_previous[row]);
		beta += m_TmpVec2[row].Dot(inv_precondition[row] * m_TmpVec2[row]);
	}

	out_beta = beta;
	out_r0z0 = r0z0;
}

template<class T>
float LMAMatrix<T>::SolveAMultU(std::vector<Vector3>& out, const std::vector<Matrix3>& constraints, const std::vector<Vector3>& u)
{
	int elen = (int)m_Elements.size();
	int nlen = (int)out.size();

	memset(&m_TmpVec1[0], 0, nlen * sizeof(Vector3));

	for (int i = 0; i < elen; ++i)
	{
		const LMAMatrixItem<T>& ele = m_Elements[i];

		Vector3 out[3];
		Vector3 rhs[3]{ u[ele.idxA], u[ele.idxB], u[ele.idxC] };

		for (int j = 0; j < 3; ++j)
		{
			int k = (j + 1) % 3;

			out[j] += ele.Kii[j] * rhs[j];
			out[j] += ele.Kij[j] * rhs[k];
			out[k] += ele.Kij[j] * rhs[j];
		}

		m_TmpVec1[ele.idxA] += out[0];
		m_TmpVec1[ele.idxB] += out[1];
		m_TmpVec1[ele.idxC] += out[2];
	}



	float accum = 0.0f;
#pragma omp parallel for reduction(+:accum)// reduction(+:r0z0)
	for (int row = 0; row < nlen; ++row)
	{
		m_TmpVec1[row] += u[row] * m_AiiMasses[row];
		InplaceMatrix3MultVector3(&out[row], constraints[row], m_TmpVec1[row]);
		accum += Vector3::Dot(u[row], out[row]);
	}

	return accum;
}