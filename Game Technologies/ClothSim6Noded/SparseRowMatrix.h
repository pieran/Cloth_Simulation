#pragma once

#include <nclgl\Vector3.h>
#include <nclgl\Matrix3.h>
#include <vector>
#include <stdio.h>
#include "SimulationDefines.h"

//Sparse Matrix implementation exploiting the fact that each row will have some data

template<class T>
struct SparseRowMatrixItem
{
	uint						  column;
	T							  value;
};

template<class T>
struct SparseRowMatrixItem_SortAscending
{
	inline bool operator()(const SparseRowMatrixItem<T>& a, const SparseRowMatrixItem<T>& b)
	{
		return a.column < b.column;
	}
};



template<class T>
class SparseRowMatrix
{
public:
	SparseRowMatrix<T>();
	~SparseRowMatrix<T>();

	void resize(uint size);
	void organise_memory();
	void clean_memory();
	void zero_memory();

	void SetItem(uint row, uint col, const T& value);
	T& operator()(uint row, uint col);
	const T& operator()(uint row, uint col) const { return *this[row, col]; }


	inline std::vector<SparseRowMatrixItem<T>>& GetRow(uint row) { return m_Rows[row]; }



	void SolveAMultX(std::vector<Vector3>& out_residual, std::vector<Vector3>& out_previous, float& out_r0z0, float& out_beta,
		const std::vector<Matrix3>& constraints, const std::vector<Matrix3>& inv_precondition, const std::vector<Vector3>& bvec, const std::vector<Vector3>& xvec);

	float SolveAMultU(std::vector<Vector3>& out, const std::vector<Matrix3>& constraints, const std::vector<Vector3>& u);




	std::vector<std::vector<SparseRowMatrixItem<T>>> m_Rows;
};

template<class T>
SparseRowMatrix<T>::SparseRowMatrix()
{

}

template<class T>
SparseRowMatrix<T>::~SparseRowMatrix()
{

}

template<class T>
void SparseRowMatrix<T>::resize(uint size)
{
	clean_memory();
	m_Rows.resize(size);
}

template<class T>
void SparseRowMatrix<T>::clean_memory()
{
	m_Rows.clear();
}

template<class T>
T& SparseRowMatrix<T>::operator()(uint row, uint col)
{
	//std::vector<SparseRowMatrixItem<T>>& trow = m_Rows[row];

	for (uint i = 0, len = m_Rows[row].size(); i < len; ++i)
	{
		if (m_Rows[row][i].column == col)
			return m_Rows[row][i].value;
	}

	SparseRowMatrixItem<T> item;
	item.column = col;
	m_Rows[row].push_back(item);

	return m_Rows[row][m_Rows[row].size() - 1].value;
}

template<class T>
void SparseRowMatrix<T>::SetItem(uint row, uint col, const T& value)
{
	std::vector<SparseRowMatrixItem<T>>& row = m_Rows[row];

	for (SparseRowMatrixItem<T>& item : row)
	{
		if (item.column == col)
		{
			item = value;
			return;
		}
	}

	SparseRowMatrixItem<T> item;
	item.coloumn = col;
	item.value = value;
	row.push_back(item);
}

template<class T>
void SparseRowMatrix<T>::organise_memory()
{
	for (uint i = 0, len = m_Rows.size(); ++i)
	{
		std::sort(m_Rows[i].begin(), m_Rows[i].end(), SparseRowMatrixItem_SortAscending());
	}
}

template<class T>
void SparseRowMatrix<T>::zero_memory()
{
	int len = (int)m_Rows.size();
	for (int i = 0; i < len; ++i)
	{
		for (auto itr = m_Rows[i].begin(), end = m_Rows[i].end(); itr != end; ++itr)
		{
			memset(&itr->value, 0, sizeof(T));
		}
	}
}





template<class T>
void SparseRowMatrix<T>::SolveAMultX(std::vector<Vector3>& out_residual, std::vector<Vector3>& out_previous, float& out_r0z0, float& out_beta,
	const std::vector<Matrix3>& constraints, const std::vector<Matrix3>& inv_precondition, const std::vector<Vector3>& bvec, const std::vector<Vector3>& xvec)
{
	int len = (int)m_Rows.size();
	float beta = 0.0f, r0z0 = 0.0f;

#pragma omp parallel for reduction(+:beta) reduction(+:r0z0)
	for (int row = 0; row < len; ++row) {
		Vector3 tmpRes = bvec[row];
		Vector3 tmpBeta = bvec[row];

		Vector3 tmpVec;
		Matrix3 tmp;

		auto itr = m_Rows[row].begin(), end = m_Rows[row].end();
		for (; itr != end; itr++)
		{
			uint col = itr->column;
			/*tmpRes -= itr->value * m_X[col];
			tmpBeta -= itr->value * ((Matrix3::Identity - m_Constraints[col]) * m_X[col]); //Vel;
			*/

			//tmp = (Matrix3::Identity - m_Constraints[col]) * m_X[col];
			tmpVec.x = (1.0f - constraints[col]._11) * xvec[col].x;
			tmpVec.y = constraints[col]._21 * xvec[col].x;
			tmpVec.z = constraints[col]._31 * xvec[col].x;

			tmpVec.x += constraints[col]._12 * xvec[col].y;
			tmpVec.y += (1.0f - constraints[col]._22) * xvec[col].y;
			tmpVec.z += constraints[col]._32 * xvec[col].y;

			tmpVec.x += constraints[col]._13 * xvec[col].z;
			tmpVec.y += constraints[col]._23 * xvec[col].z;
			tmpVec.z += (1.0f - constraints[col]._33) * xvec[col].z;

			InplaceMatrix3MultVector3Subtract(&tmpRes, itr->value, xvec[col]);
			InplaceMatrix3MultVector3Subtract(&tmpBeta, itr->value, tmpVec);

		}
		/*m_Residual[row] = m_Constraints[row] * tmpRes;
		m_Previous[row] = m_Constraints[row] * m_PreCondition[row] * m_Residual[row];*/

		InplaceMatrix3MultVector3(&out_residual[row], constraints[row], tmpRes);


		InplaceMatrix3MultMatrix3(&tmp, constraints[row], inv_precondition[row]);
		InplaceMatrix3MultVector3(&out_previous[row], tmp, out_residual[row]);

		//m_Update[row]   = m_Previous[row];

		beta += tmpBeta.Dot(inv_precondition[row] * tmpBeta);
		r0z0 += Vector3::Dot(out_residual[row], out_previous[row]);
	}

	out_beta = beta;
	out_r0z0 = r0z0;
}

template<class T>
float SparseRowMatrix<T>::SolveAMultU(std::vector<Vector3>& out, const std::vector<Matrix3>& constraints, const std::vector<Vector3>& u)
{
	float accum = 0.0f;
	int len = (int)m_Rows.size();
#pragma omp parallel for reduction(+:accum)// reduction(+:r0z0)
	for (int row = 0; row < len; ++row)
	{
		Vector3 temp(0.0f, 0.0f, 0.0f);

		auto itr = m_Rows[row].begin(), end = m_Rows[row].end();
		for (; itr != end; itr++)
		{
			InplaceMatrix3MultVector3Additve(&temp, itr->value, u[itr->column]);
		}
		InplaceMatrix3MultVector3(&out[row], constraints[row], temp);
		accum += Vector3::Dot(u[row], out[row]);

		/*float temp_x = 0.0f, temp_y = 0.0f, temp_z = 0.0f;
		float outr[3];
		auto& arr = m_Rows[row];
		int rlen = (int)arr.size();
		for (int i = 0; i < rlen; ++i)
		{
		Matrix3& mtx = arr[i].value;
		const Vector3& uvec = u[arr[i].column];

		temp_x += mtx._11 * uvec.x;
		temp_y += mtx._21 * uvec.x;
		temp_z += mtx._31 * uvec.x;

		temp_x += mtx._12 * uvec.y;
		temp_y += mtx._22 * uvec.y;
		temp_z += mtx._32 * uvec.y;

		temp_x += mtx._13 * uvec.z;
		temp_y += mtx._23 * uvec.z;
		temp_z += mtx._33 * uvec.z;
		}


		const Matrix3& mtx = constraints[row];

		outr[0] = mtx._11 * temp_x + mtx._12 * temp_y + mtx._13 * temp_z;
		outr[1] = mtx._21 * temp_x + mtx._22 * temp_y + mtx._23 * temp_z;
		outr[2] = mtx._31 * temp_x + mtx._32 * temp_y + mtx._33 * temp_z;

		memcpy(&out[row], outr, 3 * sizeof(float));
		const Vector3& urow = u[row];
		accum += urow.x * outr[0] + urow.y * outr[1] + urow.z * outr[2];*/
	}
	return accum;
}