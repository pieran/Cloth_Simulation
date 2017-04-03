#include "mpcg.h"

const float tiny = 1e-030f;       // TODO: Should be user controllable

template<class T>
MPCG<T>::MPCG()
{
	m_NumTotal = 0;
	m_Tolerence = 1e-3f;// 0.0001f;
	m_MaxIterations = 1000;
	m_EstimatedError = 0.0f;
	m_Iterations = 0;
}

template<class T>
MPCG<T>::~MPCG()
{
}

template<class T>
void MPCG<T>::AllocateMemory(uint num_total)
{
	m_NumTotal = num_total;

	m_Residual.resize(m_NumTotal);
	m_Previous.resize(m_NumTotal);
	m_Update.resize(m_NumTotal);
	m_UpdateA.resize(m_NumTotal);

	m_PreCondition.resize(m_NumTotal);
	m_Constraints.resize(m_NumTotal);

	m_B.resize(m_NumTotal);
	m_X.resize(m_NumTotal);
	m_ValueAccum1.resize(m_NumTotal);

	memset(&m_B[0], 0, m_NumTotal * sizeof(Vector3));
	memset(&m_X[0], 0, m_NumTotal * sizeof(Vector3));
	for (uint i = 0; i < m_NumTotal; ++i)
	{
		m_Constraints[i] = Matrix3::Identity;
		m_PreCondition[i] = Matrix3::Identity;
	}
}

template<class T>
void MPCG<T>::ResetMemory()
{
	memset(&m_B[0].x, 0, m_NumTotal * sizeof(Vector3));
}

template<class T>
void MPCG<T>::Solve()
{
	memset(&m_X[0].x, 0, m_NumTotal * sizeof(Vector3));
	Solve_Algorithm();

	m_ProfilingAverageIterations_Sum += m_Iterations;
	m_ProfilingAverageIterations_No++;
}

template<class T>
void MPCG<T>::SolveWithGuess(const std::vector<Vector3>& guess)
{
	memset(&m_X[0].x, 0, m_NumTotal * sizeof(Vector3));
	memcpy(&m_X[0].x, &guess[0].x, guess.size() * sizeof(Vector3));
	Solve_Algorithm();

	m_ProfilingAverageIterations_Sum += m_Iterations;
	m_ProfilingAverageIterations_No++;
}

template<class T>
void MPCG<T>::SolveWithGuess(const Vector3* guess)
{
	memcpy(&m_X[0].x, guess, m_NumTotal * sizeof(Vector3));
	Solve_Algorithm();

	m_ProfilingAverageIterations_Sum += m_Iterations;
	m_ProfilingAverageIterations_No++;
}

template<class T>
void MPCG<T>::SolveWithPreviousResult()
{
	Solve_Algorithm();

	m_ProfilingAverageIterations_Sum += m_Iterations;
	m_ProfilingAverageIterations_No++;
}

template<class T>
void MPCG<T>::Solve_Algorithm()
{
	float r0z0, d2;
	float beta = 0.0f;
	float delta = 0.0f;
	/*
	for (int i = 0; i < int(m_NumTotal); ++i)
	{
	m_PhyxelsVelChange[i] = m_Constraints[i] * m_PhyxelsVelChange[i] + (Matrix3::Identity - m_Constraints[i]) * m_PhyxelsVelChange[i]; //Vel;
	}
	*/

	//memset(&m_PhyxelsVelChange[0].x, 0, m_NumTotal * sizeof(Vector3));

	m_ProfilingInitialization.BeginTiming();
	r0z0 = 0.0f;

	/*#pragma omp parallel for reduction(+:beta) reduction(+:r0z0)
	for (int row = 0; row < m_NumTotal; ++row) {
	Vector3 tmpRes = m_B[row];
	Vector3 tmpBeta = m_B[row];
	Vector3 tmpVec;

	const auto& a_row = m_A.GetRow(row);
	auto itr = a_row.begin(), end = a_row.end();
	for (; itr != end; itr++)
	{
	uint col = itr->column;
	//tmpRes -= itr->value * m_X[col];
	//tmpBeta -= itr->value * ((Matrix3::Identity - m_Constraints[col]) * m_X[col]); //Vel;

	//tmp = (Matrix3::Identity - m_Constraints[col]) * m_X[col];
	tmpVec.x = (1.0f - m_Constraints[col]._11) * m_X[col].x;
	tmpVec.y = m_Constraints[col]._21 * m_X[col].x;
	tmpVec.z = m_Constraints[col]._31 * m_X[col].x;

	tmpVec.x += m_Constraints[col]._12 * m_X[col].y;
	tmpVec.y += (1.0f - m_Constraints[col]._22) * m_X[col].y;
	tmpVec.z += m_Constraints[col]._32 * m_X[col].y;

	tmpVec.x += m_Constraints[col]._13 * m_X[col].z;
	tmpVec.y += m_Constraints[col]._23 * m_X[col].z;
	tmpVec.z += (1.0f - m_Constraints[col]._33) * m_X[col].z;

	InplaceMatrix3MultVector3Subtract(&tmpRes, itr->value, m_X[col]);
	InplaceMatrix3MultVector3Subtract(&tmpBeta, itr->value, tmpVec);

	}
	//m_Residual[row] = m_Constraints[row] * tmpRes;
	//m_Previous[row] = m_Constraints[row] * m_PreCondition[row] * m_Residual[row];

	InplaceMatrix3MultVector3(&m_Residual[row], m_Constraints[row], tmpRes);

	Matrix3 tmp;
	InplaceMatrix3MultMatrix3(&tmp, m_Constraints[row], m_PreCondition[row]);
	InplaceMatrix3MultVector3(&m_Previous[row], tmp, m_Residual[row]);

	//m_Update[row]   = m_Previous[row];

	beta += tmpBeta.Dot(m_PreCondition[row] * tmpBeta);
	r0z0 += Vector3::Dot(m_Residual[row], m_Previous[row]);
	}*/

	m_A.SolveAMultX(m_Residual, m_Previous, r0z0, beta, m_Constraints, m_PreCondition, m_B, m_X);
	memcpy(&m_Update[0], &m_Previous[0], m_NumTotal * sizeof(Vector3));
	m_ProfilingInitialization.EndTimingAdditive();


	float tolSqBeta = m_Tolerence * m_Tolerence * beta;

	for (m_Iterations = 0; m_Iterations < m_MaxIterations; ++m_Iterations)
	{
		d2 = 0.0f;

		// alpha = Dot(r0, z0) / Dot(p0, Ap0)
		m_ProfilingUpper.BeginTiming();
		/*#pragma omp parallel for reduction(+:d2)// reduction(+:r0z0)
		for (int row = 0; row < m_NumTotal; ++row)
		{
		Vector3 temp(0.0f, 0.0f, 0.0f);
		const auto& a_row = m_A.GetRow(row);
		auto itr = a_row.begin(), end = a_row.end();
		for (; itr != end; itr++)
		{
		InplaceMatrix3MultVector3Additve(&temp, itr->value, m_Update[itr->column]);
		}
		InplaceMatrix3MultVector3(&m_UpdateA[row], m_Constraints[row], temp);

		//r0z0 += Vector3::Dot(m_Residual[row], m_Previous[row]);
		d2 += Vector3::Dot(m_Update[row], m_UpdateA[row]);
		}*/
		d2 = m_A.SolveAMultU(m_UpdateA, m_Constraints, m_Update);
		m_ProfilingUpper.EndTimingAdditive();


		if (fabs(d2) < tiny) d2 = tiny;
		float alpha = r0z0 / d2;

		// x1 = x0 + p0 * aplha
		// r1 = r0 - Ap0 * alpha

		float errorSq = 0.0f;
#pragma omp parallel for reduction(+:errorSq)
		for (int row = 0; row < m_NumTotal; ++row)
		{
			m_X[row] += m_Update[row] * alpha;
			m_Residual[row] -= m_UpdateA[row] * alpha;

			errorSq += Vector3::Dot(m_Residual[row], m_PreCondition[row] * m_Residual[row]);
		}
		m_EstimatedError = errorSq;

		// if (r1 is small) exit;
		if (m_EstimatedError < tolSqBeta)
		{
			m_EstimatedError /= beta;
			m_EstimatedError = sqrt(m_EstimatedError);
			return;
		}

		// z1 = Minv . r1
		d2 = 0.0f;
		m_ProfilingLower.BeginTiming();
#pragma omp parallel for reduction(+:d2)
		for (int row = 0; row < m_NumTotal; ++row)
		{
			//m_Previous[row] = m_PreCondition[row] * m_Residual[row];
			InplaceMatrix3MultVector3(&m_Previous[row], m_PreCondition[row], m_Residual[row]);
			d2 += Vector3::Dot(m_Previous[row], m_Residual[row]);
		}


		// change = Dot(z1, r1) / Dot(z0, r0)
		if (fabs(r0z0) < tiny) r0z0 = tiny;
		float change = d2 / r0z0;
		r0z0 = d2;

		// p1 = z1 + p0 * beta
#pragma omp parallel for
		for (int row = 0; row < m_NumTotal; ++row)
		{
			//m_Update[row] = m_Constraints[row] * (m_Previous[row] + m_Update[row] * change);
			Vector3 temp = (m_Previous[row] + m_Update[row] * change);
			InplaceMatrix3MultVector3(&m_Update[row], m_Constraints[row], temp);
		}
		m_ProfilingLower.EndTimingAdditive();
	}

	m_EstimatedError /= beta;
	m_EstimatedError = sqrt(m_EstimatedError);
}

template<class T>
void MPCG<T>::ResetProfilingData()
{
	m_ProfilingLower.ResetTotalMs();
	m_ProfilingUpper.ResetTotalMs();
	m_ProfilingInitialization.ResetTotalMs();

	m_ProfilingAverageIterations_Sum = 0;
	m_ProfilingAverageIterations_No = 0;
}
