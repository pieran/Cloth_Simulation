#include "mpcg.h"

const float tiny = 1e-030f;       // TODO: Should be user controllable


MPCG::MPCG()
{
	m_NumTotal = 0;
	m_Tolerence = 1e-4f;// 0.0001f;
	m_MaxIterations = 1000;
	m_Iterations = 0;
	m_EstimatedError = 0.0f;
}

MPCG::~MPCG()
{
}

void MPCG::AllocateMemory(uint num_total)
{
	m_NumTotal = num_total;

	m_A.resize(m_NumTotal);
	m_B.resize(m_NumTotal);
	m_X.resize(m_NumTotal);
	m_PreCondition.resize(m_NumTotal);

	m_Residual.resize(m_NumTotal);
	m_AuxVec.resize(m_NumTotal);
	m_SearhVec.resize(m_NumTotal);
	m_PreCondition.resize(m_NumTotal);

	memset(&m_B[0], 0, m_NumTotal * sizeof(float));
	memset(&m_X[0], 0, m_NumTotal * sizeof(float));
	for (uint i = 0; i < m_NumTotal; ++i)
	{
		m_PreCondition[i] = 1.0f;
	}
}

void MPCG::ResetMemory()
{
	memset(&m_B[0], 0, m_NumTotal * sizeof(float));
}

void MPCG::Solve()
{
	BuildPreconditioner();
	SolveAlgorithm();
}


void MPCG::BuildPreconditioner()
{
	for (int i = 0; i < (int)m_NumTotal; ++i)
	{
		float aii = m_A[i][i];
		m_PreCondition[i] = (fabs(aii) < 1E-6) ? 1.0f : 1.0f / m_A[i][i];
	}
}

void MPCG::SolveAlgorithm()
{
	const float tolerence_squared = m_Tolerence * m_Tolerence;
	int inumtotal = (int)m_NumTotal;

	memset(&m_X[0], 0, m_NumTotal * sizeof(float));

	//Set search vec s = precondition(r)
#pragma omp parallel for
	for (int i = 0; i < inumtotal; ++i)
	{
		m_SearhVec[i] = m_PreCondition[i] * m_B[i];
		m_Residual[i] = m_B[i];
	}


	float rsold = 0.0f;
#pragma omp parallel for reduction(+:rsold)
	for (int i = 0; i < inumtotal; ++i)
	{
		rsold += m_SearhVec[i] * m_Residual[i];
	}


	for (m_Iterations = 0; m_Iterations < m_MaxIterations; ++m_Iterations)
	{
		float rs = 0.0f;

#pragma omp parallel for reduction(+:rs)
		for (int i = 0; i < inumtotal; ++i)
		{
			float val = 0.0f;
			auto& rows = m_A[i];
			for (auto itr = rows.begin(), end = rows.end(); itr != end; ++itr)
			{
				val += itr->second * m_SearhVec[itr->first];
			}

			rs += val * m_SearhVec[i];
			m_AuxVec[i] = val;			
		}
		
		float alpha = rsold / rs;
		
		float error_squared = 0.0f;
#pragma omp parallel for reduction(+: error_squared)
		for (int i = 0; i < inumtotal; ++i)
		{
			m_X[i] += alpha * m_SearhVec[i];
			m_Residual[i] -= alpha * m_AuxVec[i];

			error_squared += m_Residual[i] * m_Residual[i];
		}

		if (error_squared <= tolerence_squared)
		{
			return;
		}

		float rsnew = 0.0f;
#pragma omp parallel for reduction(+: rsnew)
		for (int i = 0; i < inumtotal; ++i)
		{
			m_AuxVec[i] = m_PreCondition[i] * m_Residual[i];
			rsnew += m_AuxVec[i] * m_Residual[i];
		}

		float beta = rsnew / rsold;

#pragma omp parallel for
		for (int i = 0; i < inumtotal; ++i)
		{
			m_SearhVec[i] = m_AuxVec[i] + beta * m_SearhVec[i];
		}

		rsold = rsnew;
	}
}
