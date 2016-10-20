
#pragma once
#include "ProfilingTimer.h"
#include "SparseRowMatrix.h"
#include "PArray.h"
#include "common.h"
#include <nclgl\Vector3.h>
#include <nclgl\Matrix3.h>
#include <vector>
#include <map>


//typedef std::vector<std::map<uint, Matrix3>> MatrixMap;

//Modified Preconditioned Conjugate Gradient
template<class T>
class MPCG
{
public:
	MPCG<T>();
	~MPCG<T>();

	inline void AllocateMemory(uint num_total);	//A-Matrix Memory MUST be allocated seperately!!!
	inline void ResetMemory();



	inline void SetMaxIterations(uint max)	{ m_MaxIterations = max; }
	inline void SetTolerance(float tol)		{ m_Tolerence = tol; }

	inline uint GetMaxIterations() const	{ return m_MaxIterations; }
	inline float GetTolerance() const		{ return m_Tolerence; }



	T							m_A;
	std::vector<Vector3>		m_B;
	std::vector<Vector3>		m_X;
	std::vector<Matrix3>		m_Constraints;
	std::vector<Matrix3>		m_PreCondition;	


	void Solve();
	void SolveWithGuess(const std::vector<Vector3>& guess);
	void SolveWithGuess(const Vector3* guess);
	void SolveWithPreviousResult();	//SolveWithGuess(this->m_X)

	inline void ResetProfilingData();

	ProfilingTimer m_ProfilingInitialization;
	ProfilingTimer m_ProfilingUpper;
	ProfilingTimer m_ProfilingLower;
	float		   GetAverageIterations() { return ((float)m_ProfilingAverageIterations_Sum / (float)m_ProfilingAverageIterations_No); }

protected:
	uint		   m_ProfilingAverageIterations_Sum;
	uint		   m_ProfilingAverageIterations_No;

	void Solve_Algorithm();

protected:
	uint					m_MaxIterations;
	uint					m_Iterations;
	float					m_Tolerence;
	float					m_EstimatedError;
	uint					m_NumTotal;

	std::vector<Vector3>	m_Residual;
	std::vector<Vector3>	m_Previous;
	std::vector<Vector3>	m_Update;
	std::vector<Vector3>	m_UpdateA;
	std::vector<float>		m_ValueAccum1;
};

#include "mpcg.inl"