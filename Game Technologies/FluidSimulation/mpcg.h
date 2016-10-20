
#pragma once
#include <nclcore\common.h>
#include <vector>
#include <map>


//typedef std::vector<std::map<uint, Matrix3>> MatrixMap;

//Modified Preconditioned Conjugate Gradient
class MPCG
{
public:
	MPCG();
	~MPCG();

	void AllocateMemory(uint num_total);	//A-Matrix Memory MUST be allocated seperately!!!
	void ResetMemory();



	void SetMaxIterations(uint max) { m_MaxIterations = max; }
	void SetTolerance(float tol) { m_Tolerence = tol; }

	uint GetMaxIterations() const { return m_MaxIterations; }
	float GetTolerance() const { return m_Tolerence; }


	std::vector<float>					m_X;
	std::vector<std::map<uint, float>>	m_A;
	std::vector<float>					m_B;

	void Solve();

protected:

	void BuildPreconditioner();
	void SolveAlgorithm();

protected:
	uint					m_NumTotal;
	uint					m_MaxIterations;
	uint					m_Iterations;
	float					m_Tolerence;
	float					m_EstimatedError;


	std::vector<float>					m_PreCondition; 

	std::vector<float>	m_Residual;
	std::vector<float>	m_AuxVec;
	std::vector<float>	m_SearhVec;
};
