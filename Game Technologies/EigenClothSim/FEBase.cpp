#include "FEBase.h"

FEBase::FEBase() : ClothBase()
{
//	m_UseRungeKutta4Integration = false;
//	m_NumPhyxels = 0;
	m_SolverAllocated = 0;
}

FEBase::~FEBase()
{
}

uint FEBase::AddPoint(const uint& triIdx, const Vector3& barycentric_coords)
{
	uint nIdx = ClothBase::AddPoint(triIdx, barycentric_coords);

	uint idxA = m_TriIndicies[triIdx];
	uint idxB = m_TriIndicies[triIdx + 1];
	uint idxC = m_TriIndicies[triIdx + 2];

	m_PhyxelPosInitial[nIdx] = barycentric_interpolate<Vector3>(m_PhyxelPosInitial[idxA], m_PhyxelPosInitial[idxB], m_PhyxelPosInitial[idxC], barycentric_coords);
	return nIdx;
}

void FEBase::UpdateTriangleMaterial(const uint& triIdx)
{
	
}

void FEBase::InitializeClothDesign(const ClothDesignEntity<3>* design)
{
	ClothBase::InitializeClothDesign(design);

	uint i;

	//Phyxel Memory
	m_PhyxelPosInitial.resize(m_NumAllocatedPhyxels);		
	for (i = 0; i < m_NumPhyxels; ++i)
	{
		const Vertex& v = design->Vertices()[i];
		m_PhyxelPosInitial[i] = Vector3(v.initial_position.x(), v.initial_position.y(), v.initial_position.z());
	}

	//Tri Memory
	m_NumTriMappings = m_NumTris;
	m_TriB.resize(m_NumAllocatedTris * 3);
	m_TriAii.resize(m_NumAllocatedTris * 3);
	m_TriAij.resize(m_NumAllocatedTris * 3);
	BuildPhyxelTriMappings();


	//Solver
	ReAllocateSolverMemory();
	m_Solver.ResetMemory();
	m_Solver.m_A.zero_memory();
}

void FEBase::ReAllocateSolverMemory()
{
	if (m_SolverAllocated != m_NumPhyxels)
	{
		m_SolverAllocated = m_NumPhyxels;
		m_Solver.AllocateMemory(m_NumPhyxels);
		m_Solver.m_A.resize(m_NumPhyxels);
		for (int i = 0; i < (int)m_NumPhyxels; ++i)
		{
			m_Solver.m_A(i, i) = Matrix3();

			if (m_PhyxelIsStatic[i])
				m_Solver.m_Constraints[i].ToZero();
			else
				m_Solver.m_Constraints[i].ToIdentity();

			const auto& tris = m_PhyxelTriLookupIJ[i];
			for (auto itr = tris.begin(), end = tris.end(); itr != end; ++itr)
			{
				uint idxB = itr->first;

				m_Solver.m_A(i, idxB) = Matrix3();
				m_Solver.m_A(idxB, i) = Matrix3();
			}
		}
	}

	if (m_NumTriMappings != m_NumTris)
	{
		m_NumTriMappings = m_NumTris;
		BuildPhyxelTriMappings();
	}
}

void FEBase::ComputeVelocity(float subtimestep, const Vector3* in_clothpos, const Vector3* in_clothvel, Vector3* out_clothvel)
{
	m_Profiler.BeginTiming(PROFILERID_CLOTH_BUILDROTATIONANDNORMALS);
	GenRotations(in_clothpos);
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_BUILDROTATIONANDNORMALS);

	//TODO: Loop until valid timestep achieved
	m_Profiler.BeginTiming(PROFILERID_CLOTH_SOLVER);
	ReAllocateSolverMemory();
	m_Solver.ResetMemory();
	m_Solver.m_A.zero_memory();
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_SOLVER);


	m_Profiler.BeginTiming(PROFILERID_CLOTH_COMPUTEFORCES);
	memset(&m_TriB[0].x, 0, m_NumTris * 3 * sizeof(Vector3));
	GenStressStrainForces(subtimestep, in_clothpos, in_clothvel);
	ConstructGlobalMatrix(subtimestep, in_clothpos, in_clothvel);
	m_Profiler.EndTiming(PROFILERID_CLOTH_COMPUTEFORCES);


	m_Profiler.BeginTiming(PROFILERID_CLOTH_SOLVER);
	m_Solver.SolveWithGuess(in_clothvel);
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_SOLVER);

	memcpy(out_clothvel, &m_Solver.m_X[0], m_NumPhyxels * sizeof(Vector3));
}

void FEBase::UpdateStaticPhyxelConstraints()
{
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumPhyxels; ++i)
	{
		if (m_PhyxelIsStatic[i])
		{
			m_Solver.m_Constraints[i] = Matrix3::ZeroMatrix;
		}
		else
		{
			m_Solver.m_Constraints[i] = Matrix3::Identity;
		}
	}
}

void FEBase::BuildPhyxelTriMappings()
{
	m_PhyxelTriLookup.clear();
	m_PhyxelTriLookupIJ.clear();

	m_PhyxelTriLookup.resize(m_NumAllocatedPhyxels);
	m_PhyxelTriLookupIJ.resize(m_NumAllocatedPhyxels);

	for (unsigned int i = 0; i < m_NumTris; ++i)
	{
		unsigned int i3 = i * 3;
		uint idxA = m_TriIndicies[i3];
		uint idxB = m_TriIndicies[i3 + 1];
		uint idxC = m_TriIndicies[i3 + 2];

		m_PhyxelTriLookup[idxA].push_back(i3 + 0);
		m_PhyxelTriLookup[idxB].push_back(i3 + 1);
		m_PhyxelTriLookup[idxC].push_back(i3 + 2);

		bool AltB = idxA < idxB;
		bool BltC = idxB < idxC;
		bool CltA = idxC < idxA;

		//AB
		m_PhyxelTriLookupIJ[idxA][idxB].push_back(i3 + 0);
		m_PhyxelTriLookupIJ[idxB][idxA].push_back(i3 + 0);

		//BC
		m_PhyxelTriLookupIJ[idxB][idxC].push_back(i3 + 1);
		m_PhyxelTriLookupIJ[idxC][idxB].push_back(i3 + 1);

		//CA
		m_PhyxelTriLookupIJ[idxC][idxA].push_back(i3 + 2);
		m_PhyxelTriLookupIJ[idxA][idxC].push_back(i3 + 2);
	}
}

void FEBase::ConstructGlobalMatrix(float subtimestep, const Vector3* in_clothpos, const Vector3* in_phyxelVel)
{
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumPhyxels; ++i)
	{
		float mass = m_PhyxelMass[i];

		Vector3 BVec = Vector3(0.0f, 0.0, 0.0f);
		Matrix3 Aii = Matrix3(mass, 0.0f, 0.0f,
			0.0f, mass, 0.0f,
			0.0f, 0.0f, mass);

		const auto& tris = m_PhyxelTriLookup[i];
		for (uint itr = 0, end = tris.size(); itr < end; ++itr)
		{
			BVec += m_TriB[tris[itr]];
			Aii += m_TriAii[tris[itr]];
		}

		float invmass = 1.0f / mass;
		m_RenderValue[i] = BVec.Length() * invmass;


		BVec += m_PhyxelForces[i];
		BVec += m_PhyxelExtForces[i];
		BVec *= subtimestep;
		BVec += in_phyxelVel[i] * mass;

		m_Solver.m_A(i, i) = Aii;
		m_Solver.m_PreCondition[i] = Matrix3::Inverse(Aii);
		m_Solver.m_B[i] = BVec;


		uint idxA = i;
		const auto& edges = m_PhyxelTriLookupIJ[i];
		for (auto itr = edges.begin(), end = edges.end(); itr != end; ++itr)
		{
			uint idxB = itr->first;
			Matrix3 Aij;
			for (uint k = 0, k_end = itr->second.size(); k < k_end; ++k)
			{
				Aij += m_TriAij[itr->second[k]];
			}

			m_Solver.m_A(idxA, idxB) = Aij;
		}
	}
}


void FEBase::SetSolverConstraint(int idx, const Matrix3& constraint)
{
	m_Solver.m_Constraints[idx] = constraint;
}