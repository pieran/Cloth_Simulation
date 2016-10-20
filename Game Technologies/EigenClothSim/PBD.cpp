#include "PBD.h"
#include <list>

PBD::PBD() : ClothBase()
{
	m_SolverIterations = 50;
}

PBD::~PBD()
{

}

bool HACKY_IsContainedInDistanceList(uint idxA, uint idxB, const std::list<PBDDistanceConstraint>& list)
{
	for (auto itr = list.begin(), end = list.end(); itr != end; itr++)
	{
		if (itr->p1 == idxA && itr->p2 == idxB)
			return true;
	}
	return false;
}

void PBD::InitializeClothDesign(const ClothDesignEntity<3>* design)
{
	ClothBase::InitializeClothDesign(design);
	InitializePBDConstraints();
}

void PBD::InitializePBDConstraints()
{
	m_PhyxelPosNew.resize(m_NumPhyxels);
	m_PhyxelInvMass.resize(m_NumPhyxels);

#pragma omp parallel for
	for (int i = 0; i < (int)m_NumPhyxels; ++i)
	{
		if (m_PhyxelIsStatic[i])
			m_PhyxelInvMass[i] = 0.0f;
		else
			m_PhyxelInvMass[i] = 1.0f / m_PhyxelMass[i];
	}

	//For each quad, create:
	//	- 4 perimeter distance constraints (REQUIRES CHECKING FOR DUPLICATES!!!)
	//	- 2 shear distance constraints
	//	- 2 bending constraints
	std::list<PBDDistanceConstraint> d_constraints;
	std::list<PBDBendingConstraint>  b_constraints;

	uint numBendingConstraints = m_NumQuads * 2;
	uint numDistanceConstraints = 0;

	PBDDistanceConstraint dc;
	PBDBendingConstraint bc;
	for (int i = 0; i < (int)m_NumQuads; ++i)
	{
		int i4 = i * 4;
		uint idxA = m_QuadIndicies[i4];
		uint idxB = m_QuadIndicies[i4 + 1];
		uint idxC = m_QuadIndicies[i4 + 2];
		uint idxD = m_QuadIndicies[i4 + 3];
		//Quad Layout
		//  A -- B
		//  |    |
		//  D -- C

		//Distance TL-BL
		if (!HACKY_IsContainedInDistanceList(idxA, idxD, d_constraints))
		{
			if (BuildDistanceConstraint(idxA, idxD, m_MaterialStretchWarp, dc))
				d_constraints.push_back(dc);
		}

		//Distance TR-BR
		if (!HACKY_IsContainedInDistanceList(idxB, idxC, d_constraints))
		{
			if (BuildDistanceConstraint(idxB, idxC, m_MaterialStretchWarp, dc))
				d_constraints.push_back(dc);
		}

		//Distance TL-TR
		if (!HACKY_IsContainedInDistanceList(idxA, idxB, d_constraints))
		{
			if (BuildDistanceConstraint(idxA, idxB, m_MaterialStretchWeft, dc))
				d_constraints.push_back(dc);
		}

		//Distance BL-BR
		if (!HACKY_IsContainedInDistanceList(idxD, idxC, d_constraints))
		{
			if (BuildDistanceConstraint(idxD, idxC, m_MaterialStretchWeft, dc))
				d_constraints.push_back(dc);
		}

		//Shear TL-BR
		if (BuildDistanceConstraint(idxA, idxC, m_MaterialShear, dc))
			d_constraints.push_back(dc);

		//Shear TR-BL
		if (BuildDistanceConstraint(idxB, idxC, m_MaterialShear, dc))
			d_constraints.push_back(dc);


		//Bending A-C-B
		if (BuildBendingConstraint(idxA, idxC, idxB, m_MaterialBend, bc))
			b_constraints.push_back(bc);

		//Bending A-D-C
		if (BuildBendingConstraint(idxA, idxD, idxC, m_MaterialBend, bc))
			b_constraints.push_back(bc);
	}


	int i = 0;
	m_ConstraintsDistance.resize(d_constraints.size());
	for (auto itr = d_constraints.begin(), end = d_constraints.end(); itr != end; ++itr)
	{
		m_ConstraintsDistance[i++] = *itr;
	}

	i = 0;
	m_ConstraintsBending.resize(b_constraints.size());
	for (auto itr = b_constraints.begin(), end = b_constraints.end(); itr != end; ++itr)
	{
		m_ConstraintsBending[i++] = *itr;
	}

	uint large_size_d = d_constraints.size() * 2;
	uint large_size_b = b_constraints.size() * 3;
	m_LargeArrayParallisation.resize(max(large_size_d, large_size_b));

	BuildLookupArrays();
}

void PBD::BuildLookupArrays()
{
	m_LookupsDistance.clear();
	m_LookupsBending.clear();

	m_LookupsDistance.resize(m_NumPhyxels);
	m_LookupsBending.resize(m_NumPhyxels);

	int i, ilen;
	for (i = 0, ilen = (int)m_ConstraintsDistance.size(); i < ilen; ++i) {
		PBDDistanceConstraint& c = m_ConstraintsDistance[i];
		m_LookupsDistance[c.p1].push_back(i * 2);
		m_LookupsDistance[c.p2].push_back(i * 2 + 1);
	}

	for (i = 0, ilen = (int)m_ConstraintsBending.size(); i < ilen; ++i) {
		PBDBendingConstraint& c = m_ConstraintsBending[i];
		m_LookupsBending[c.p1].push_back(i * 3);
		m_LookupsBending[c.p2].push_back(i * 3 + 1);
		m_LookupsBending[c.p3].push_back(i * 3 + 2);
	}
}

void PBD::ComputeVelocity(float subtimestep, const Vector3* in_clothpos, const Vector3* in_clothvel, Vector3* out_clothvel)
{
	float dampFactor = 1.0f - m_MaterialDamp;

#pragma omp parallel for
	for (int i = 0; i < m_NumPhyxels; ++i)
	{
		out_clothvel[i] = in_clothvel[i] * dampFactor + (m_PhyxelForces[i] + m_PhyxelExtForces[i]) * subtimestep;
		out_clothvel[i] = in_clothvel[i] * dampFactor + (m_PhyxelForces[i] + m_PhyxelExtForces[i]) * subtimestep * m_PhyxelInvMass[i];

		m_PhyxelPosNew[i] = in_clothpos[i] + out_clothvel[i] * subtimestep;
	}

	float weighting = 1.0f;
	for (uint si = 0; si < m_SolverIterations; ++si) {
		weighting = 1.0f - (float)(si * si) / (float)(m_SolverIterations * m_SolverIterations);
		SolveDistanceConstraints(weighting);
		SolveBendingConstraints(weighting);
	}

	float inv_timestep = 1.0f / subtimestep;
#pragma omp parallel for
	for (int i = 0; i < m_NumPhyxels; ++i)
	{
		out_clothvel[i] = (m_PhyxelPosNew[i] - in_clothpos[i]) * inv_timestep;
	}
}

bool PBD::BuildDistanceConstraint(uint idxA, uint idxB, float k, PBDDistanceConstraint& out_constraint)
{
	float invmass = m_PhyxelInvMass[idxA] + m_PhyxelInvMass[idxB];
	if (invmass < FLT_EPSILON)
		return false;

	out_constraint.p1 = idxA;
	out_constraint.p2 = idxB;
	out_constraint.k = k;
	out_constraint.k_prime = (1.0f - pow((1.0f - k), 1.0f / m_SolverIterations));  //1.0f-pow((1.0f-c.k), 1.0f/ns);

	if (out_constraint.k_prime > 1.0f)
		out_constraint.k_prime = 1.0f;

	Vector3 deltaP = m_ClothStatesPos[idxA] - m_ClothStatesPos[idxB];
	out_constraint.rest_length = deltaP.Length();

	return true;
}

bool PBD::BuildBendingConstraint(uint idxA, uint idxB, uint idxC, float k, PBDBendingConstraint& out_constraint)
{
	out_constraint.p1 = idxA;
	out_constraint.p2 = idxB;
	out_constraint.p3 = idxC;

	out_constraint.w = m_PhyxelInvMass[idxA] + m_PhyxelInvMass[idxB] + m_PhyxelInvMass[idxC];

	if (out_constraint.w < FLT_EPSILON)
		return false;

	out_constraint.w = 1.0f / out_constraint.w;

	out_constraint.k = k;
	out_constraint.k_prime = 1.0f - pow((1.0f - k), 1.0f / m_SolverIterations);  //1.0f-pow((1.0f-c.k), 1.0f/ns);

	if (out_constraint.k_prime > 1.0f)
		out_constraint.k_prime = 1.0f;

	Vector3 center = (m_ClothStatesPos[idxA] + m_ClothStatesPos[idxB] + m_ClothStatesPos[idxC]) * 0.3333333f;
	out_constraint.rest_length = (m_ClothStatesPos[idxC] - center).Length();

	return true;
}

void PBD::SolveDistanceConstraints(float weighting)
{
	Vector3 dir, dP;
	float len;

	int i, ilen = m_ConstraintsDistance.size();
#pragma omp parallel for private(i, len, dir, dP)
	for (i = 0; i < ilen; ++i) {
		PBDDistanceConstraint& c = m_ConstraintsDistance[i];
		dir = m_PhyxelPosNew[c.p1] - m_PhyxelPosNew[c.p2];

		float invmassA = m_PhyxelIsStatic[c.p1] ? 0.0f : m_PhyxelInvMass[c.p1];
		float invmassB = m_PhyxelIsStatic[c.p2] ? 0.0f : m_PhyxelInvMass[c.p2];
		float w = invmassA + invmassB;

		if (w > 0.0f)
		{
			len = dir.Length();
			//if (len <= 0.00001f)
			//	return;

			dP = dir / len * (c.k_prime * (len - c.rest_length));
			dP = dP / w;

			m_LargeArrayParallisation[i * 2] = -dP * invmassA;
			m_LargeArrayParallisation[i * 2 + 1] = dP * invmassB;
			//m_PhyxelPosNew[c.p1] -= dP * m_PhyxelInvMass[c.p1];
			//m_PhyxelPosNew[c.p2] += dP * m_PhyxelInvMass[c.p2];
		}
		else
		{
			m_LargeArrayParallisation[i * 2] = Vector3(0.0f, 0.0f, 0.0f);
			m_LargeArrayParallisation[i * 2 + 1] = Vector3(0.0f, 0.0f, 0.0f);
		}
	}

	ilen = (int)m_NumPhyxels;

	Vector3 sum;
#pragma omp parallel for private(i, sum)
	for (i = 0; i < ilen; ++i)
	{
		auto& arr = m_LookupsDistance[i];
		sum = Vector3(0, 0, 0);
		for (uint idx : arr)
		{
			sum += m_LargeArrayParallisation[idx];
		}
		sum *= weighting;
		m_PhyxelPosNew[i] += sum;
	}
}

void PBD::SolveBendingConstraints(float weighting)
{
	Vector3 center, dir_center, dir_force;
	float dist_center, diff;
	int i, ilen;
	ilen = (int)m_ConstraintsBending.size();

#pragma omp parallel for private(i, dist_center, diff, center, dir_center, dir_force)
	for (i = 0; i < ilen; ++i) {
		PBDBendingConstraint& c = m_ConstraintsBending[i];

		//Using the paper suggested by DevO
		//http://image.diku.dk/kenny/download/kelager.niebe.ea10.pdf
		center = (m_PhyxelPosNew[c.p1] + m_PhyxelPosNew[c.p2] + m_PhyxelPosNew[c.p3]) * 0.3333333f;
		dir_center = m_PhyxelPosNew[c.p3] - center;
		dist_center = dir_center.Length();

		diff = 1.0f - (c.rest_length / dist_center);
		dir_force = dir_center * diff;

		float invmassA = m_PhyxelIsStatic[c.p1] ? 0.0f : m_PhyxelInvMass[c.p1];
		float invmassB = m_PhyxelIsStatic[c.p2] ? 0.0f : m_PhyxelInvMass[c.p2];
		float invmassC = m_PhyxelIsStatic[c.p3] ? 0.0f : m_PhyxelInvMass[c.p3];

		float w = invmassA + invmassB + invmassC;
		if (w > 0.0f)
		{
			w = 1.0f / w;
			m_LargeArrayParallisation[i * 3] = dir_force * (c.k_prime * ((2.0f* invmassA) * w));
			m_LargeArrayParallisation[i * 3 + 1] = dir_force * (c.k_prime * ((2.0f* invmassB) * w));
			m_LargeArrayParallisation[i * 3 + 2] = dir_force  * (-c.k_prime * ((4.0f* invmassC) * w));
		}
		else
		{
			m_LargeArrayParallisation[i * 3] = Vector3(0.0f, 0.0f, 0.0f);
			m_LargeArrayParallisation[i * 3 + 1] = Vector3(0.0f, 0.0f, 0.0f);
			m_LargeArrayParallisation[i * 3 + 2] = Vector3(0.0f, 0.0f, 0.0f);
		}


		//m_PhyxelPosNew[c.p1] += dir_force * (c.k_prime * ((2.0f* m_PhyxelInvMass[c.p1]) * c.w));
		//m_PhyxelPosNew[c.p2] += dir_force * (c.k_prime * ((2.0f* m_PhyxelInvMass[c.p2]) * c.w));
		//m_PhyxelPosNew[c.p3] += dir_force  * (-c.k_prime * ((4.0f* m_PhyxelInvMass[c.p3]) * c.w));
	}

	ilen = (int)m_NumPhyxels;
	Vector3 sum;
#pragma omp parallel for private(i, sum)
	for (i = 0; i < ilen; ++i)
	{
		auto& arr = m_LookupsBending[i];
		sum = Vector3(0, 0, 0);
		for (uint idx : arr)
		{
			sum += m_LargeArrayParallisation[idx];
		}
		sum *= weighting;
		m_PhyxelPosNew[i] += sum;
	}
}
