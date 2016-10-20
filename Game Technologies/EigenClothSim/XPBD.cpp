#include "XPBD.h"
#include <list>

XPBD::XPBD() : ClothBase()
{
	m_SolverIterations = 560;
}

XPBD::~XPBD()
{

}

bool HACKY_IsContainedInDistanceList(uint idxA, uint idxB, const std::list<XPBDDistanceConstraint>& list)
{
	for (auto itr = list.begin(), end = list.end(); itr != end; itr++)
	{
		if (itr->p1 == idxA && itr->p2 == idxB)
			return true;
	}
	return false;
}

void XPBD::InitializeClothDesign(const ClothDesignEntity<3>* design)
{
	ClothBase::InitializeClothDesign(design);
	InitializePBDConstraints(design->GetWidth(), design->GetHeight());
}

void XPBD::InitializePBDConstraints(int num_width, int num_height)
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
	std::list<XPBDDistanceConstraint> d_constraints;
	std::list<XPBDBendingConstraint>  b_constraints;

	uint numBendingConstraints = m_NumQuads * 2;
	uint numDistanceConstraints = 0;

	XPBDDistanceConstraint dc;
	XPBDBendingConstraint bc;

	// Horizontal
	int u = num_width;
	int v = num_height;
	int l1, l2;
	for (l1 = 0; l1 < v; l1++)	// v
		for (l2 = 0; l2 < (u - 1); l2++) {
			if (BuildDistanceConstraint((l1 * u) + l2, (l1 * u) + l2 + 1, m_MaterialShear, dc))
				d_constraints.push_back(dc);
		}

	// Vertical
	for (l1 = 0; l1 < (u); l1++)
		for (l2 = 0; l2 < (v - 1); l2++) {
			if (BuildDistanceConstraint((l2 * u) + l1, ((l2 + 1) * u) + l1, m_MaterialShear, dc))
				d_constraints.push_back(dc);
		}


	// Shearing distance constraint
	for (l1 = 0; l1 < (v - 1); l1++)
		for (l2 = 0; l2 < (u - 1); l2++) {
			if (BuildDistanceConstraint((l1 * u) + l2, ((l1 + 1) * u) + l2 + 1, m_MaterialShear, dc))
				d_constraints.push_back(dc);
			if (BuildDistanceConstraint(((l1 + 1) * u) + l2, (l1 * u) + l2 + 1, m_MaterialShear, dc))
				d_constraints.push_back(dc);
		}



#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
	auto get_idx = [&](int i, int j) {
		return i * num_width + j;
	};

	for (int i = 0; i < num_width; i++) {
		for (int j = 0; j<num_height - 2; j++) {
			if (BuildBendingConstraint(get_idx(i, j), get_idx(i, (j + 1)), get_idx(i, j + 2), m_MaterialBend, bc))
				b_constraints.push_back(bc);
		}
	}
	//add horizontal constraints
	for (int i = 0; i< num_width - 2; i++) {
		for (int j = 0; j < num_height; j++) {
			if (BuildBendingConstraint(get_idx(i, j), get_idx(i + 1, j), get_idx(i + 2, j), m_MaterialBend, bc))
				b_constraints.push_back(bc);
		}
}

	//Bending A-C-B


	//Bending A-D-C

#else
	for (int i = 0; i < num_height - 1; i++)
	{
		for (int j = 0; j < num_width - 1; j++)
		{
			int p1 = i * (num_width) + j;
			int p2 = p1 + 1;
			int p3 = p1 + (num_width);
			int p4 = p3 + 1;

			if ((j + i) % 2) {
				if (BuildBendingConstraint(p3, p2, p1, p4, m_MaterialBend, bc))
					b_constraints.push_back(bc);
			}
			else {
				if (BuildBendingConstraint(p4, p1, p3, p2, m_MaterialBend, bc))
					b_constraints.push_back(bc);
			}

		}
	}
#endif




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

void XPBD::BuildLookupArrays()
{
	/*m_LookupsDistance.clear();
	m_LookupsBending.clear();

	m_LookupsDistance.resize(m_NumPhyxels);
	m_LookupsBending.resize(m_NumPhyxels);

	int i, ilen;
	for (i = 0, ilen = (int)m_ConstraintsDistance.size(); i < ilen; ++i) {
		XPBDDistanceConstraint& c = m_ConstraintsDistance[i];
		m_LookupsDistance[c.p1].push_back(i * 2);
		m_LookupsDistance[c.p2].push_back(i * 2 + 1);
	}

	for (i = 0, ilen = (int)m_ConstraintsBending.size(); i < ilen; ++i) {
		XPBDBendingConstraint& c = m_ConstraintsBending[i];
		m_LookupsBending[c.p1].push_back(i * 3);
		m_LookupsBending[c.p2].push_back(i * 3 + 1);
		m_LookupsBending[c.p3].push_back(i * 3 + 2);
	}*/
}

void XPBD::ComputeVelocity(float subtimestep, const Vector3* in_clothpos, const Vector3* in_clothvel, Vector3* out_clothvel)
{
	float dampFactor = 1.0f - m_MaterialDamp;

#pragma omp parallel for
	for (int i = 0; i < m_NumPhyxels; ++i)
	{
		out_clothvel[i] = in_clothvel[i] * dampFactor + (m_PhyxelForces[i] + m_PhyxelExtForces[i]) * subtimestep;
		out_clothvel[i] = in_clothvel[i] * dampFactor + (m_PhyxelForces[i] + m_PhyxelExtForces[i]) * subtimestep * m_PhyxelInvMass[i];

		m_PhyxelPosNew[i] = in_clothpos[i] + out_clothvel[i] * subtimestep;
	}	
	
	int i, ilen = m_ConstraintsDistance.size();
	for (i = 0; i < ilen; ++i) {
		XPBDDistanceConstraint& c = m_ConstraintsDistance[i];
		c.lamdaij = 0.0f;
	}
	ilen = m_ConstraintsBending.size();
	for (i = 0; i < ilen; ++i) {
		XPBDBendingConstraint& c = m_ConstraintsBending[i];
		c.lambdaij = 0.0f;
	}
	float weighting = 1.0f;
	for (uint si = 0; si < m_SolverIterations; ++si) {
		//weighting = 1.0f - (float)(si * si) / (float)(m_SolverIterations * m_SolverIterations);
		SolveDistanceConstraints(weighting);
		SolveBendingConstraints(weighting);
		SolveSphereConstraints(weighting);
	}

	float inv_timestep = 1.0f / subtimestep;
#pragma omp parallel for
	for (int i = 0; i < m_NumPhyxels; ++i)
	{
		out_clothvel[i] = (m_PhyxelPosNew[i] - in_clothpos[i]) * inv_timestep;
	}
}

bool XPBD::BuildDistanceConstraint(uint idxA, uint idxB, float k, XPBDDistanceConstraint& out_constraint)
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

#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
bool XPBD::BuildBendingConstraint(uint idxA, uint idxB, uint idxC, float k, XPBDBendingConstraint& out_constraint)
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
#else
bool XPBD::BuildBendingConstraint(uint idxA, uint idxB, uint idxC, uint idxD, float k, XPBDBendingConstraint& out_constraint)
{
	out_constraint.p1 = idxA;
	out_constraint.p2 = idxB;
	out_constraint.p3 = idxC;
	out_constraint.p4 = idxD;
	out_constraint.w1 = m_PhyxelInvMass[idxA] + m_PhyxelInvMass[idxB] + 2 * m_PhyxelInvMass[idxC];
	out_constraint.w2 = m_PhyxelInvMass[idxA] + m_PhyxelInvMass[idxB] + 2 * m_PhyxelInvMass[idxD];
	
	Vector3 center1 = (m_ClothStatesPos[idxA] + m_ClothStatesPos[idxB] + m_ClothStatesPos[idxC]) / 3.0f;
	Vector3 center2 = (m_ClothStatesPos[idxA] + m_ClothStatesPos[idxB] + m_ClothStatesPos[idxD]) / 3.0f;

	out_constraint.rest_length1 = (m_ClothStatesPos[idxC] - center1).Length();
	out_constraint.rest_length2 = (m_ClothStatesPos[idxD] - center2).Length();
	out_constraint.k = k;

	out_constraint.k_prime = 1.0f - pow((1.0f - out_constraint.k), 1.0f / m_SolverIterations);  //1.0f-pow((1.0f-c.k), 1.0f/ns);
	if (out_constraint.k_prime>1.0)
		out_constraint.k_prime = 1.0;



	Vector3 n1 = Vector3::Cross(m_ClothStatesPos[idxA] - m_ClothStatesPos[idxB], m_ClothStatesPos[idxC] - m_ClothStatesPos[idxB]); n1.Normalise();
	Vector3 n2 = Vector3::Cross(m_ClothStatesPos[idxA] - m_ClothStatesPos[idxB], m_ClothStatesPos[idxD] - m_ClothStatesPos[idxB]); n2.Normalise();;
	out_constraint.phi0 = acos(Vector3::Dot(n1, n2));

	return true;
}
#endif 

void XPBD::SolveDistanceConstraints(float weighting)
{
	Vector3 dir, dP;
	float len;

	int i, ilen = m_ConstraintsDistance.size();
//#pragma omp parallel for private(i, len, dir, dP)
	for (i = 0; i < ilen; ++i) {
		XPBDDistanceConstraint& c = m_ConstraintsDistance[i];
		dir = m_PhyxelPosNew[c.p1] - m_PhyxelPosNew[c.p2];

		float invmassA = m_PhyxelIsStatic[c.p1] ? 0.0f : m_PhyxelInvMass[c.p1];
		float invmassB = m_PhyxelIsStatic[c.p2] ? 0.0f : m_PhyxelInvMass[c.p2];
		float w = invmassA + invmassB;

		if (w > 0.0f)
		{
			len = dir.Length();
			//if (len <= 0.00001f)
			//	return;

			//dP = dir / len * (c.k_prime * (len - c.rest_length));
			//dP = dP / w;

			//m_LargeArrayParallisation[i * 2] = -dP * invmassA;
			//m_LargeArrayParallisation[i * 2 + 1] = dP * invmassB;
			//m_PhyxelPosNew[c.p1] -= dP * m_PhyxelInvMass[c.p1];
			//m_PhyxelPosNew[c.p2] += dP * m_PhyxelInvMass[c.p2];

			float lambda = (len - c.rest_length - c.k * c.lamdaij) / (w + c.k);
			dP = dir / len * lambda;


			c.lamdaij += lambda;

			
			m_PhyxelPosNew[c.p1] -= dP * invmassA;
			m_PhyxelPosNew[c.p2] += dP * invmassB;
		}
		else
		{
		//	m_LargeArrayParallisation[i * 2] = Vector3(0.0f, 0.0f, 0.0f);
		//	m_LargeArrayParallisation[i * 2 + 1] = Vector3(0.0f, 0.0f, 0.0f);
		}
	}

	ilen = (int)m_NumPhyxels;

/*	Vector3 sum;
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
	}*/
}


#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
void XPBD::SolveBendingConstraints(float weighting)
{
	Vector3 center, dir_center, dir_force;
	float dist_center, diff;
	int i, ilen;
	ilen = (int)m_ConstraintsBending.size();

//#pragma omp parallel for private(i, dist_center, diff, center, dir_center, dir_force)
	for (i = 0; i < ilen; ++i) {
		XPBDBendingConstraint& c = m_ConstraintsBending[i];

		//Using the paper suggested by DevO
		//http://image.diku.dk/kenny/download/kelager.niebe.ea10.pdf
		center = (m_PhyxelPosNew[c.p1] + m_PhyxelPosNew[c.p2] + m_PhyxelPosNew[c.p3]) * 0.3333333f;
		dir_center = m_PhyxelPosNew[c.p3] - center;
		dist_center = dir_center.Length();

		float invmassA = m_PhyxelIsStatic[c.p1] ? 0.0f : m_PhyxelInvMass[c.p1];
		float invmassB = m_PhyxelIsStatic[c.p2] ? 0.0f : m_PhyxelInvMass[c.p2];
		float invmassC = m_PhyxelIsStatic[c.p3] ? 0.0f : m_PhyxelInvMass[c.p3];

		float w = invmassA + invmassB + invmassC;

		diff = (1.0f - (c.rest_length / dist_center) - c.k * c.lambdaij) / (w + c.k);
		c.lambdaij += diff;
		dir_force = dir_center * diff;


		/*if (w > 0.0f)
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
		}*/


			m_PhyxelPosNew[c.p1] += dir_force * (c.k * 2.0f* invmassA);
			m_PhyxelPosNew[c.p2] += dir_force * (c.k * 2.0f* invmassB);
			m_PhyxelPosNew[c.p3] += dir_force  * (-c.k * 4.0f* invmassC);

	}

	/*	ilen = (int)m_NumPhyxels;
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
	}*/
}
#else
void XPBD::SolveBendingConstraints(float weighting)
{
	Vector3 center, dir_center, dir_force;
	float dist_center, diff;
	int i, ilen;
	ilen = (int)m_ConstraintsBending.size();

	//#pragma omp parallel for private(i, dist_center, diff, center, dir_center, dir_force)
	for (i = 0; i < ilen; ++i) {
		XPBDBendingConstraint& c = m_ConstraintsBending[i];

		//Using the dihedral angle approach of the position based dynamics		
		float d = 0, phi = 0, i_d = 0;
		Vector3 n1 = Vector3(0.0f, 0.0f, 0.0f), n2 = Vector3(0.0f, 0.0f, 0.0f);

		Vector3 p1 = m_PhyxelPosNew[c.p1];
		Vector3 p2 = m_PhyxelPosNew[c.p2] - p1;
		Vector3 p3 = m_PhyxelPosNew[c.p3] - p1;
		Vector3 p4 = m_PhyxelPosNew[c.p4] - p1;

		Vector3 p2p3 = Vector3::Cross(p2, p3);
		Vector3 p2p4 = Vector3::Cross(p2, p4);

		float lenp2p3 = (p2p3).Length();

		if (lenp2p3 == 0.0) { return; } //need to handle this case.

		float lenp2p4 = (p2p4).Length();

		if (lenp2p4 == 0.0) { return; } //need to handle this case.

		n1 = (p2p3); n1.Normalise();
		n2 = (p2p4); n2.Normalise();

		d = Vector3::Dot(n1, n2);
		phi = acos(d);

		//try to catch invalid values that will return NaN.
		// sqrt(1 - (1.0001*1.0001)) = NaN 
		// sqrt(1 - (-1.0001*-1.0001)) = NaN 
		if (d<-1.0)
			d = -1.0;
		else if (d>1.0)
			d = 1.0; //d = clamp(d,-1.0,1.0);

					 //in both case sqrt(1-d*d) will be zero and nothing will be done.
					 //0° case, the triangles are facing in the opposite direction, folded together.
		if (d == -1.0) {
			phi = PI;  //acos(-1.0) == PI
			if (phi == c.phi0)
				return; //nothing to do 

						//in this case one just need to push 
						//vertices 1 and 2 in n1 and n2 directions, 
						//so the constrain will do the work in second iterations.
			if (c.p1 != 0)
				m_PhyxelPosNew[c.p3] += n1 / 100.0f;

			if (c.p2 != 0)
				m_PhyxelPosNew[c.p4] += n2 / 100.0f;

			return;
		}
		if (d == 1.0) { //180° case, the triangles are planar
			phi = 0.0;  //acos(1.0) == 0.0
			if (phi == c.phi0)
				return; //nothing to do 
		}

		i_d = sqrt(1 - (d*d))*(phi - c.phi0);

		Vector3 p2n1 = Vector3::Cross(p2, n1);
		Vector3 p2n2 = Vector3::Cross(p2, n2);
		Vector3 p3n2 = Vector3::Cross(p3, n2);
		Vector3 p4n1 = Vector3::Cross(p4, n1);
		Vector3 n1p2 = -p2n1;
		Vector3 n2p2 = -p2n2;
		Vector3 n1p3 = Vector3::Cross(n1, p3);
		Vector3 n2p4 = Vector3::Cross(n2, p4);

		Vector3 q3 = (p2n2 + n1p2*d) / lenp2p3;
		Vector3 q4 = (p2n1 + n2p2*d) / lenp2p4;
		Vector3 q2 = (-(p3n2 + n1p3*d) / lenp2p3) - ((p4n1 + n2p4*d) / lenp2p4);

		Vector3 q1 = -q2 - q3 - q4;

		float q1_len2 = Vector3::Dot(q1, q1);// glm::length(q1)*glm::length(q1);
		float q2_len2 = Vector3::Dot(q2, q2);// glm::length(q2)*glm::length(q1);
		float q3_len2 = Vector3::Dot(q3, q3);// glm::length(q3)*glm::length(q1);
		float q4_len2 = Vector3::Dot(q4, q4);// glm::length(q4)*glm::length(q1); 

		float sum = m_PhyxelInvMass[c.p1] * (q1_len2)+
			m_PhyxelInvMass[c.p2] * (q2_len2)+
			m_PhyxelInvMass[c.p3] * (q3_len2)+
			m_PhyxelInvMass[c.p4] * (q4_len2);

		Vector3 dP1 = q1 * -((m_PhyxelInvMass[c.p1] * i_d) / sum);
		Vector3 dP2 = q2 * -((m_PhyxelInvMass[c.p2] * i_d) / sum);
		Vector3 dP3 = q3 * -((m_PhyxelInvMass[c.p3] * i_d) / sum);
		Vector3 dP4 = q4  * -((m_PhyxelInvMass[c.p4] * i_d) / sum);

		if (m_PhyxelInvMass[c.p1] > 0.0) {
			m_PhyxelPosNew[c.p1] += dP1*c.k;
		}
		if (m_PhyxelInvMass[c.p2] > 0.0) {
			m_PhyxelPosNew[c.p2] += dP2*c.k;
		}
		if (m_PhyxelInvMass[c.p3] > 0.0) {
			m_PhyxelPosNew[c.p3] += dP3*c.k;
		}
		if (m_PhyxelInvMass[c.p4] > 0.0) {
			m_PhyxelPosNew[c.p4] += dP4*c.k;
		}
	}
}
#endif

void XPBD::SolveSphereConstraints(float weighting)
{
	for (uint j = 0; j < m_SphereConstraints.size(); ++j)
	{
		XPBDSphereConstraint& sphere = m_SphereConstraints[j];

#pragma omp parallel for
		for (int i = 0; i < m_NumPhyxels; ++i)
		{
			Vector3 pos = m_PhyxelPosNew[i];

			Vector3 axis = pos - sphere.centre;
			float radiusSq = sphere.radius * sphere.radius;
			float distSquared = Vector3::Dot(axis, axis);

			if (distSquared < radiusSq)
			{
				float dist = sqrt(distSquared);

				float excess = 1.0f - dist / sphere.radius;
				pos += axis * excess;

				m_PhyxelPosNew[i] = pos;
			}
		}
	}
}