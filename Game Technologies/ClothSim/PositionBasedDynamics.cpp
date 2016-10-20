#include "PositionBasedDynamics.h"
#include  <nclgl\Vector2.h>

const float global_dampening = 0.999f;
const float kBend = 0.9f;
const float kStretch = 0.75f;
const float kShear = 0.5f;
const float kDamp = 0.001f;
const float mass_density = 0.26f;

PositionBasedDynamicsMS::PositionBasedDynamicsMS()
{
}

PositionBasedDynamicsMS::~PositionBasedDynamicsMS()
{
}


#pragma region SimulationFunctionality

void PositionBasedDynamicsMS::simulation_OnClothDesignChanged(ClothDesignEntity* design)
{
	//Gen Vertices
	const std::vector<Vertex>& vertices = design->Vertices();
	m_NumTotal = vertices.size();

	m_PhyxelTexCoords.resize(m_NumTotal);
	m_PhyxelsPos.resize(m_NumTotal);
	m_PhyxelsPosNew.resize(m_NumTotal);
	m_PhyxelsVel.resize(m_NumTotal);
	m_PhyxelForces.resize(m_NumTotal);
	m_PhyxelIsStatic.resize(m_NumTotal);
	m_PhyxelsMass.resize(m_NumTotal);
	m_PhyxelsInvMass.resize(m_NumTotal);
	m_PhyxelNormals.resize(m_NumTotal);
	m_PhyxelRi.resize(m_NumTotal);

	memset(&m_PhyxelsVel[0], 0, m_NumTotal * sizeof(Vector3));
	memset(&m_PhyxelForces[0], 0, m_NumTotal * sizeof(Vector3));

	float uniform_mass = mass_density / m_NumTotal;
	float inv_uniform_mass = 1.0f / uniform_mass;
	for (unsigned int i = 0; i < m_NumTotal; ++i)
	{
		m_PhyxelsMass[i] = uniform_mass;
		m_PhyxelsInvMass[i] = inv_uniform_mass;
		//m_PhyxelForces[i] = GRAVITY;

		const Vertex& v = vertices[i];
		m_PhyxelTexCoords[i] = v.texCoord;
		m_PhyxelsPos[i] = v.position;
		m_PhyxelsPosNew[i] = v.position;
		m_PhyxelIsStatic[i] = v.flags & VERTEXFLAGS_IS_STATIC;
		m_PhyxelForces[i] = v.force;

		if (m_PhyxelIsStatic[i])
			m_PhyxelsInvMass[i] = 0.0f;
	}


	//Copy Triangles (For generic use not with simulation)
	m_NumTriangles = design->Triangles().size();
	m_Triangles.resize(m_NumTriangles);
	memcpy(&m_Triangles[0], &design->Triangles()[0], sizeof(Triangle) * m_NumTriangles);

	//Setup Constraints
	unsigned int l1, l2, divisor = sqrt(m_NumTotal);;
	// Horizontal
	for (l1 = 0; l1 < divisor; ++l1)	// v
		for (l2 = 0; l2 < (divisor - 1); ++l2) {
		AddDistanceConstraint((l1 * divisor) + l2, (l1 * divisor) + l2 + 1, kStretch);
		}

	// Vertical
	for (l1 = 0; l1 < (divisor); ++l1)
		for (l2 = 0; l2 < (divisor - 1); ++l2) {
		AddDistanceConstraint((l2 * divisor) + l1, ((l2 + 1) * divisor) + l1, kStretch);
		}


	// Shearing distance constraint
	for (l1 = 0; l1 < (divisor - 1); ++l1)
		for (l2 = 0; l2 < (divisor - 1); ++l2) {
		AddDistanceConstraint((l1 * divisor) + l2, ((l1 + 1) * divisor) + l2 + 1, kShear);
		AddDistanceConstraint(((l1 + 1) * divisor) + l2, (l1 * divisor) + l2 + 1, kShear);
		}


	// create bending constraints	
	//add vertical constraints
	for (l1 = 0; l1 < divisor; ++l1) {
		for (l2 = 0; l2 < divisor - 2; ++l2) {
			unsigned int a = (l1 * divisor) + l2;
			unsigned int b = (l1 * divisor) + l2 + 1;
			unsigned int c = (l1 * divisor) + l2 + 2;
			AddBendingConstraint(a, b, c, kBend);
		}
	}
	//add horizontal constraints
	for (l1 = 0; l1 < divisor - 2; ++l1) {
		for (l2 = 0; l2 < divisor; ++l2) {
			unsigned int a = (l1 * divisor) + l2;
			unsigned int b = ((l1 + 1) * divisor) + l2;
			unsigned int c = ((l1 + 2) * divisor) + l2;
			AddBendingConstraint(a, b, c, kBend);
		}
	}
}

void PositionBasedDynamicsMS::AddDistanceConstraint(int a, int b, float k)
{
	float invmass = m_PhyxelsInvMass[a] + m_PhyxelsInvMass[b];;
	if (invmass < 0.0001f)
		return;

	PBDDistanceConstraint c;
	c.p1 = a;
	c.p2 = b;
	c.k = k;
	c.k_prime = (1.0f - pow((1.0f - c.k), 1.0f / PBDSOLVER_ITERATIONS));  //1.0f-pow((1.0f-c.k), 1.0f/ns);

	if (c.k_prime>1.0)
		c.k_prime = 1.0;

	//c.k_prime = c.k_prime;

	Vector3 deltaP = m_PhyxelsPos[c.p1] - m_PhyxelsPos[c.p2];
	c.rest_length = deltaP.Length();

	d_constraints.push_back(c);
}

void PositionBasedDynamicsMS::AddBendingConstraint(int pa, int pb, int pc, float k)
{
	float invmass = (m_PhyxelsInvMass[pa] + m_PhyxelsInvMass[pb] + 2 * m_PhyxelsInvMass[pc]);
	if (invmass < 0.0000001f)
		return;

	PBDBendingConstraint c;
	c.p1 = pa;
	c.p2 = pb;
	c.p3 = pc;

	c.w = 1.0f / invmass;
	Vector3 center = (m_PhyxelsPos[pa] + m_PhyxelsPos[pb] + m_PhyxelsPos[pc]) * 0.3333333f;
	c.rest_length = (m_PhyxelsPos[pc] - center).Length();
	c.k = k;
	c.k_prime = 1.0f - pow((1.0f - c.k), 1.0f / PBDSOLVER_ITERATIONS);  //1.0f-pow((1.0f-c.k), 1.0f/ns);
	if (c.k_prime>1.0)
		c.k_prime = 1.0;
	b_constraints.push_back(c);
}

void PositionBasedDynamicsMS::Simulation_StepSimulation(float dt)
{
	unsigned int i;

	//Explicitly integrate positions
	Vector3 Xcm = Vector3(0, 0, 0);
	Vector3 Vcm = Vector3(0, 0, 0);
	float sumM = 0.0f;
	for (i = 0; i < m_NumTotal; ++i) {

		m_PhyxelsVel[i] = m_PhyxelsVel[i] * global_dampening; //global velocity dampening !!!		
		m_PhyxelsVel[i] = m_PhyxelsVel[i] + (m_PhyxelForces[i] * dt) / m_PhyxelsMass[i] + GRAVITY * dt;

		//calculate the center of mass's position 
		//and velocity for damping calc
		Xcm += (m_PhyxelsPos[i] * m_PhyxelsMass[i]);
		Vcm += (m_PhyxelsVel[i] * m_PhyxelsMass[i]);
		sumM += m_PhyxelsMass[i];
	}
	Xcm = Xcm / sumM;
	Vcm = Vcm / sumM;

	Matrix3 I;
	Vector3 L = Vector3(0, 0, 0);
	Vector3 w = Vector3(0, 0, 0);//angular velocity


	for (i = 0; i < m_NumTotal; ++i) {
		m_PhyxelRi[i] = (m_PhyxelsPos[i] - Xcm);

		L += Vector3::Cross(m_PhyxelRi[i], m_PhyxelsVel[i] * m_PhyxelsMass[i]);

		//thanks to DevO for pointing this and these notes really helped.
		//http://www.sccg.sk/~onderik/phd/ca2010/ca10_lesson11.pdf

		Matrix3 tmp = Matrix3(0, -m_PhyxelRi[i].z, m_PhyxelRi[i].y,
			m_PhyxelRi[i].z, 0, -m_PhyxelRi[i].x,
			-m_PhyxelRi[i].y, m_PhyxelRi[i].x, 0);
		I += (tmp*Matrix3::Transpose(tmp))*m_PhyxelsMass[i];
	}

	w = Matrix3::Inverse(I)*L;

	//apply center of mass damping
	for (i = 0; i < m_NumTotal; ++i) {
		Vector3 delVi = Vcm + Vector3::Cross(w, m_PhyxelRi[i]) - m_PhyxelsVel[i];
		m_PhyxelsVel[i] += delVi * kDamp;
	}

	//calculate predicted position
	for (i = 0; i < m_NumTotal; ++i) {
		if (m_PhyxelIsStatic[i])
		{
			m_PhyxelsPosNew[i] = m_PhyxelsPos[i]; //fixed points
		}
		else
		{
			m_PhyxelsPosNew[i] = m_PhyxelsPos[i] + (m_PhyxelsVel[i] * dt);
		}
	}

	//Compute Rotation Matrix for each triangle
	UpdateInternalConstraints(dt);

	//Update external constraints
	for (uint i = 0; i < m_colEllipsoids.size(); ++i)
	{
		EllipsoidCollision(m_colEllipsoids[i]);
	}

	BuildNormals();

	//Update Phyxel Positions
	float inv_dt = 1.0f / dt;
	for (i = 0; i < m_NumTotal; i++) {
		m_PhyxelsVel[i] = (m_PhyxelsPosNew[i] - m_PhyxelsPos[i])*inv_dt;
		m_PhyxelsPos[i] = m_PhyxelsPosNew[i];
	}

	//memset(&m_PhyxelForces[0], 0, m_NumTotal * sizeof(Vector3));
}

void PositionBasedDynamicsMS::BuildNormals()
{
	unsigned int i;
	memset(&m_PhyxelNormals[0], 0, m_NumTotal * sizeof(Vector3));

	Vector3 normal;
	//Build Rotations and Sum Up Normals
	for (i = 0; i < m_NumTriangles; ++i)
	{
		Triangle& t = m_Triangles[i];

		normal = Vector3::Cross(m_PhyxelsPos[t.vert_a] - m_PhyxelsPos[t.vert_c],
			m_PhyxelsPos[t.vert_b] - m_PhyxelsPos[t.vert_c]);

		normal = normal.Normalise();
		m_PhyxelNormals[t.vert_a] += normal;
		m_PhyxelNormals[t.vert_b] += normal;
		m_PhyxelNormals[t.vert_c] += normal;
	}

	//Normalise Vertex Normals
	for (i = 0; i < m_NumTotal; ++i)
	{
		m_PhyxelNormals[i].Normalise();
	}
}

void PositionBasedDynamicsMS::UpdateInternalConstraints(float dt)
{
	unsigned int si;

	for (si = 0; si < PBDSOLVER_ITERATIONS; ++si) {
		UpdateDistanceConstraints();
		UpdateBendingConstraints();
		GroundCollision();
	}
}
void PositionBasedDynamicsMS::UpdateDistanceConstraints()
{
	Vector3 dir, dP;
	float len;
	uint i, ilen;
	for (i = 0, ilen = d_constraints.size(); i < ilen; ++i) {
		PBDDistanceConstraint& c = d_constraints[i];
		dir = m_PhyxelsPosNew[c.p1] - m_PhyxelsPosNew[c.p2];

		len = dir.Length();
		//if (len <= 0.00001f)
		//	return;

		dP = dir / len * (c.k_prime * (len - c.rest_length));
		dP = dP / (m_PhyxelsInvMass[c.p1] + m_PhyxelsInvMass[c.p2]);

		m_PhyxelsPosNew[c.p1] -= dP * m_PhyxelsInvMass[c.p1];
		m_PhyxelsPosNew[c.p2] += dP * m_PhyxelsInvMass[c.p2];
	}
}

void PositionBasedDynamicsMS::UpdateBendingConstraints()
{
	Vector3 center, dir_center, dir_force;
	float dist_center, diff;
	uint i, ilen;
	for (i = 0, ilen = b_constraints.size(); i < ilen; ++i) {
		PBDBendingConstraint& c = b_constraints[i];

		//Using the paper suggested by DevO
		//http://image.diku.dk/kenny/download/kelager.niebe.ea10.pdf
		center = (m_PhyxelsPosNew[c.p1] + m_PhyxelsPosNew[c.p2] + m_PhyxelsPosNew[c.p3]) * 0.3333333f;
		dir_center = m_PhyxelsPosNew[c.p3] - center;
		dist_center = dir_center.Length();

		diff = 1.0f - (c.rest_length / dist_center);
		dir_force = dir_center * diff;

		m_PhyxelsPosNew[c.p1] += dir_force * (c.k_prime * ((2.0f*m_PhyxelsInvMass[c.p1]) * c.w));
		m_PhyxelsPosNew[c.p2] += dir_force * (c.k_prime * ((2.0f*m_PhyxelsInvMass[c.p2]) * c.w));
		m_PhyxelsPosNew[c.p3] += dir_force  * (-c.k_prime * ((4.0f*m_PhyxelsInvMass[c.p3]) * c.w));
	}
}

#pragma endregion SimulationFunctionality


#pragma region RenderFunctionality
void PositionBasedDynamicsMS::Render_DrawingToVisualiser()
{
	uint i, ilen;
	for (i = 0, ilen = d_constraints.size(); i < ilen; ++i) {
		PBDDistanceConstraint& c = d_constraints[i];

		NCLDebug::DrawThickLine(m_PhyxelsPosNew[c.p1], m_PhyxelsPosNew[c.p2], 0.001f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	}

}
#pragma endregion RenderFunctionality

void PositionBasedDynamicsMS::GroundCollision()
{
	for (uint i = 0; i< m_NumTotal; i++) {
		if (m_PhyxelsPosNew[i].y < 0) //collision with ground
		{
			m_PhyxelsPosNew[i].y = 0;
			m_PhyxelsPos[i] = m_PhyxelsPosNew[i]; //Completely in-ellastic collision
			m_PhyxelsVel[i] = Vector3(0, 0, 0);
		}
		//Friction needed??
		// - Very hacky plane collision anyway, should really use proper generic collision response
	}
}

void PositionBasedDynamicsMS::EllipsoidCollision(const Ellipsoid& e)
{
	for (uint i = 0; i < m_NumTotal; i++) {
		if (m_PhyxelIsStatic[i])
			continue;

		Vector4 X_0 = (e.invTransform * Vector4(m_PhyxelsPosNew[i], 1));
		Vector3 delta0 = Vector3(X_0.x, X_0.y, X_0.z);
		float distance = delta0.Length();

		if (distance < e.radius) {
			delta0 = delta0 * ((e.radius - distance) / distance);

			// Transform the delta back to original space
			Vector3 delta;
			Vector3 transformInv;

			transformInv = Vector3(e.Transform[0], e.Transform[1], e.Transform[2]);
			transformInv = transformInv / Vector3::Dot(transformInv, transformInv);
			delta.x = Vector3::Dot(delta0, transformInv);
			
			transformInv = Vector3(e.Transform[4], e.Transform[5], e.Transform[6]);
			transformInv = transformInv / Vector3::Dot(transformInv, transformInv);
			delta.y = Vector3::Dot(delta0, transformInv);
			
			transformInv = Vector3(e.Transform[8], e.Transform[9], e.Transform[10]);
			transformInv = transformInv / Vector3::Dot(transformInv, transformInv);
			delta.z = Vector3::Dot(delta0, transformInv);


			m_PhyxelsPosNew[i] += delta;

			Vector3 reboundVec = m_PhyxelsPosNew[i] - m_PhyxelsPos[i];
			m_PhyxelsPos[i] = m_PhyxelsPosNew[i] - reboundVec * 0.9;

			m_PhyxelsPos[i] = m_PhyxelsPosNew[i];// = Vector3(0, 0, 0); //Completely in-ellastic collision
			m_PhyxelsVel[i] = (m_PhyxelsPosNew[i] - m_PhyxelsPos[i]);
		}
	}
}