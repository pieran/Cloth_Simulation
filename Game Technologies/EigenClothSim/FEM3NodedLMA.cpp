#include "FEM3NodedLMA.h"
#include <ncltech\NCLDebug.h>




FEM3NodedLMA::FEM3NodedLMA()
{
	m_TimeStep = max_timestep;
	m_StepCounter = 0;
	m_StepsBeforeIncrease = 10;
}

FEM3NodedLMA::~FEM3NodedLMA()
{
}

void FEM3NodedLMA::simulation_OnClothDesignChanged(ClothDesignEntity<3>* design)
{

	//Gen Vertices
	const std::vector<Vertex>& vertices = design->Vertices();
	m_NumTotal = vertices.size();


	m_PhyxelsPos.resize(m_NumTotal);
	m_PhyxelsPosTemp.resize(m_NumTotal);
	m_PhyxelsPosInitial.resize(m_NumTotal);
	m_PhyxelsVel.resize(m_NumTotal);
	m_PhyxelsVelChange.resize(m_NumTotal);
	m_PhyxelForces.resize(m_NumTotal);
	m_PhyxelExtForces.resize(m_NumTotal);

	m_PhyxelNormals.resize(m_NumTotal);
	m_PhyxelIsStatic.resize(m_NumTotal);
	m_PhyxelsMass.resize(m_NumTotal);
	m_PhyxelTexCoords.resize(m_NumTotal);

	memset(&m_PhyxelsVel[0].x, 0, m_NumTotal * sizeof(Vector3));
	memset(&m_PhyxelsVelChange[0].x, 0, m_NumTotal * sizeof(Vector3));
	memset(&m_PhyxelExtForces[0].x, 0, m_NumTotal * sizeof(Vector3));
	memset(&m_PhyxelNormals[0].x, 0, m_NumTotal * sizeof(Vector3));

	for (unsigned int i = 0; i < m_NumTotal; ++i)
	{
		const Vertex& v = vertices[i];

		m_PhyxelsPos[i] = Vector3(v.position.x(), v.position.y(), v.position.z());
		m_PhyxelsPosInitial[i] = Vector3(v.initial_position.x(), v.initial_position.y(), v.initial_position.z());
		m_PhyxelIsStatic[i] = v.flags & VERTEXFLAGS_IS_STATIC;

		if (!m_PhyxelIsStatic[i])
		{
			m_PhyxelForces[i] = GRAVITY / float(m_NumTotal);
		}


		m_PhyxelTexCoords[i] = Vector2(v.texCoord.x(), v.texCoord.y());
	}




	//Gen Triangles
	const std::vector<Triangle>& triangles = design->Triangles();
	m_NumTriangles = triangles.size();
	m_RenderIndices.resize(m_NumTriangles * 3);
	m_Triangles.resize(m_NumTriangles);
	m_Quads.resize(design->Quads().size());

	memcpy(&m_Quads[0].vert_a, &(design->Quads()[0].vert_a), m_Quads.size() * sizeof(Quad));

	m_PhyxelsTriangleLookup.resize(m_NumTotal);
	m_TrianglesB.resize(m_NumTriangles * 3);
	m_TrianglesAii.resize(m_NumTriangles * 3);
	m_TrianglesAij.resize(m_NumTriangles * 3);

	float totalArea = 0.0f;
	for (unsigned int i = 0; i < m_NumTriangles; ++i)
	{
		unsigned int i3 = i * 3;
		uint idxA = triangles[i].verts[0];
		uint idxB = triangles[i].verts[1];
		uint idxC = triangles[i].verts[2];

		m_Triangles[i].phyxels[0] = idxA;
		m_Triangles[i].phyxels[1] = idxB;
		m_Triangles[i].phyxels[2] = idxC;

		m_RenderIndices[i3] = idxA;
		m_RenderIndices[i3 + 1] = idxB;
		m_RenderIndices[i3 + 2] = idxC;

		m_PhyxelsTriangleLookup[idxA].push_back(i3 + 0);
		m_PhyxelsTriangleLookup[idxB].push_back(i3 + 1);
		m_PhyxelsTriangleLookup[idxC].push_back(i3 + 2);

		GenFEMTriangle(i, m_Triangles[i]);

		totalArea += m_Triangles[i].Area;
	}

	float uniform_mass = (totalArea * mass_density) / m_NumTotal;
	for (unsigned int i = 0; i < m_NumTotal; ++i)
	{
		m_PhyxelsMass[i] = uniform_mass;
	}


	m_Solver.AllocateMemory(m_NumTotal);
	m_Solver.m_A.resize(m_NumTriangles, m_NumTotal);

	//Build Global Matricies (To be deleted at the start of the simulation)
	//SimpleCorotatedBuildAMatrix(0.0f);
	for (int i = 0; i < (int)m_NumTriangles; ++i)
	{
		m_Solver.m_A[i].idxA = m_Triangles[i].phyxels[0];
		m_Solver.m_A[i].idxB = m_Triangles[i].phyxels[1];
		m_Solver.m_A[i].idxC = m_Triangles[i].phyxels[2];
		m_Solver.m_A.SetMass(i, m_PhyxelsMass[i]);
	}



	for (uint i = 0; i < m_NumTotal; ++i)
	{
		m_PhyxelExtForces[i] = Vector3(0, 0, 0);
	}
	m_Solver.ResetMemory();
	UpdateConstraints();
}

void  FEM3NodedLMA::UpdateConstraints()
{
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTotal; ++i)
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

void FEM3NodedLMA::Simulation_StepSimulation(float dt)
{
	m_ProfilingTotalTime.BeginTiming();
	m_ProfilingBuildMatricies.ResetTotalMs();
	m_ProfilingConjugateGradient.ResetTotalMs();
	m_ProfilingExternalCollisions.ResetTotalMs();
	m_ProfilingRotationNormals.ResetTotalMs();
	m_Solver.ResetProfilingData();

	//Handle External Constraints
	{
		m_ProfilingExternalCollisions.BeginTiming();
		UpdateConstraints();
		for (const Ellipsoid& e : m_Ellipsoids)
		{
			CollideEllipsoid(e);
		}

		const Vector3 col_norm = Vector3(0, 1, 0);
#pragma omp parallel for
		for (int i = 0; i < (int)m_NumTotal; ++i)
		{
			if (m_PhyxelsPos[i].y < 0.1f)
			{
				m_PhyxelsVel[i] -= col_norm * Vector3::Dot(m_PhyxelsVel[i], col_norm);
				m_PhyxelsVel[i] = m_PhyxelsVel[i] * 0.9f;
				m_PhyxelsVel[i].y = (0.1f - m_PhyxelsPos[i].y) * 60.0f;

				m_Solver.m_Constraints[i] = Matrix3::Identity - Matrix3::OuterProduct(col_norm, col_norm);
			}
		}
		m_ProfilingExternalCollisions.EndTimingAdditive();
	}


	//Compute Rotation Matrix for each triangle
	dt = min(dt, 1.0f / 60.0f);

	uint iterations = 0;
	for (float timestep_accum = 0.0f; timestep_accum < dt;)
	{
		//iterations++;
		//Handle External Constraints
		//if (iterations == 10)
		{
			m_ProfilingExternalCollisions.BeginTiming();
			iterations = 0;
			UpdateConstraints();
			for (const Ellipsoid& e : m_Ellipsoids)
			{
				CollideEllipsoid(e);
			}

			const Vector3 col_norm = Vector3(0, 1, 0);
#pragma omp parallel for
			for (int i = 0; i < (int)m_NumTotal; ++i)
			{
				if (m_PhyxelsPos[i].y < 0)
				{
					m_PhyxelsVel[i] -= col_norm * Vector3::Dot(m_PhyxelsVel[i], col_norm);
					m_PhyxelsVel[i].y -= m_PhyxelsPos[i].y * 60.0f;

					m_Solver.m_Constraints[i] = Matrix3::Identity - Matrix3::OuterProduct(col_norm, col_norm);
					m_Solver.m_X[i] = m_PhyxelsVel[i];
				}
			}
			m_ProfilingExternalCollisions.EndTimingAdditive();
		}



		for (uint i = 0; i < m_NumTotal; ++i)
		{
			if (m_PhyxelIsStatic[i])
			{
				m_PhyxelsVel[i].z = sin(angle) * 0.01f;// *(static_itr < 33 ? 1.0f : -1.0f);
				m_PhyxelsVel[i].y = sin(angle) * 5.0f;// *(static_itr++ < 33 ? 1.0f : -1.0f);

				m_Solver.m_X[i] = m_PhyxelsVel[i];
			}
		}

		m_ProfilingRotationNormals.BeginTiming();
		BuildRotationAndNormals();
		m_ProfilingRotationNormals.EndTimingAdditive();

		bool valid_timestep = false;
		while (!valid_timestep)
		{
			m_ProfilingConjugateGradient.BeginTiming();
			m_Solver.ResetMemory();
			m_ProfilingConjugateGradient.EndTimingAdditive();

			m_ProfilingBuildMatricies.BeginTiming();
			SimpleCorotatedBuildAMatrix(m_TimeStep);
			m_ProfilingBuildMatricies.EndTimingAdditive();

			m_ProfilingConjugateGradient.BeginTiming();
			//m_Solver.SolveWithPreviousResult();
			m_Solver.SolveWithGuess(m_PhyxelsVel);
			m_ProfilingConjugateGradient.EndTimingAdditive();

			for (int i = 0; i < (int)m_NumTotal; ++i)
			{
				m_PhyxelExtForces[i] = Vector3(0, 0, 0);
				m_PhyxelsVel[i] = m_Solver.m_X[i];
				m_PhyxelsPosTemp[i] = m_PhyxelsPos[i] + m_PhyxelsVel[i] * m_TimeStep;
			}

			valid_timestep = true;// ValidateVelocityTimestep();
								  /*if (!valid_timestep)
								  {
								  float new_timestep = max(min_timestep, m_TimeStep / 2.0f);
								  if (new_timestep != m_TimeStep)
								  {
								  printf("Unstable Timestep\t%f\t->\t%f\n", m_TimeStep, new_timestep);
								  m_TimeStep = new_timestep;
								  m_StepCounter = 0;
								  }
								  else
								  {
								  printf("ERROR: UNSTABLE CLOTH! Smallest timestep is not small enough for simulation.\n\n");
								  valid_timestep = true;
								  }
								  }*/
		}

		timestep_accum += m_TimeStep;
		angle += m_TimeStep;

		//Reset All Accumulative Data & Update Positions	
		memcpy(&m_PhyxelsPos[0].x, &m_PhyxelsPosTemp[0].x, m_NumTotal * sizeof(Vector3));
	}

	//m_StepCounter++;
	//if (m_StepCounter == m_StepsBeforeIncrease)
	//{
	//	float new_timestep = min(max_timestep, m_TimeStep * 2.0f);
	//	if (new_timestep != m_TimeStep)
	//	{
	//		printf("Trial Timestep   \t%f\t->\t%f\n", m_TimeStep, new_timestep);
	//		m_TimeStep = new_timestep;
	//		m_StepCounter = 0;
	//		m_StepsBeforeIncrease += m_StepsBeforeIncrease * 0.1f;
	//		m_StepsBeforeIncrease = min(m_StepsBeforeIncrease, 120);
	//	}
	//}

	m_ProfilingTotalTime.EndTiming();
}


bool FEM3NodedLMA::ValidateVelocityTimestep()
{
	//Check whether any of the edges of the triangles change by more than 10%
	for (uint i = 0; i < m_NumTriangles; ++i)
	{
		uint idxA = m_Triangles[i].phyxels[0];
		uint idxB = m_Triangles[i].phyxels[1];
		uint idxC = m_Triangles[i].phyxels[2];

		Vector3 ABnew = m_PhyxelsPosTemp[idxB] - m_PhyxelsPosTemp[idxA];
		Vector3 ACnew = m_PhyxelsPosTemp[idxC] - m_PhyxelsPosTemp[idxA];
		Vector3 BCnew = m_PhyxelsPosTemp[idxC] - m_PhyxelsPosTemp[idxB];

		Vector3 ABold = m_PhyxelsPos[idxB] - m_PhyxelsPos[idxA];
		Vector3 ACold = m_PhyxelsPos[idxC] - m_PhyxelsPos[idxA];
		Vector3 BCold = m_PhyxelsPos[idxC] - m_PhyxelsPos[idxB];

		float ABLenDif = 1.0f - ABnew.Length() / ABold.Length();
		float ACLenDif = 1.0f - ACnew.Length() / ACold.Length();
		float BCLenDif = 1.0f - BCnew.Length() / BCold.Length();

		const float tolerence = 0.1f;
		if (fabs(ABLenDif) > tolerence) return false;
		if (fabs(ACLenDif) > tolerence) return false;
		if (fabs(BCLenDif) > tolerence) return false;
	}

	return true;
}

void FEM3NodedLMA::CollideEllipsoid(const Ellipsoid& e)
{
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTotal; i++) {
		if (m_PhyxelIsStatic[i])
			continue;

		Vector4 X_0 = (e.invTransform * Vector4(m_PhyxelsPos[i], 1));
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


			//m_PhyxelsPos[i] += delta;
			m_PhyxelsPos[i] += delta;

			Vector3 col_norm = delta; col_norm.Normalise();
			if (Vector3::Dot(col_norm, m_PhyxelsVel[i]) < 0.0f)
			{

				m_PhyxelsVel[i] -= col_norm * Vector3::Dot(m_PhyxelsVel[i], col_norm);

				m_Solver.m_X[i] = m_PhyxelsVel[i];
				m_Solver.m_Constraints[i] = Matrix3::Identity - Matrix3::OuterProduct(col_norm, col_norm);// Matrix3::OuterProduct(tan1, tan1) + Matrix3::OuterProduct(tan2, tan2);
			}
		}
	}
}

void FEM3NodedLMA::GenFEMTriangle(unsigned int idx, FETriangle& tri)
{
	const unsigned int idxA = tri.phyxels[0];
	const unsigned int idxB = tri.phyxels[1];
	const unsigned int idxC = tri.phyxels[2];

	Vector3 a = m_PhyxelsPosInitial[idxA];
	Vector3 b = m_PhyxelsPosInitial[idxB];
	Vector3 c = m_PhyxelsPosInitial[idxC];

	Vector3 e1 = a - c;
	Vector3 e2 = b - c;

	//Area + Mass of Triangle/Phyxels
	tri.Area = 0.5f * abs(e1.x * e2.y - e1.y * e2.x);

	//Rotation
	tri.rotBase = Matrix3::Transpose(BuildRotationMatrix(a, b, c));

	//Computes Shape Function and Derivitives (Assumes triangle only exists in the x/y planes)
	Vector2 r1 = Vector2(e1.x, e2.x);
	Vector2 r2 = Vector2(e1.y, e2.y);

	float detJ = 1.f / (r1.x * r2.y - r1.y * r2.x);

	Vector2 invJ1 = Vector2(r2.y, -r1.y) * detJ;
	Vector2 invJ2 = Vector2(-r2.x, r1.x) * detJ;
	Vector2 invJ3 = Vector2(0, 0) - invJ1 - invJ2;

	Vector3 dN[3]{
		Vector3(invJ1.x, invJ1.y, 0.f),
		Vector3(invJ2.x, invJ2.y, 0.f),
		Vector3(invJ3.x, invJ3.y, 0.f)
	};

	/*float det = Vector3::Dot(e1, e1) * Vector3::Dot(e2, e2) - Vector3::Dot(e1, e2) * Vector3::Dot(e1, e2);
	tri.dN[0] = (e1 * Vector3::Dot(e2, e2) - e2 * Vector3::Dot(e1, e2)) / det;
	tri.dN[1] = (e2 * Vector3::Dot(e1, e1) - e1 * Vector3::Dot(e1, e2)) / det;
	tri.dN[2] = -tri.dN[0] - tri.dN[1];*/


	//Stiffness Matrix for each vertex pair
	for (int i_itr = 0; i_itr < 3; ++i_itr)
	{
		uint i = i_itr;
		uint j = (i + 1) % 3;

		const Vector3& dNa = dN[i];
		const Vector3& dNb = dN[j];

		Matrix3& Ke = tri.Ke[i];

		Ke = BuildStiffnessMatrix(dNa, dNb, tri.Area);
	}

	//Bending Matrix for Each Vertex Pair
	for (int i = 0; i < 3; ++i)
	{
		int j = (i + 1) % 3;

		Matrix3& B = tri.B[i];
		const Vector3& dNa = dN[i];
		const Vector3& dNb = dN[j];

		B(0, 0) = dNa.x * dNb.x * B1 * tri.Area;
		B(1, 1) = dNa.y * dNb.y * B2 * tri.Area;
		B(2, 2) = 0.0f;
	}
}

Matrix3 FEM3NodedLMA::BuildStiffnessMatrix(const Vector3& dNa, const Vector3& dNb, float triArea)
{
	Matrix3 Ke = Matrix3::ZeroMatrix;

	//Cloth Method
	/*Ke(0, 0) = (dNa[0] * C1111 * dNb[0])
	+ (dNa[0] * C1211 * dNb[0])
	+ (dNa[0] * C1112 * dNb[1])
	+ (dNa[0] * C1212 * dNb[1]);

	Ke(0, 1) = (dNa[0] * C1121 * dNb[0])
	+ (dNa[0] * C1221 * dNb[0])
	+ (dNa[0] * C1122 * dNb[1])
	+ (dNa[0] * C1222 * dNb[1]);

	Ke(1, 0) = (dNa[1] * C2111 * dNb[0])
	+ (dNa[1] * C2211 * dNb[0])
	+ (dNa[1] * C2112 * dNb[1])
	+ (dNa[1] * C2212 * dNb[1]);

	Ke(1, 1) = (dNa[1] * C2121 * dNb[0])
	+ (dNa[1] * C2221 * dNb[0])
	+ (dNa[1] * C2122 * dNb[1])
	+ (dNa[1] * C2222 * dNb[1]);*/


	/*Ke(0, 0) = (dNa[0] * dNb[0] * C1111)
	//+ (dNa[1] * dNb[0] * C1121)
	//+ (dNa[0] * dNb[1] * C1112)
	+ (dNa[1] * dNb[1] * C1122);

	Ke(1, 0) = 0.0f//(dNa[0] * dNb[0] * C2111)
	+ (dNa[1] * dNb[0] * C2121)
	+ (dNa[0] * dNb[1] * C2112);
	//+ (dNa[1] * dNb[1] * C2122);

	Ke(0, 1) = 0.0f//(dNa[0] * dNb[0] * C1211)
	+ (dNa[1] * dNb[0] * C1221)
	+ (dNa[0] * dNb[1] * C1212);
	//+ (dNa[1] * dNb[1] * C1222);

	Ke(1, 1) = (dNa[0] * dNb[0] * C2211)
	//+ (dNa[1] * dNb[0] * C2221)
	//+ (dNa[0] * dNb[1] * C2212)
	+ (dNa[1] * dNb[1] * C2222);*/


	/*float ab = dNa[0] * dNb[1] + dNa[1] * dNb[0];
	ab *= 0.5f;

	float trace = (dNa[0] * dNb[0]) + (dNa[1] * dNb[1]);
	Ke._11 = (dNa[0] * C1111 * dNb[0])
	//+ (dNa[0] * C1212 * dNb[1])
	+ab * C1122;
	//+ trace * C1111;

	Ke._12 = ab * C1122;
	Ke._21 = ab * C1122;

	Ke._22 = (dNa[1] * C2222 * dNb[1])
	//+ (dNa[1] * C2222 * dNb[1])
	+ ab * C1122;
	//+ trace * C2222;*/



	Ke(0, 0) = (dNa[0] * C1111 * dNb[0])
		+ (dNa[0] * C1211 * dNb[0])
		+ (dNa[0] * C1112 * dNb[1])
		+ (dNa[0] * C1212 * dNb[1]);

	Ke(0, 1) = (dNa[0] * C1121 * dNb[0])
		+ (dNa[0] * C1221 * dNb[0])
		+ (dNa[0] * C1122 * dNb[1])
		+ (dNa[0] * C1222 * dNb[1]);

	Ke(1, 0) = (dNa[1] * C2111 * dNb[0])
		+ (dNa[1] * C2211 * dNb[0])
		+ (dNa[1] * C2112 * dNb[1])
		+ (dNa[1] * C2212 * dNb[1]);

	Ke(1, 1) = (dNa[1] * C2121 * dNb[0])
		+ (dNa[1] * C2221 * dNb[0])
		+ (dNa[1] * C2122 * dNb[1])
		+ (dNa[1] * C2222 * dNb[1]);

	float trace = (dNa[0] * dNb[0] + dNa[1] * dNb[1]);
	Ke(0, 0) += trace * C2211 * 0.5f;
	Ke(1, 1) += trace * C2211 * 0.5f;

	Matrix3 E;
	E(0, 0) = 756;
	E(0, 1) = 311;
	E(0, 2) = 0;
	E(1, 0) = 460;
	E(1, 1) = 975;
	E(1, 2) = 0;
	E(2, 0) = 0;
	E(2, 1) = 0;
	E(2, 2) = 30;

	Matrix3 dnab = Matrix3::OuterProduct(dNa, dNb);

	Ke = dnab * E * Matrix3::Transpose(dnab);


	//Ke = Ke * 0.5f;
	float lambda = 300.f;// 1366.f;
	float mu = 300.f;// 800.f;

	Matrix3 tmp2 = (Matrix3::OuterProduct(dNa, dNb) + Matrix3::OuterProduct(dNb, dNa)) * 0.5f;
	//trace = tmp2.Trace();
	/*	Ke = tmp2 * (mu * 2.0f);

	float mu2 = mu * 2.0f;
	Ke._11 = (dNa[0] * C1111 * dNb[0]);
	Ke._21 = (dNa[1] * C2211 * dNb[0]);
	Ke._12 = (dNa[0] * C2211 * dNb[1]);
	Ke._22 = (dNa[1] * C2222 * dNb[1]);

	Ke._11 += trace * C2211 * 0.5f;
	Ke._22 += trace * C2211 * 0.5f;
	*/
	Ke = (Matrix3::Identity * (lambda * tmp2.Trace())) + tmp2 * (mu * 2.0f);
	Ke(2, 2) = 0.f;
	Ke *= triArea;
	//Ke = (Ke + Matrix3::Transpose(Ke)) * 0.5f;

	/*Matrix3 tmp = (Matrix3::OuterProduct(dNa, dNb) + Matrix3::OuterProduct(dNb, dNa)) * 0.5f;
	Ke = tmp;
	Ke(2, 2) = 0.0f;

	trace = tmp.Trace() * 0.1f;
	Ke._12 *= C1122;
	Ke._21 *= C1122;

	Ke._11 = trace * C2222 + Ke._11 * C1111;
	Ke._22 = trace * C2222 + Ke._22 * C2222;

	Ke *= triArea;*/

	//Mat33 meh = (Ke + Ke.transpose()) * 0.5f;
	//Mat33 meh2 = (Mat33::Identity() * (lambda * tmp.trace()) + tmp * (mu * 2.0f)) * triArea;
	return Ke;
}

void FEM3NodedLMA::BuildRotationAndNormals()
{

	memset(&m_PhyxelNormals[0], 0, m_NumTotal * sizeof(Vec3));


	//Build Rotations and Sum Up Normals
#pragma omp parallel for
	for (int i = 0; i < m_NumTriangles; ++i)
	{
		FETriangle& t = m_Triangles[i];
		int idxA = t.phyxels[0];
		int idxB = t.phyxels[1];
		int idxC = t.phyxels[2];

		Vector3 a = m_PhyxelsPos[idxA];
		Vector3 b = m_PhyxelsPos[idxB];
		Vector3 c = m_PhyxelsPos[idxC];

		t.R = BuildRotationMatrix(a, b, c) * t.rotBase;
		//R2[i] = BuildRotationMatrix33(a, b, c) * rotBase2[i].transpose();

	}

#if !USE_TRIANGLE_NORMALS
	for (const Quad& quad : m_Quads)
	{
		Vector3 a = m_PhyxelsPos[quad.vert_a];
		Vector3 b = m_PhyxelsPos[quad.vert_b];
		Vector3 c = m_PhyxelsPos[quad.vert_c];
		Vector3 d = m_PhyxelsPos[quad.vert_d];

		Vector3 normalA = Vector3::Cross(a - c, b - c);
		Vector3 normalB = Vector3::Cross(d - c, a - c);
		normalA.Normalise();
		normalB.Normalise();

		Vector3 normal = normalA + normalB;
		normal.Normalise();

		m_PhyxelNormals[quad.vert_a] += normal;
		m_PhyxelNormals[quad.vert_b] += normal;
		m_PhyxelNormals[quad.vert_c] += normal;
		m_PhyxelNormals[quad.vert_d] += normal;
	}
#endif
	//Normalise Vertex Normals
#pragma omp parallel for
	for (int i = 0; i < m_NumTotal; ++i)
	{
		m_PhyxelNormals[i].Normalise();
	}
}

Matrix3 FEM3NodedLMA::BuildRotationMatrix(const Vector3& a, const Vector3& b, const Vector3& c)
{
	Vector3 t_rows[3];
	t_rows[0] = a - c;
	t_rows[1] = b - c;
	t_rows[2] = Vector3::Cross(t_rows[0], t_rows[1]);

	//Build Rotation Matrtix
	t_rows[0].Normalise();
	t_rows[1] -= t_rows[0] * Vector3::Dot(t_rows[0], t_rows[1]);
	t_rows[2] -= t_rows[0] * Vector3::Dot(t_rows[0], t_rows[2]);

	t_rows[1].Normalise();
	t_rows[2] -= t_rows[1] * Vector3::Dot(t_rows[1], t_rows[2]);

	t_rows[2].Normalise();

	Matrix3 rot(t_rows[0], t_rows[1], t_rows[2]);
	return rot;
}

void FEM3NodedLMA::SimpleCorotatedBuildAMatrix(float dt)
{
#define LOCK_FREE_AMATRIX 0

#pragma omp parallel
{
	const float dt2 = dt * dt;
	const float dt2_visos = dt2 + V_SCALAR * dt;


#pragma omp for
	for (int i = 0; i < (int)m_NumTotal; ++i)
	{
		m_Solver.m_B[i] = m_PhyxelsVel[i] * m_PhyxelsMass[i] + m_PhyxelForces[i] * dt;
	}




	Matrix3 RK, RKR, B_rot, Q, ReT;
	Vector3 posAB, posIAB;

#pragma omp for private(RK, RKR, ReT, B_rot, Q, posAB, posIAB)
	for (int i = 0; i < (int)m_NumTriangles; ++i)
	{
		const FETriangle& t = m_Triangles[i];
		uint idxA = t.phyxels[0];
		uint idxB = t.phyxels[1];
		uint idxC = t.phyxels[2];

		Matrix3* A_ii = m_Solver.m_A[i].Kii;
		Matrix3* A_ij = m_Solver.m_A[i].Kij;

		memset(A_ii, 0, 3 * sizeof(Matrix3));
		memset(A_ij, 0, 3 * sizeof(Matrix3));

		Vector3 B_i[3];



		Vector3 triNormals[3] = { m_PhyxelNormals[idxA], m_PhyxelNormals[idxB], m_PhyxelNormals[idxC] };

		Matrix3 Pa = Matrix3::OuterProduct(triNormals[0], triNormals[0]);
		Matrix3 Pb = Matrix3::OuterProduct(triNormals[1], triNormals[1]);
		Matrix3 Pc = Matrix3::OuterProduct(triNormals[2], triNormals[2]);

		Vector3 face_normal = (triNormals[0] + triNormals[1] + triNormals[2]) / 3.0f;
		face_normal.Normalise();

		//Air Drag
		{
			const float Cd = 0.55f; //Drag Coefficient
			const float Cl = 0.55f; //Lift Coefficient
			const float p = 1.293f; //Air Density

									//Vrel = Vi - u (Where u is velocity field of the wind)
									//F = 0.5f * Cd * p * |Vrel|^2 * A * (n . Vrel) . (-Vrel)
			Vector3 vrel = (m_PhyxelsVel[idxA] + m_PhyxelsVel[idxB] + m_PhyxelsVel[idxC]) / 3.0f - WIND;
			Vector3 vrel_norm = vrel; vrel_norm.Normalise();

			if (Vector3::Dot(face_normal, vrel) < 0) face_normal = -face_normal;

			Vector3 lift_dir = Vector3::Cross(Vector3::Cross(face_normal, vrel_norm), vrel_norm);

			float vrelsq = Vector3::Dot(vrel, vrel);
			float drag_coef = 0.5f * Cd * p * vrelsq * t.Area * Vector3::Dot(face_normal, vrel_norm);
			float lift_coef = 0.5f * Cl * p * vrelsq * t.Area * Vector3::Dot(lift_dir, vrel_norm);

			Vector3 Fid = (-vrel_norm * drag_coef) / 3.0f;		//Drag
			Vector3 Fil = (lift_dir * lift_coef) / 3.0f;		//Lift

			Vector3 totalForce = (Fid + Fil);
			B_i[0] = totalForce;
			B_i[1] = totalForce;
			B_i[2] = totalForce;
		}




		Vector3 tmpMult;
		for (uint j = 0; j < 3; ++j)
		{
			uint k = (j + 1) % 3;

			const int idxA = t.phyxels[j];
			const int idxB = t.phyxels[k];

			const Matrix3& ke = t.Ke[j];
			const Matrix3& Re = t.R;
			ReT = Matrix3::Transpose(Re);


			posAB = m_PhyxelsPos[idxA] - m_PhyxelsPos[idxB];
			posIAB = m_PhyxelsPosInitial[idxA] - m_PhyxelsPosInitial[idxB];

			//RK = Re * ke;
			//RKR = Re * ke * ReT;
			InplaceMatrix3MultMatrix3(&RK, Re, ke);
			InplaceMatrix3MultMatrix3(&RKR, RK, ReT);

			//tmpMult = RKR * posAB;
			InplaceMatrix3MultVector3(&tmpMult, RKR, posAB);
			B_i[j] += tmpMult;
			B_i[k] -= tmpMult;

			InplaceMatrix3MultVector3(&tmpMult, RK, posIAB);
			B_i[j] -= tmpMult;
			B_i[k] += tmpMult;

			RKR *= dt2_visos;
			A_ii[j] -= RKR;
			A_ii[k] -= RKR;
			A_ij[j] += RKR;


			//Bending			
			/*B_rot = Re * t.B[j] * ReT;
			Q = Pa * B_rot * Pa
			+ Pb * B_rot * Pb
			+ Pc * B_rot * Pc;
			Q /= 3.f;*/
			InplaceMatrix3MultMatrix3(&RK, Re, t.B[j]);
			InplaceMatrix3MultMatrix3(&B_rot, RK, ReT);

			InplaceMatrix3MultMatrix3(&RK, Pa, B_rot);
			InplaceMatrix3MultMatrix3(&Q, RK, Pa);

			InplaceMatrix3MultMatrix3(&RK, Pb, B_rot);
			InplaceMatrix3MultMatrix3Additive(&Q, RK, Pb);
			InplaceMatrix3MultMatrix3(&RK, Pc, B_rot);
			InplaceMatrix3MultMatrix3Additive(&Q, RK, Pc);
			Q *= (1.0f / 3.0f);

			InplaceMatrix3MultVector3(&tmpMult, Q, posAB);
			B_i[j] += tmpMult;
			B_i[k] -= tmpMult;


			Q *= dt2_visos;
			A_ii[j] -= Q;
			A_ii[k] -= Q;
			A_ij[j] += Q;
		}


		B_i[0] *= dt;
		B_i[1] *= dt;
		B_i[2] *= dt;

		int i3 = i * 3;
		m_TrianglesB[i3 + 0] = B_i[0];
		m_TrianglesB[i3 + 1] = B_i[1];
		m_TrianglesB[i3 + 2] = B_i[2];
	}

#pragma omp for
	for (int i = 0; i < (int)m_NumTotal; ++i)
	{
		float mass = m_PhyxelsMass[i];

		Vector3 BVec = m_PhyxelsVel[i] * mass + m_PhyxelForces[i] * dt;
		Matrix3 Aii = Matrix3(mass, 0.0f, 0.0f,
			0.0f, mass, 0.0f,
			0.0f, 0.0f, mass);

		const auto& tris = m_PhyxelsTriangleLookup[i];
		for (uint itr = 0, end = tris.size(); itr < end; ++itr)
		{
			uint idx = tris[itr];

			BVec += m_TrianglesB[idx];
			Aii += m_Solver.m_A[idx / 3].Kii[idx % 3];
		}

		m_Solver.m_PreCondition[i] = Matrix3::Inverse(Aii);
		m_Solver.m_B[i] = BVec;
	}
	
}
}

void FEM3NodedLMA::Render_DrawingToVisualiser()
{
	/*unsigned int i;
	for (i = 0; i < m_NumTriangles; ++i)
	{
	const FETriangle& t = m_Triangles[i];

	Vector3 a = m_PhyxelsPos[t.phyxels[0]];
	Vector3 b = m_PhyxelsPos[t.phyxels[1]];
	Vector3 c = m_PhyxelsPos[t.phyxels[2]];

	Vector3 centre = (a + b + c) * 0.3333333f;

	//	memcpy(&rotation._11, &t.R(0,0), 9 * sizeof(float));
	//NCLDebug::DrawMatrix(t.R, centre, 0.02f);

	NCLDebug::DrawThickLine(a, b, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	NCLDebug::DrawThickLine(b, c, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	NCLDebug::DrawThickLine(a, c, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	}*/
}