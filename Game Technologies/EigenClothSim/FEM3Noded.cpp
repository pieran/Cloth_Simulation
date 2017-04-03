#include "FEM3Noded.h"
#include <ncltech\NCLDebug.h>




FEM3Noded::FEM3Noded()
{
	m_TimeStep = max_timestep;
	m_StepCounter = 0;
	m_StepsBeforeIncrease = 10;
}

FEM3Noded::~FEM3Noded()
{
}

void FEM3Noded::simulation_OnClothDesignChanged(ClothDesignEntity<3>* design)
{

	//Gen Vertices
	const std::vector<Vertex>& vertices = design->Vertices();
	m_NumTotal = vertices.size();

	m_Solver.AllocateMemory(m_NumTotal);
	m_Solver.m_A.resize(m_NumTotal);

	m_PhyxelsPos.resize(m_NumTotal);
	m_PhyxelsPosTemp.resize(m_NumTotal);
	m_PhyxelsPosInitial.resize(m_NumTotal);
	m_PhyxelsVel.resize(m_NumTotal);
	m_PhyxelsVelChange.resize(m_NumTotal);
	m_PhyxelForces.resize(m_NumTotal);
	m_PhyxelExtForces.resize(m_NumTotal);

	m_RenderValue.resize(m_NumTotal);
	m_PhyxelNormals.resize(m_NumTotal);
	m_PhyxelIsStatic.resize(m_NumTotal);
	m_PhyxelsMass.resize(m_NumTotal);
	m_PhyxelTexCoords.resize(m_NumTotal);

	memset(&m_PhyxelsVel[0].x, 0, m_NumTotal * sizeof(Vector3));
	memset(&m_PhyxelsVelChange[0].x, 0, m_NumTotal * sizeof(Vector3));
	memset(&m_PhyxelExtForces[0].x, 0, m_NumTotal * sizeof(Vector3));
	memset(&m_PhyxelNormals[0].x, 0, m_NumTotal * sizeof(Vector3));

	m_RenderValueMinMax.x = 0.0f;
	m_RenderValueMinMax.y = 500.5f;

	for (unsigned int i = 0; i < m_NumTotal; ++i)
	{
		m_RenderValue[i] = (float(i * 5.0f) / float(m_NumTotal));// ((rand() % 100) / 100.0f) * 5.0f;

		const Vertex& v = vertices[i];

		m_PhyxelsPos[i] = Vector3(v.position.x(), v.position.y(), v.position.z());
		m_PhyxelsPosInitial[i] = Vector3(v.initial_position.x(), v.initial_position.y(), v.initial_position.z());
		m_PhyxelIsStatic[i] = v.flags & VERTEXFLAGS_IS_STATIC;

		if (!m_PhyxelIsStatic[i])
		{
			m_PhyxelForces[i] = Vector3(0.f, 0.f, 0.f);// GRAVITY / float(m_NumTotal);
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

		GenFEMTriangle(i, m_Triangles[i]);

		totalArea += m_Triangles[i].Area;
	}

	float uniform_mass = (totalArea * mass_density) / m_NumTotal;
	for (unsigned int i = 0; i < m_NumTotal; ++i)
	{
		m_PhyxelsMass[i] = uniform_mass;
	}

	//Build Global Matricies (To be deleted at the start of the simulation)
	SimpleCorotatedBuildAMatrix(0.0f);
	m_Solver.m_A.zero_memory();



	for (uint i = 0; i < m_NumTotal; ++i)
	{
		m_PhyxelExtForces[i] = Vector3(0, 0, 0);
	}
	m_Solver.ResetMemory();
	m_Solver.m_A.zero_memory();
	UpdateConstraints();
}

void  FEM3Noded::UpdateConstraints()
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

void FEM3Noded::Simulation_StepSimulation(float dt)
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
				float spd = (i < (m_NumTotal >> 1)) ? 0.85f : 0.0f;

				//m_PhyxelsVel[i].z = sin(angle) * 0.01f;// *(static_itr < 33 ? 1.0f : -1.0f);
				m_PhyxelsVel[i].y = sin(angle) * spd;// *(static_itr++ < 33 ? 1.0f : -1.0f);

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
			m_Solver.m_A.zero_memory();
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


bool FEM3Noded::ValidateVelocityTimestep()
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

void FEM3Noded::CollideEllipsoid(const Ellipsoid& e)
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

void FEM3Noded::GenFEMTriangle(unsigned int idx, FETriangleLite& tri)
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
}




// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0
typedef double Real;
void Diagonalize(const Matrix3&  A, Matrix3& Q, Matrix3& D)
{
	// A must be a symmetric matrix.
	// returns Q and D such that 
	// Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
	const int maxsteps = 24;  // certainly wont need that many.
	int k0, k1, k2;
	Real o[3], m[3];
	Real q[4] = { 0.0,0.0,0.0,1.0 };
	Real jr[4];
	Real sqw, sqx, sqy, sqz;
	Real tmp1, tmp2, mq;
	Real AQ[3][3];
	Real thet, sgn, t, c;
	for (int i = 0; i < maxsteps; ++i)
	{
		// quat to matrix
		sqx = q[0] * q[0];
		sqy = q[1] * q[1];
		sqz = q[2] * q[2];
		sqw = q[3] * q[3];
		Q[0] = (sqx - sqy - sqz + sqw);
		Q[3 + 1] = (-sqx + sqy - sqz + sqw);
		Q[6 + 2] = (-sqx - sqy + sqz + sqw);
		tmp1 = q[0] * q[1];
		tmp2 = q[2] * q[3];
		Q[3 + 0] = 2.0 * (tmp1 + tmp2);
		Q[0 + 1] = 2.0 * (tmp1 - tmp2);
		tmp1 = q[0] * q[2];
		tmp2 = q[1] * q[3];
		Q[6 + 0] = 2.0 * (tmp1 - tmp2);
		Q[0 + 2] = 2.0 * (tmp1 + tmp2);
		tmp1 = q[1] * q[2];
		tmp2 = q[0] * q[3];
		Q[6 + 1] = 2.0 * (tmp1 + tmp2);
		Q[3 + 2] = 2.0 * (tmp1 - tmp2);

		// AQ = A * Q
		AQ[0][0] = Q[0 + 0] * A[0 + 0] + Q[3 + 0] * A[0 + 1] + Q[6 + 0] * A[0 + 2];
		AQ[0][1] = Q[0 + 1] * A[0 + 0] + Q[3 + 1] * A[0 + 1] + Q[6 + 1] * A[0 + 2];
		AQ[0][2] = Q[0 + 2] * A[0 + 0] + Q[3 + 2] * A[0 + 1] + Q[6 + 2] * A[0 + 2];
		AQ[1][0] = Q[0 + 0] * A[0 + 1] + Q[3 + 0] * A[3 + 1] + Q[6 + 0] * A[3 + 2];
		AQ[1][1] = Q[0 + 1] * A[0 + 1] + Q[3 + 1] * A[3 + 1] + Q[6 + 1] * A[3 + 2];
		AQ[1][2] = Q[0 + 2] * A[0 + 1] + Q[3 + 2] * A[3 + 1] + Q[6 + 2] * A[3 + 2];
		AQ[2][0] = Q[0 + 0] * A[0 + 2] + Q[3 + 0] * A[3 + 2] + Q[6 + 0] * A[6 + 2];
		AQ[2][1] = Q[0 + 1] * A[0 + 2] + Q[3 + 1] * A[3 + 2] + Q[6 + 1] * A[6 + 2];
		AQ[2][2] = Q[0 + 2] * A[0 + 2] + Q[3 + 2] * A[3 + 2] + Q[6 + 2] * A[6 + 2];
		// D = Qt * AQ
		D[0 + 0] = AQ[0][0] * Q[0 + 0] + AQ[1][0] * Q[3 + 0] + AQ[2][0] * Q[6 + 0];
		D[0 + 1] = AQ[0][0] * Q[0 + 1] + AQ[1][0] * Q[3 + 1] + AQ[2][0] * Q[6 + 1];
		D[0 + 2] = AQ[0][0] * Q[0 + 2] + AQ[1][0] * Q[3 + 2] + AQ[2][0] * Q[6 + 2];
		D[3 + 0] = AQ[0][1] * Q[0 + 0] + AQ[1][1] * Q[3 + 0] + AQ[2][1] * Q[6 + 0];
		D[3 + 1] = AQ[0][1] * Q[0 + 1] + AQ[1][1] * Q[3 + 1] + AQ[2][1] * Q[6 + 1];
		D[3 + 2] = AQ[0][1] * Q[0 + 2] + AQ[1][1] * Q[3 + 2] + AQ[2][1] * Q[6 + 2];
		D[6 + 0] = AQ[0][2] * Q[0 + 0] + AQ[1][2] * Q[3 + 0] + AQ[2][2] * Q[6 + 0];
		D[6 + 1] = AQ[0][2] * Q[0 + 1] + AQ[1][2] * Q[3 + 1] + AQ[2][2] * Q[6 + 1];
		D[6 + 2] = AQ[0][2] * Q[0 + 2] + AQ[1][2] * Q[3 + 2] + AQ[2][2] * Q[6 + 2];
		o[0] = D[3 + 2];
		o[1] = D[0 + 2];
		o[2] = D[0 + 1];
		m[0] = fabs(o[0]);
		m[1] = fabs(o[1]);
		m[2] = fabs(o[2]);

		k0 = (m[0] > m[1] && m[0] > m[2]) ? 0 : (m[1] > m[2]) ? 1 : 2; // index of largest element of offdiag
		k1 = (k0 + 1) % 3;
		k2 = (k0 + 2) % 3;
		if (o[k0] == 0.0)
		{
			break;  // diagonal already
		}
		thet = (D[k2 * 3 + k2] - D[k1 * 3 + k1]) / (2.0*o[k0]);
		sgn = (thet > 0.0) ? 1.0 : -1.0;
		thet *= sgn; // make it positive
		t = sgn / (thet + ((thet < 1.E6) ? sqrt(thet*thet + 1.0) : thet)); // sign(T)/(|T|+sqrt(T^2+1))
		c = 1.0 / sqrt(t*t + 1.0); //  c= 1/(t^2+1) , t=s/c 
		if (c == 1.0)
		{
			break;  // no room for improvement - reached machine precision.
		}
		jr[0] = jr[1] = jr[2] = jr[3] = 0.0;
		jr[k0] = sgn*sqrt((1.0 - c) / 2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)  
		jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
		jr[3] = sqrt(1.0f - jr[k0] * jr[k0]);
		if (jr[3] == 1.0)
		{
			break; // reached limits of floating point precision
		}
		q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
		q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
		q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
		q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
		mq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		q[0] /= mq;
		q[1] /= mq;
		q[2] /= mq;
		q[3] /= mq;
	}
}

void FEM3Noded::BuildRotationAndNormals()
{

	memset(&m_PhyxelNormals[0], 0, m_NumTotal * sizeof(Vec3));


	//Build Rotations and Sum Up Normals
#pragma omp parallel for
	for (int i = 0; i < m_NumTriangles; ++i)
	{
		FETriangleLite& t = m_Triangles[i];
		int idxA = t.phyxels[0];
		int idxB = t.phyxels[1];
		int idxC = t.phyxels[2];

		Vector3 a = m_PhyxelsPos[idxA];
		Vector3 b = m_PhyxelsPos[idxB];
		Vector3 c = m_PhyxelsPos[idxC];

		t.R = BuildRotationMatrix(a, b, c) * t.rotBase;

		//SHAPE MATCHING!!!
		Vector3 pos[3]{
			m_PhyxelsPos[idxA],
			m_PhyxelsPos[idxB],
			m_PhyxelsPos[idxC]
		};

		Vector3 posi[3]{
			m_PhyxelsPosInitial[idxA],
			m_PhyxelsPosInitial[idxB],
			m_PhyxelsPosInitial[idxC]
		};

		Vector3 center = (pos[0] + pos[1] + pos[2]) / 3.0f;
		Vector3 centeri = (posi[0] + posi[1] + posi[2]) / 3.0f;


		Matrix3 Apq = Matrix3::ZeroMatrix;
		for (int i = 0; i < 3; ++i)
		{
			Apq += Matrix3::OuterProduct(pos[i] - center, posi[i] - centeri);
		}
		Apq /= 3.0f;
		Apq._33 = 1.0f;
		

		Matrix3 S2 = Matrix3::Transpose(Apq) * Apq;

		Matrix3 Q, D;
		Diagonalize(S2, Q, D);
		const float epsilon = 1E-6f;
		Matrix3 sD(sqrtf(D._11), 0.0f, 0.0f,
			0.0f, sqrtf(D._22), 0.0f,
			0.0f, 0.0f, sqrtf(D._33));

		Matrix3 S = Matrix3::Transpose(Q) * sD * Q;
		Matrix3 rot = Apq * Matrix3::Inverse(S);
		
		Quaternion q_apq = Quaternion::FromMatrix(Matrix3::Inverse(Apq));
		q_apq.Normalise();
		Matrix3 rot2 = q_apq.ToMatrix3();

		if (Window::GetKeyboard()->KeyDown(KEYBOARD_B))
		t.R = rot;
	//	printf("moo");

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

Matrix3 FEM3Noded::BuildRotationMatrix(const Vector3& a, const Vector3& b, const Vector3& c)
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

void FEM3Noded::SimpleCorotatedBuildAMatrix(float dt)
{
	const float dt2 = dt * dt;
	const float dt2_viscos = dt2 + V_SCALAR * dt;


//#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTotal; ++i)
	{
		m_Solver.m_A(i, i) = Matrix3::Identity * m_PhyxelsMass[i];
		m_Solver.m_B[i] = m_PhyxelsVel[i] * m_PhyxelsMass[i] +m_PhyxelForces[i] * dt;
	}


	Matrix3 RK, RKR, ReT;
	Vector3 posAB, posIAB;


	for (int i = 0; i < (int)m_NumTriangles; ++i)
	{
		const FETriangleLite& t = m_Triangles[i];
		uint idxA = t.phyxels[0];
		uint idxB = t.phyxels[1];
		uint idxC = t.phyxels[2];


		Vector3 triNormals[3] = { m_PhyxelNormals[idxA], m_PhyxelNormals[idxB], m_PhyxelNormals[idxC] };

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
			m_Solver.m_B[idxA] += totalForce * dt;
			m_Solver.m_B[idxB] += totalForce * dt;
			m_Solver.m_B[idxC] += totalForce * dt;
		}





		Vector3 a = m_PhyxelsPosInitial[idxA];
		Vector3 b = m_PhyxelsPosInitial[idxB];
		Vector3 c = m_PhyxelsPosInitial[idxC];

		Vector3 e1 = a - c;
		Vector3 e2 = b - c;
		Vector3 e3 = a - b;

		float detJ = e1.x * e2.y - e1.y - e2.x;
		Eigen::MatrixXf B(3, 6);
			B.setZero();
			B(0, 0) = e2.y;
			B(2, 0) = -e2.x;
			B(1, 1) = -e2.x;
			B(2, 1) = e2.y;
			B(0, 2) = -e1.y;
			B(2, 2) = e1.x;
			B(1, 3) = e1.x;
			B(2, 3) = -e1.y;
			B(0, 4) = e3.y;
			B(2, 4) = -e3.x;
			B(1, 5) = -e3.x;
			B(2, 5) = e3.y;
			B *= 1.0f / detJ;

		float area = 0.5f * detJ;


		const float E = 1000.0f;		//Youngs Modulus
		const float v = 0.3f;		//Poisson coefficient
		Eigen::Matrix3f D;
			D.setZero();
			D(0, 0) = 1.0f;
			D(1, 0) = v;
			D(0, 1) = v;
			D(1, 1) = 1.0f;
			D(2, 2) = (1.0f - v) * 0.5f;
			D *= E / (1.0f - v * v);



		Eigen::MatrixXf ke = B.transpose() * D * B * area;


		const Matrix3& Re = t.R;
		ReT = Matrix3::Transpose(Re);

		Vector3 displacements[3]
		{
			(ReT * m_PhyxelsPos[idxA]) - m_PhyxelsPosInitial[idxA],
			(ReT * m_PhyxelsPos[idxB]) - m_PhyxelsPosInitial[idxB],
			(ReT * m_PhyxelsPos[idxC]) - m_PhyxelsPosInitial[idxC],
		};

		Eigen::VectorXf Q(6);
		for (int i = 0; i < 3; ++i)
		{
			Q(i * 2 + 0) = displacements[i].x;
			Q(i * 2 + 1) = displacements[i].y;
		}

		Eigen::VectorXf force = ke * Q;


		for (int j = 0; j < 3; ++j)
		{
			int idxJ = t.phyxels[j];

			Vector3 B = Re * Vector3(force(j * 2 + 0), force(j * 2 + 1), 0.0f);
			m_Solver.m_B[idxJ] += B * dt;

			for (int k = 0; k < 3; ++k)
			{

				int idxK = t.phyxels[k];

				Matrix3 sub_ke = Matrix3::ZeroMatrix;
				sub_ke(0, 0) = ke(j * 2, k * 2);
				sub_ke(1, 0) = ke(j * 2 + 1, k * 2);
				sub_ke(0, 1) = ke(j * 2, k * 2 + 1);
				sub_ke(1, 1) = ke(j * 2 + 1, k * 2 + 1);


				Matrix3 transformed =  (Re * sub_ke * ReT) * dt2_viscos;		
				m_Solver.m_A(idxJ, idxK) -= transformed;

			}
		}

	}


}

void FEM3Noded::Render_DrawingToVisualiser()
{
	unsigned int i;
	for (i = 0; i < m_NumTriangles; ++i)
	{
	const FETriangleLite& t = m_Triangles[i];

	Vector3 a = m_PhyxelsPos[t.phyxels[0]];
	Vector3 b = m_PhyxelsPos[t.phyxels[1]];
	Vector3 c = m_PhyxelsPos[t.phyxels[2]];

	Vector3 centre = (a + b + c) * 0.3333333f;

	//	memcpy(&rotation._11, &t.R(0,0), 9 * sizeof(float));
	//NCLDebug::DrawMatrix(t.R, centre, 0.02f);

	/*NCLDebug::DrawThickLine(a, b, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	NCLDebug::DrawThickLine(b, c, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	NCLDebug::DrawThickLine(a, c, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));*/

	NCLDebug::DrawHairLine(a, b, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	NCLDebug::DrawHairLine(b, c, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	NCLDebug::DrawHairLine(a, c, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
	}
}
