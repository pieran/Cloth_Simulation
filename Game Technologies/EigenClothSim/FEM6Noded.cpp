#include "FEM6Noded.h"
#include <ncltech\NCLDebug.h>




FEM6Noded::FEM6Noded()
{
	m_TimeStep = max_timestep;
	m_StepCounter = 0;
	m_StepsBeforeIncrease = 10;
}

FEM6Noded::~FEM6Noded()
{
}

void FEM6Noded::simulation_OnClothDesignChanged(ClothDesignEntity<6>* design)
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
			m_PhyxelForces[i] = GRAVITY / float(m_NumTotal);
		}


		m_PhyxelTexCoords[i] = Vector2(v.texCoord.x(), v.texCoord.y());
	}



	//Populate Triangles & Add Additional Mid-Points
	const std::vector<TriangleG<6>>& triangles = design->Triangles();
	m_NumTriangles = triangles.size();
	m_TriangleStiffness.resize(m_NumTriangles);
	m_RenderIndices.resize(m_NumTriangles * 3);
	m_Triangles.resize(m_NumTriangles);

	float totalArea = 0.0f;
	for (unsigned int i = 0; i < m_NumTriangles; ++i)
	{
		unsigned int i3 = i * 3;
		memcpy(m_Triangles[i].phyxels, triangles[i].verts, 6 * sizeof(uint));
		memcpy(&m_RenderIndices[i3], triangles[i].verts, 3 * sizeof(uint));

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

void  FEM6Noded::UpdateConstraints()
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

void FEM6Noded::Simulation_StepSimulation(float dt)
{
	m_ProfilingTotalTime.BeginTiming();
	m_ProfilingBuildMatricies.ResetTotalMs();
	m_ProfilingConjugateGradient.ResetTotalMs();
	m_ProfilingExternalCollisions.ResetTotalMs();
	m_ProfilingRotationNormals.ResetTotalMs();
	m_Solver.ResetProfilingData();

	//Handle External Constraints
	/*{
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
	}*/


	//Compute Rotation Matrix for each triangle
	dt = min(dt, 1.0f / 60.0f);

	uint iterations = 0;
	for (float timestep_accum = 0.0f; timestep_accum < dt;)
	{
		//iterations++;
		//Handle External Constraints
		//if (iterations == 10)
		/*{
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
		}*/



		for (uint i = 0; i < m_NumTotal; ++i)
		{
			if (m_PhyxelIsStatic[i])
			{
				float spd = (i < (m_NumTotal >> 1)) ? 0.5f : 0.0f;

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


bool FEM6Noded::ValidateVelocityTimestep()
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

void FEM6Noded::CollideEllipsoid(const Ellipsoid& e)
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

void FEM6Noded::GenFEMTriangle(unsigned int idx, FETriangleGeneric<6>& tri)
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


	//Stiffness Matrix
	BuildStiffnessMatrix(idx, tri);
}

void TrigGaussRuleInfo(int rule, int point_idx, float* tcoor, float* weight)
{
	//{zeta, p = rule, i = point, g1, g2, info = { { Null,Null,Null },0 } },

	switch (rule)
	{
		case 1:
		{
			tcoor[0] = 1.0f / 3.0f;
			tcoor[1] = 1.0f / 3.0f;
			tcoor[2] = 1.0f / 3.0f;
			*weight = 1.0f;
		}
		break;

		case 3:
		{
			tcoor[0] = 1.0f / 6.0f;
			tcoor[1] = 1.0f / 6.0f;
			tcoor[2] = 1.0f / 6.0f;
			tcoor[point_idx] = 2.0f / 3.0f;
			*weight = 1.0f / 3.0f;
		}
		break;

		case -3:
		{
			tcoor[0] = 1.0f / 2.0f;
			tcoor[1] = 1.0f / 2.0f;
			tcoor[2] = 1.0f / 2.0f;
			tcoor[point_idx] = 0.0f;
			*weight = 1.0f / 3.0f;
		}
		break;

		case 6:
		{		
			if (point_idx < 3)
			{
				float g1 = (float)((8.0 - sqrt(10.0) + sqrt(38.0 - 44.0 * sqrt(0.4))) / 18.0);
				tcoor[0] = g1;
				tcoor[1] = g1;
				tcoor[2] = g1;
				tcoor[point_idx] = 1.0f - 2.0f * g1;
				*weight = (float)((620.0 + sqrt(213125.0 - 53320.0 * sqrt(10.0))) / 3720.0);
			}
			else
			{
				float g2 = (float)((8.0 - sqrt(10.0) - sqrt(38.0 - 44.0 * sqrt(0.4))) / 18.0);
				tcoor[0] = g2;
				tcoor[1] = g2;
				tcoor[2] = g2;
				tcoor[point_idx-3] = 1.0f - 2.0f * g2;
				*weight = (float)((620.0 - sqrt(213125.0 - 53320.0 * sqrt(10.0))) / 3720.0);
			}
		}
		break;

		case 7:
		{
			float g1 = (float)((6.0 - sqrt(15.0)) / 21.0); 
			float g2 = (float)((6.0 + sqrt(15.0)) / 21.0);

			if (point_idx < 3)
			{
				tcoor[0] = g1;
				tcoor[1] = g1;
				tcoor[2] = g1;
				tcoor[point_idx] = 1.0f - 2.0f * g1;
				*weight = (float)((155.0 - sqrt(15.0)) / 1200.0);
			}
			else if (point_idx == 6)
			{
				tcoor[0] = 1.0f / 3.0f;
				tcoor[1] = 1.0f / 3.0f;
				tcoor[2] = 1.0f / 3.0f;
				*weight = 9.0f / 40.0f;
			}
			else
			{
				tcoor[0] = g2;
				tcoor[1] = g2;
				tcoor[2] = g2;
				tcoor[point_idx - 3] = 1.0f - 2.0f * g2;
				*weight = (float)((155.0 + sqrt(15.0)) / 1200.0);
			}
		}
		break;

		default:
			assert(FALSE);
			break;
	}
}

void Trig6IsoPShapeFunDer(const Vector2* v, const float* tcoor, Vector2* dN, float* jDet)
{
	float dx4 = v[3].x - (v[0].x + v[1].x) / 2;
	float dx5 = v[4].x - (v[1].x + v[2].x) / 2;
	float dx6 = v[5].x - (v[2].x + v[0].x) / 2;
	float dy4 = v[3].y - (v[0].y + v[1].y) / 2;
	float dy5 = v[4].y - (v[1].y + v[2].y) / 2;
	float dy6 = v[5].y - (v[2].y + v[0].y) / 2;

	float Nf[6]{
		tcoor[0]*(2 * tcoor[0] - 1),
		tcoor[1]*(2 * tcoor[1] - 1),
		tcoor[2]*(2 * tcoor[2] - 1),
		4 * tcoor[0]*tcoor[1],
		4 * tcoor[1]*tcoor[2],
		4 * tcoor[2]*tcoor[0]
	};

	float Jx21 = v[1].x - v[0].x + 4 * (dx4*(tcoor[0] - tcoor[1]) + (dx5 - dx6)*tcoor[2]);
	float Jx32 = v[2].x - v[1].x + 4 * (dx5*(tcoor[1] - tcoor[2]) + (dx6 - dx4)*tcoor[0]);
	float Jx13 = v[0].x - v[2].x + 4 * (dx6*(tcoor[2] - tcoor[0]) + (dx4 - dx5)*tcoor[1]);
	float Jy12 = v[0].y - v[1].y + 4 * (dy4*(tcoor[1] - tcoor[0]) + (dy6 - dy5)*tcoor[2]);
	float Jy23 = v[1].y - v[2].y + 4 * (dy5*(tcoor[2] - tcoor[1]) + (dy4 - dy6)*tcoor[0]);
	float Jy31 = v[2].y - v[0].y + 4 * (dy6*(tcoor[0] - tcoor[2]) + (dy5 - dy4)*tcoor[1]);

	*jDet = Jx21*Jy31 - Jy12*Jx13;

	dN[0] = Vector2((4 * tcoor[0] - 1)*Jy23, (4 * tcoor[0] - 1)*Jx32);
	dN[1] = Vector2((4 * tcoor[1] - 1)*Jy31, (4 * tcoor[1] - 1)*Jx13);
	dN[2] = Vector2((4 * tcoor[2] - 1)*Jy12, (4 * tcoor[2] - 1)*Jx21);
	dN[3] = Vector2(4 * (tcoor[1]*Jy23 + tcoor[0]*Jy31), 4 * (tcoor[1]*Jx32 + tcoor[0]*Jx13));
	dN[4] = Vector2(4 * (tcoor[2]*Jy31 + tcoor[1]*Jy12), 4 * (tcoor[2]*Jx13 + tcoor[1]*Jx21));
	dN[5] = Vector2(4 * (tcoor[0]*Jy12 + tcoor[2]*Jy23), 4 * (tcoor[0]*Jx21 + tcoor[2]*Jx32));


	float invJDet = 1.0f / *jDet;
	for (int i = 0; i < 6; ++i)
	{
		dN[i].x *= invJDet;
		dN[i].y *= invJDet;
	}
}


void FEM6Noded::BuildStiffnessMatrix(unsigned int idx, FETriangleGeneric<6>& tri)
{
	Vector3 a = m_PhyxelsPosInitial[tri.phyxels[0]];
	Vector3 b = m_PhyxelsPosInitial[tri.phyxels[1]];
	Vector3 c = m_PhyxelsPosInitial[tri.phyxels[2]];

	const float E = 1000.0f;	//Youngs Modulus
	const float v = 0.3f;		//Poisson coefficient
	Eigen::Matrix3f D;
	D.setZero();
	D(0, 0) = 1.0f;
	D(1, 0) = v;
	D(0, 1) = v;
	D(1, 1) = 1.0f;
	D(2, 2) = (1.0f - v) * 0.5f;
	D *= E / (1.0f - v * v);

	const int gauss_formation = -3;
	//If[!MemberQ[{1, -3, 3, 6, 7}, p],
	//	Print["Illegal p"];

	float weight, jDet;
	float tcoor[3];
	Vector2 pos[6], dN[6];

	for (int i = 0; i < 6; ++i)
	{
		Vector3& pos3 = m_PhyxelsPosInitial[tri.phyxels[i]];
		pos[i] = Vector2(pos3.x, pos3.y);
	}

	Eigen::Matrix<float, 3, 12> Be;
	Be.setZero();

	m_TriangleStiffness[idx].setZero();
	for (int k = 0; k < abs(gauss_formation); ++k)
	{
		TrigGaussRuleInfo(gauss_formation, k, tcoor, &weight);
		Trig6IsoPShapeFunDer(pos, tcoor, dN, &jDet);

		for (int i = 0; i < 6; ++i)
		{
			Be(0, i * 2) = dN[i].x;
			Be(1, i * 2 + 1) = dN[i].y;

			Be(2, i * 2) = dN[i].y;
			Be(2, i * 2 + 1) = dN[i].x;
		}


		float c = weight * jDet * 0.5f;
		m_TriangleStiffness[idx] += c * Be.transpose() * D * Be;
	}
}

void FEM6Noded::BuildRotationAndNormals()
{

	memset(&m_PhyxelNormals[0], 0, m_NumTotal * sizeof(Vec3));


	//Build Rotations and Sum Up Normals
#pragma omp parallel for
	for (int i = 0; i < m_NumTriangles; ++i)
	{
		auto& t = m_Triangles[i];
		int idxA = t.phyxels[0];
		int idxB = t.phyxels[1];
		int idxC = t.phyxels[2];

		Vector3 a = m_PhyxelsPos[idxA];
		Vector3 b = m_PhyxelsPos[idxB];
		Vector3 c = m_PhyxelsPos[idxC];

		t.R = BuildRotationMatrix(a, b, c) * t.rotBase;


	}

	for (auto& tri : m_Triangles)
	{
		Vector3 a = m_PhyxelsPos[tri.phyxels[0]];
		Vector3 b = m_PhyxelsPos[tri.phyxels[1]];
		Vector3 c = m_PhyxelsPos[tri.phyxels[2]];

		Vector3 normal = Vector3::Cross(a - c, b - c);
		normal.Normalise();

		m_PhyxelNormals[tri.phyxels[0]] += normal;
		m_PhyxelNormals[tri.phyxels[1]] += normal;
		m_PhyxelNormals[tri.phyxels[2]] += normal;
		m_PhyxelNormals[tri.phyxels[3]] += normal;
		m_PhyxelNormals[tri.phyxels[4]] += normal;
		m_PhyxelNormals[tri.phyxels[5]] += normal;
	}

	//Normalise Vertex Normals
#pragma omp parallel for
	for (int i = 0; i < m_NumTotal; ++i)
	{
		m_PhyxelNormals[i].Normalise();
	}
}

Matrix3 FEM6Noded::BuildRotationMatrix(const Vector3& a, const Vector3& b, const Vector3& c)
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

void FEM6Noded::SimpleCorotatedBuildAMatrix(float dt)
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
		const FETriangleGeneric<6>& t = m_Triangles[i];
		uint idxA = t.phyxels[0];
		uint idxB = t.phyxels[1];
		uint idxC = t.phyxels[2];
		uint idxD = t.phyxels[3];
		uint idxE = t.phyxels[4];
		uint idxF = t.phyxels[5];

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





		const StiffnessMatrix& ke = m_TriangleStiffness[i];

		const Matrix3 Re = t.R;
		ReT = Matrix3::Transpose(Re);

		Vector3 displacements[6]
		{
			(ReT * m_PhyxelsPos[idxA]) - m_PhyxelsPosInitial[idxA],
			(ReT * m_PhyxelsPos[idxB]) - m_PhyxelsPosInitial[idxB],
			(ReT * m_PhyxelsPos[idxC]) - m_PhyxelsPosInitial[idxC],
			(ReT * m_PhyxelsPos[idxD]) - m_PhyxelsPosInitial[idxD],
			(ReT * m_PhyxelsPos[idxE]) - m_PhyxelsPosInitial[idxE],
			(ReT * m_PhyxelsPos[idxF]) - m_PhyxelsPosInitial[idxF],
		};

		Eigen::VectorXf Q(12);
		for (int i = 0; i < 6; ++i)
		{
			Q(i * 2 + 0) = displacements[i].x;
			Q(i * 2 + 1) = displacements[i].y;
		}

		Eigen::VectorXf force = ke * Q;

		m_Solver.SetMaxIterations(10);

		for (int j = 0; j < 6; ++j)
		{
			int idxJ = t.phyxels[j];

			Vector3 B = Re * Vector3(force(j * 2 + 0), force(j * 2 + 1), 0.0f);
			m_Solver.m_B[idxJ] -= B * dt;

			for (int k = 0; k < 6; ++k)
			{

				int idxK = t.phyxels[k];

				Matrix3 sub_ke = Matrix3::ZeroMatrix;
				sub_ke(0, 0) = ke(j * 2, k * 2);
				sub_ke(1, 0) = ke(j * 2 + 1, k * 2);
				sub_ke(0, 1) = ke(j * 2, k * 2 + 1);
				sub_ke(1, 1) = ke(j * 2 + 1, k * 2 + 1);


				Matrix3 transformed = (Re * sub_ke * ReT) * dt2_viscos;
				m_Solver.m_A(idxJ, idxK) += transformed;
			}
		}

	}


}

void FEM6Noded::Render_DrawingToVisualiser()
{
	unsigned int i;
	for (i = 0; i < m_NumTriangles; ++i)
	{
		const auto& t = m_Triangles[i];

		Vector3 a = m_PhyxelsPos[t.phyxels[0]];
		Vector3 b = m_PhyxelsPos[t.phyxels[1]];
		Vector3 c = m_PhyxelsPos[t.phyxels[2]];

		Vector3 ab = m_PhyxelsPos[t.phyxels[3]];
		Vector3 bc = m_PhyxelsPos[t.phyxels[4]];
		Vector3 ac = m_PhyxelsPos[t.phyxels[5]];


		//Vector3 centre = (a + b + c) * 0.3333333f;

		//	memcpy(&rotation._11, &t.R(0,0), 9 * sizeof(float));
		//NCLDebug::DrawMatrix(t.R, centre, 0.02f);

		/*NCLDebug::DrawThickLine(a, b, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
		NCLDebug::DrawThickLine(b, c, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
		NCLDebug::DrawThickLine(a, c, 0.002f, Vector4(0.0f, 0.0f, 0.0f, 1.0f));*/

		NCLDebug::DrawHairLine(a, ab, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
		NCLDebug::DrawHairLine(ab, b, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
		NCLDebug::DrawHairLine(b, bc, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
		NCLDebug::DrawHairLine(bc, c, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
		NCLDebug::DrawHairLine(c, ac, Vector4(0.0f, 0.0f, 0.0f, 1.0f));
		NCLDebug::DrawHairLine(ac, a, Vector4(0.0f, 0.0f, 0.0f, 1.0f));

		NCLDebug::DrawHairLine(ab, bc, Vector4(0.0f, 0.0f, 0.0f, 0.5f));
		NCLDebug::DrawHairLine(ab, ac, Vector4(0.0f, 0.0f, 0.0f, 0.5f));
		NCLDebug::DrawHairLine(ac, bc, Vector4(0.0f, 0.0f, 0.0f, 0.5f));
	}
}
