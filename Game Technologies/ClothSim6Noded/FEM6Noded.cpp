#include "FEM6Noded.h"
#include <ncltech\NCLDebug.h>




FEM6Noded::FEM6Noded()
{
	m_TimeStep = max_timestep;
	m_StepCounter = 0;
	m_StepsBeforeIncrease = 10;

	InitGaussWeights();
}

FEM6Noded::~FEM6Noded()
{
}

void FEM6Noded::simulation_OnClothDesignChanged(uint nTris, FETriangle* tris, uint nVerts, FEVertDescriptor* verts)
{


	//Gen Vertices
	m_NumTotal = nVerts;

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

		const FEVertDescriptor& v = verts[i];

		m_PhyxelsPos[i] = v.pos;
		m_PhyxelsPosInitial[i] = v.ipos;
		m_PhyxelIsStatic[i] = v.isStatic;

		if (!m_PhyxelIsStatic[i])
		{
			m_PhyxelForces[i] = GRAVITY / float(m_NumTotal);
		}


		m_PhyxelTexCoords[i] = v.tCoord;
	}



	//Populate Triangles & Add Additional Mid-Points
	m_NumTriangles = nTris;
	m_TriangleStiffness.resize(m_NumTriangles);
	m_RenderIndices.resize(m_NumTriangles * 3);
	m_Triangles.resize(m_NumTriangles);

	float totalArea = 0.0f;
	for (unsigned int i = 0; i < m_NumTriangles; ++i)
	{
		unsigned int i3 = i * 3;
		memcpy(m_Triangles[i].phyxels, tris[i].phyxels, 6 * sizeof(uint));
	
		totalArea += GenFEMTriangle(i, m_Triangles[i]);
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


	//Compute Rotation Matrix for each triangle
	dt = min(dt, 1.0f / 60.0f);

	uint iterations = 0;
	for (float timestep_accum = 0.0f; timestep_accum < dt;)
	{
		//iterations++;
		for (uint i = 0; i < m_NumTotal; ++i)
		{
			if (m_PhyxelIsStatic[i])
			{
				float spd = (i < (m_NumTotal >> 1)) ? 0.2f : 0.0f;

				m_PhyxelsVel[i].y = -cos(angle) * spd;// *(static_itr < 33 ? 1.0f : -1.0f);
				m_PhyxelsVel[i].z = sin(angle) * spd;// *(static_itr++ < 33 ? 1.0f : -1.0f);

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
		angle += m_TimeStep * 2.f;

		//Reset All Accumulative Data & Update Positions	
		memcpy(&m_PhyxelsPos[0].x, &m_PhyxelsPosTemp[0].x, m_NumTotal * sizeof(Vector3));
	}

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

float FEM6Noded::GenFEMTriangle(unsigned int idx, FETriangle& tri)
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
	//tri.Area = 0.5f * abs(e1.x * e2.y - e1.y * e2.x);

	//Rotation
	//tri.rotBase = Matrix3::Transpose(BuildRotationMatrix(a, b, c));


	//Stiffness Matrix
	BuildStiffnessMatrix(idx, tri);

	return 0.5f * abs(e1.x * e2.y - e1.y * e2.x); //return area of triangle
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
			tcoor[point_idx - 3] = 1.0f - 2.0f * g2;
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
		tcoor[0] * (2 * tcoor[0] - 1),
		tcoor[1] * (2 * tcoor[1] - 1),
		tcoor[2] * (2 * tcoor[2] - 1),
		4 * tcoor[0] * tcoor[1],
		4 * tcoor[1] * tcoor[2],
		4 * tcoor[2] * tcoor[0]
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
	dN[3] = Vector2(4 * (tcoor[1] * Jy23 + tcoor[0] * Jy31), 4 * (tcoor[1] * Jx32 + tcoor[0] * Jx13));
	dN[4] = Vector2(4 * (tcoor[2] * Jy31 + tcoor[1] * Jy12), 4 * (tcoor[2] * Jx13 + tcoor[1] * Jx21));
	dN[5] = Vector2(4 * (tcoor[0] * Jy12 + tcoor[2] * Jy23), 4 * (tcoor[0] * Jx21 + tcoor[2] * Jx32));


	float invJDet = 1.0f / *jDet;
	for (int i = 0; i < 6; ++i)
	{
		dN[i].x *= invJDet;
		dN[i].y *= invJDet;
	}
}


void FEM6Noded::BuildStiffnessMatrix(unsigned int idx, FETriangle& tri)
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
/*#pragma omp parallel for
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
	}*/

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



void FEM6Noded::InitGaussWeights()
{

	//Define Gauss points for numerical intergration
	GaussPoint12(0, 0) = 0.873821971016996;
	GaussPoint12(0, 1) = 0.063089014491502;
	GaussPoint12(0, 2) = 0.063089014491502;
	GaussPoint12(1, 0) = 0.063089014491502;
	GaussPoint12(1, 1) = 0.063089014491502;
	GaussPoint12(1, 2) = 0.873821971016996;
	GaussPoint12(2, 0) = 0.063089014491502;
	GaussPoint12(2, 1) = 0.873821971016996;
	GaussPoint12(2, 2) = 0.063089014491502;

	for (int i = 0; i < 3; ++i)
		GaussWeight12(i) = 0.050844906370207;

	GaussPoint12(3, 0) = 0.501426509658179;
	GaussPoint12(3, 1) = 0.249286745170910;
	GaussPoint12(3, 2) = 0.249286745170910;
	GaussPoint12(4, 0) = 0.249286745170910;
	GaussPoint12(4, 1) = 0.249286745170910;
	GaussPoint12(4, 2) = 0.501426509658179;
	GaussPoint12(5, 0) = 0.249286745170910;
	GaussPoint12(5, 1) = 0.501426509658179;
	GaussPoint12(5, 2) = 0.249286745170910;

	for (int i = 3; i < 6; ++i)
		GaussWeight12(i) = 0.116786275726379;


	GaussPoint12(6, 0) = 0.636502499121399;
	GaussPoint12(6, 1) = 0.310352451033784;
	GaussPoint12(6, 2) = 0.053145049844817;
	GaussPoint12(7, 0) = 0.636502499121399;
	GaussPoint12(7, 1) = 0.053145049844817;
	GaussPoint12(7, 2) = 0.310352451033784;
	GaussPoint12(8, 0) = 0.310352451033784;
	GaussPoint12(8, 1) = 0.636502499121399;
	GaussPoint12(8, 2) = 0.053145049844817;

	GaussPoint12(9, 0) = 0.310352451033784;
	GaussPoint12(9, 1) = 0.053145049844817;
	GaussPoint12(9, 2) = 0.636502499121399;
	GaussPoint12(10, 0) = 0.053145049844817;
	GaussPoint12(10, 1) = 0.310352451033784;
	GaussPoint12(10, 2) = 0.636502499121399;
	GaussPoint12(11, 0) = 0.053145049844817;
	GaussPoint12(11, 1) = 0.636502499121399;
	GaussPoint12(11, 2) = 0.310352451033784;

	for (int i = 6; i < 12; ++i)
		GaussWeight12(i) = 0.082851075618374;

}


void FEM6Noded::BuildTransformationMatrix(const FETriangle& tri, const std::vector<Vector3>& pos, const Vec3& gaussPoint, const float warp_angle, Mat33& out_t, JaMatrix& out_ja, Mat26& out_dn)
{
	//Build local fill direction (local warp direction)
	Vec3 vFill(-sin(warp_angle), cos(warp_angle), 0.0f);

	//Build natural coordinate basis vectors
	out_dn(0, 0) = 4 * gaussPoint.x() - 1.f;
	out_dn(0, 1) = 0.0f;
	out_dn(0, 2) = -4 * gaussPoint.z() + 1.f;
	out_dn(0, 3) = 4 * gaussPoint.y();
	out_dn(0, 4) = -4 * gaussPoint.y();
	out_dn(0, 5) = 4 * (gaussPoint.z() - gaussPoint.x());

	out_dn(1, 0) = 0.f;
	out_dn(1, 1) = 4 * gaussPoint.y() - 1.f;
	out_dn(1, 2) = -4 * gaussPoint.z() + 1.f;
	out_dn(1, 3) = 4 * gaussPoint.x();
	out_dn(1, 4) = 4 * (gaussPoint.z() - gaussPoint.y());
	out_dn(1, 5) = -4 * gaussPoint.x();


	Vec3 V_xi(0, 0, 0), V_eta(0, 0, 0);
	Vector3 tmp;
	for (int i = 0; i < 6; ++i)
	{
		tmp = pos[tri.phyxels[i]] * out_dn(0, i);
		V_xi += Vec3(tmp.x, tmp.y, tmp.z);

		tmp = pos[tri.phyxels[i]] * out_dn(1, i);
		V_eta += Vec3(tmp.x, tmp.y, tmp.z);
	}

	Vec3 V_z = V_xi.cross(V_eta);
	V_z.normalize();

	if (V_z.z() < 0.0f)
	{
		V_z = V_eta.cross(V_xi);
		V_z.normalize();
	}

	Vec3 V_x = vFill.cross(V_z);
	V_x.normalize();

	Vec3 V_y = V_z.cross(V_x);
	V_y.normalize();

	//Build the transformation matrix
	out_t(0, 0) = V_x(0); out_t(0, 1) = V_x(1); out_t(0, 2) = V_x(2);
	out_t(1, 0) = V_y(0); out_t(1, 1) = V_y(1); out_t(1, 2) = V_y(2);
	out_t(2, 0) = V_z(0); out_t(2, 1) = V_z(1); out_t(2, 2) = V_z(2);

	//Build the Jacobian Matrix
	out_ja(0, 0) = V_xi.dot(V_x);
	out_ja(0, 1) = V_xi.dot(V_y);
	out_ja(1, 0) = V_eta.dot(V_x);
	out_ja(1, 1) = V_eta.dot(V_y);
}

void FEM6Noded::CalcBMatrix(const FETriangle& tri, const std::vector<Vector3>& pos, const Vec3& gaussPoint, const float warp_angle, const VDisplacement& displacements, BMatrix& out_b, JaMatrix& out_ja, GMatrix& out_g)
{
	Mat26 DN;
	Mat33 T;

	BuildTransformationMatrix(tri, pos, gaussPoint, warp_angle, T, out_ja, DN);

	Mat22 Ja_inv = out_ja.inverse();

	auto a = Ja_inv * DN;

	BMatrix B_l, B_nl;

	for (int i = 0; i < 6; ++i)
	{
		B_l(0, i * 3) = T(0, 0) * a(0, i);
		B_l(0, i * 3 + 1) = T(0, 1) * a(0, i);
		B_l(0, i * 3 + 2) = T(0, 2) * a(0, i);

		B_l(1, i * 3) = T(1, 0) * a(1, i);
		B_l(1, i * 3 + 1) = T(1, 1) * a(1, i);
		B_l(1, i * 3 + 2) = T(1, 2) * a(1, i);

		B_l(2, i * 3) = T(0, 0) * a(1, i) + T(1, 0) * a(0, i);
		B_l(2, i * 3 + 1) = T(0, 1) * a(1, i) + T(1, 1) * a(0, i);
		B_l(2, i * 3 + 2) = T(0, 2) * a(1, i) + T(1, 2) * a(0, i);
	}

	//Build G Matrix
	for (int n = 0; n < 6; ++n)
	{
		for (int i = 0; i < 3; ++i)
		{
			out_g(i, n * 3) = T(i, 0) * a(0, n);
			out_g(i, n * 3 + 1) = T(i, 1) * a(0, n);
			out_g(i, n * 3 + 2) = T(i, 2) * a(0, n);
		}

		for (int i = 3; i < 6; ++i)
		{
			out_g(i, n * 3) = T(i-3, 0) * a(1, n);
			out_g(i, n * 3 + 1) = T(i-3, 1) * a(1, n);
			out_g(i, n * 3 + 2) = T(i-3, 2) * a(1, n);
		}
	}

	//Find displacement dependant terms of local Bmatrix
	auto delta = out_g * displacements;

	Eigen::Matrix<float, 3, 6> derivatives;
	derivatives.setZero();
	derivatives(0, 0) = delta[0];
	derivatives(0, 1) = delta[1];
	derivatives(0, 2) = delta[2];
	derivatives(1, 3) = delta[3];
	derivatives(1, 4) = delta[4];
	derivatives(1, 5) = delta[5];

	derivatives(2, 0) = delta[3];
	derivatives(2, 1) = delta[4];
	derivatives(2, 2) = delta[5];
	derivatives(2, 3) = delta[0];
	derivatives(2, 4) = delta[1];
	derivatives(2, 5) = delta[2];

	B_nl = 0.5 * derivatives * out_g;

	out_b = B_l + B_nl;
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


	Mat33 E;
	const float Y = 2500.0f;	//Youngs Modulus
	const float v = 0.3f;		//Poisson coefficient
	E.setZero();
	E(0, 0) = 1.0f;
	E(1, 0) = v;
	E(0, 1) = v;
	E(1, 1) = 1.0f;
	E(2, 2) = (1.0f - v) * 0.5f;
	E *= Y / (1.0f - v * v);

	BMatrix B_nl, B_0;
	JaMatrix Ja;
	GMatrix G;

	VDisplacement d_g, d_0, force;
	d_0.setZero();

	Eigen::Matrix<float, 18, 18> K_E, K_S, K_T;
	for (uint i = 0; i < m_NumTriangles; ++i)
	{
		FETriangle& tri = m_Triangles[i];

		K_E.setZero();
		K_S.setZero();

		for (int j = 0; j < 6; ++j)
		{
			d_g(j * 3 + 0, 0) = m_PhyxelsPos[tri.phyxels[j]].x - m_PhyxelsPosInitial[tri.phyxels[j]].x;
			d_g(j * 3 + 1, 0) = m_PhyxelsPos[tri.phyxels[j]].y - m_PhyxelsPosInitial[tri.phyxels[j]].y;
			d_g(j * 3 + 2, 0) = m_PhyxelsPos[tri.phyxels[j]].z - m_PhyxelsPosInitial[tri.phyxels[j]].z;
		}

		 force.setZero();
		for (uint j = 0; j < 12; ++j)
		{
			CalcBMatrix(tri, m_PhyxelsPosInitial, GaussPoint12.row(j), 0.0f, d_g, B_nl, Ja, G);

			Vec3 strain = B_nl * d_g;
			Vec3 stress = E * strain;

			CalcBMatrix(tri, m_PhyxelsPos, GaussPoint12.row(j), 0.0f, d_0, B_0, Ja, G);
		
			float area = Ja.determinant() * 0.5f;
			if (area < 0)
			{
				printf("ERROR:: Element %d:%d has a negative area!!!\n", i, j);
			}
			
			float tfactor = GaussWeight12(j) * area;


			Eigen::Matrix<float, 6, 6> M; M.setZero();
			M(0, 0) = stress.x();
			M(1, 1) = stress.x();
			M(2, 2) = stress.x();

			M(0, 3) = stress.z();
			M(1, 4) = stress.z();
			M(2, 5) = stress.z();
			M(3, 0) = stress.z();
			M(4, 1) = stress.z();
			M(5, 2) = stress.z();

			M(3, 3) = stress.y();
			M(4, 4) = stress.y();
			M(5, 5) = stress.y();



			K_E += B_0.transpose() * E * B_nl * tfactor;
			K_S += G.transpose() * M * G * tfactor;

			force += B_0.transpose() * stress * tfactor;
		}

		K_T = K_E + K_S;


		//Convert Eigen Matrices back to global A Matrix + B Vectors for global solver
		for (uint j = 0; j < 6; ++j)
		{
			m_Solver.m_B[tri.phyxels[j]] -= Vector3(force(j * 3), force(j * 3 + 1), force(j * 3 + 2)) * dt;

		//	m_PhyxelsPos[tri.phyxels[j]] += Vector3(force(j * 3), force(j * 3 + 1), force(j * 3 + 2)) * dt * 0.001f;

			for (uint k = 0; k < 6; ++k)
			{
				uint idxJ = tri.phyxels[j];
				uint idxK = tri.phyxels[k];

				Matrix3 submtx;
				submtx._11 = K_T(j * 3 + 0, k * 3 + 0);
				submtx._12 = K_T(j * 3 + 0, k * 3 + 1);
				submtx._13 = K_T(j * 3 + 0, k * 3 + 2);
				submtx._21 = K_T(j * 3 + 1, k * 3 + 0);
				submtx._22 = K_T(j * 3 + 1, k * 3 + 1);
				submtx._23 = K_T(j * 3 + 1, k * 3 + 2);
				submtx._31 = K_T(j * 3 + 2, k * 3 + 0);
				submtx._32 = K_T(j * 3 + 2, k * 3 + 1);
				submtx._33 = K_T(j * 3 + 2, k * 3 + 2);

				m_Solver.m_A(idxJ, idxK) += submtx * dt * dt;
			}
		}


	}


#if 0
	//for (int i = 0; i < (int)m_NumTriangles; ++i)
	//{
	//	const FETriangle& t = m_Triangles[i];
	//	uint idxA = t.phyxels[0];
	//	uint idxB = t.phyxels[1];
	//	uint idxC = t.phyxels[2];
	//	uint idxD = t.phyxels[3];
	//	uint idxE = t.phyxels[4];
	//	uint idxF = t.phyxels[5];

	//	Vector3 triNormals[3] = { m_PhyxelNormals[idxA], m_PhyxelNormals[idxB], m_PhyxelNormals[idxC] };

	//	Vector3 face_normal = (triNormals[0] + triNormals[1] + triNormals[2]) / 3.0f;
	//	face_normal.Normalise();

	//	/*//Air Drag
	//	{
	//		const float Cd = 0.55f; //Drag Coefficient
	//		const float Cl = 0.55f; //Lift Coefficient
	//		const float p = 1.293f; //Air Density

	//								//Vrel = Vi - u (Where u is velocity field of the wind)
	//								//F = 0.5f * Cd * p * |Vrel|^2 * A * (n . Vrel) . (-Vrel)
	//		Vector3 vrel = (m_PhyxelsVel[idxA] + m_PhyxelsVel[idxB] + m_PhyxelsVel[idxC]) / 3.0f - WIND;
	//		Vector3 vrel_norm = vrel; vrel_norm.Normalise();

	//		if (Vector3::Dot(face_normal, vrel) < 0) face_normal = -face_normal;

	//		Vector3 lift_dir = Vector3::Cross(Vector3::Cross(face_normal, vrel_norm), vrel_norm);

	//		float vrelsq = Vector3::Dot(vrel, vrel);
	//		float drag_coef = 0.5f * Cd * p * vrelsq * t.Area * Vector3::Dot(face_normal, vrel_norm);
	//		float lift_coef = 0.5f * Cl * p * vrelsq * t.Area * Vector3::Dot(lift_dir, vrel_norm);

	//		Vector3 Fid = (-vrel_norm * drag_coef) / 3.0f;		//Drag
	//		Vector3 Fil = (lift_dir * lift_coef) / 3.0f;		//Lift

	//		Vector3 totalForce = (Fid + Fil);
	//		m_Solver.m_B[idxA] += totalForce * dt;
	//		m_Solver.m_B[idxB] += totalForce * dt;
	//		m_Solver.m_B[idxC] += totalForce * dt;
	//	}*/





	//	/*const StiffnessMatrix& ke = m_TriangleStiffness[i];

	//	const Matrix3 Re = t.R;
	//	ReT = Matrix3::Transpose(Re);

	//	Vector3 displacements[6]
	//	{
	//		(ReT * m_PhyxelsPos[idxA]) - m_PhyxelsPosInitial[idxA],
	//		(ReT * m_PhyxelsPos[idxB]) - m_PhyxelsPosInitial[idxB],
	//		(ReT * m_PhyxelsPos[idxC]) - m_PhyxelsPosInitial[idxC],
	//		(ReT * m_PhyxelsPos[idxD]) - m_PhyxelsPosInitial[idxD],
	//		(ReT * m_PhyxelsPos[idxE]) - m_PhyxelsPosInitial[idxE],
	//		(ReT * m_PhyxelsPos[idxF]) - m_PhyxelsPosInitial[idxF],
	//	};

	//	Eigen::VectorXf Q(12);
	//	for (int i = 0; i < 6; ++i)
	//	{
	//		Q(i * 2 + 0) = displacements[i].x;
	//		Q(i * 2 + 1) = displacements[i].y;
	//	}

	//	Eigen::VectorXf force = ke * Q;

	//	m_Solver.SetMaxIterations(10);

	//	for (int j = 0; j < 6; ++j)
	//	{
	//		int idxJ = t.phyxels[j];

	//		Vector3 B = Re * Vector3(force(j * 2 + 0), force(j * 2 + 1), 0.0f);
	//		m_Solver.m_B[idxJ] -= B * dt;

	//		for (int k = 0; k < 6; ++k)
	//		{

	//			int idxK = t.phyxels[k];

	//			Matrix3 sub_ke = Matrix3::ZeroMatrix;
	//			sub_ke(0, 0) = ke(j * 2, k * 2);
	//			sub_ke(1, 0) = ke(j * 2 + 1, k * 2);
	//			sub_ke(0, 1) = ke(j * 2, k * 2 + 1);
	//			sub_ke(1, 1) = ke(j * 2 + 1, k * 2 + 1);


	//			Matrix3 transformed = (Re * sub_ke * ReT) * dt2_viscos;
	//			m_Solver.m_A(idxJ, idxK) += transformed;
	//		}
	//	}*/

	//}

#endif
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

		NCLDebug::DrawThickLine(a, ab, 0.02f, Vector4(0.0f, 0.0f, 0.0f, 0.3f));
		NCLDebug::DrawThickLine(ab, b, 0.02f, Vector4(0.0f, 0.0f, 0.0f, 0.3f));
		NCLDebug::DrawThickLine(b, bc, 0.02f, Vector4(0.0f, 0.0f, 0.0f, 0.3f));
		NCLDebug::DrawThickLine(bc, c, 0.02f, Vector4(0.0f, 0.0f, 0.0f, 0.3f));
		NCLDebug::DrawThickLine(c, ac, 0.02f, Vector4(0.0f, 0.0f, 0.0f, 0.3f));
		NCLDebug::DrawThickLine(ac, a, 0.02f, Vector4(0.0f, 0.0f, 0.0f, 0.3f));

		//NCLDebug::DrawHairLine(ab, bc, Vector4(0.0f, 0.0f, 0.0f, 0.5f));
		//NCLDebug::DrawHairLine(ab, ac, Vector4(0.0f, 0.0f, 0.0f, 0.5f));
		//NCLDebug::DrawHairLine(ac, bc, Vector4(0.0f, 0.0f, 0.0f, 0.5f));







		Vector3 sp = Vector3(1.f, 1.f, 1.f) / 3.f;
		Vector3 centre = (a + b + c + ab + ac + bc) / 6.f;

		Vector3 vxi, veta;

		vxi = a * 4.f * sp.x - 1.f;
		vxi += b * 0.0f;
		vxi += c * (-4.f * sp.z + 1.f);
		vxi += ab * 4.f * sp.y;
		vxi += bc * (-4.f * sp.y);
		vxi += ac * 4.f * (sp.z - sp.x);

		veta = a * 0.0f;
		veta += b * 4.f * sp.y - 1.f;
		veta += c * (-4.f * sp.z + 1.f);
		veta += ab * 4.f * sp.x;
		veta += bc * 4.f * (sp.z - sp.y);
		veta += ac * (-4.f * sp.x);

		veta.Normalise();
		vxi.Normalise();

		int nsteps = 16;
		float step = 1.f / float(nsteps);

		auto get_point = [&](int ix, int iy)
		{
			float x = ix * step;
			float y = iy * step;

			Vector3 gp = Vector3(x, y, 0.f);
			gp.z = 1.f - (x + y);


			Vector3 p;
			p = a * (gp.x * (2.f * gp.x - 1.f));
			p += b * (gp.y * (2.f * gp.y - 1.f));
			p += c * (gp.z * (2.f * gp.z - 1.f));
			p += ab * (4.f * gp.x * gp.y);
			p += bc * (4.f * gp.y * gp.z);
			p += ac * (4.f * gp.x * gp.z);

			return p;
		};


		//out_dn(0, 0) = 4 * gaussPoint.x() - 1.f;
		//out_dn(0, 1) = 0.0f;
		//out_dn(0, 2) = -4 * gaussPoint.z() + 1.f;
		//out_dn(0, 3) = 4 * gaussPoint.y();
		//out_dn(0, 4) = -4 * gaussPoint.y();
		//out_dn(0, 5) = 4 * (gaussPoint.z() - gaussPoint.x());

		//out_dn(1, 0) = 0.f;
		//out_dn(1, 1) = 4 * gaussPoint.y() - 1.f;
		//out_dn(1, 2) = -4 * gaussPoint.z() + 1.f;
		//out_dn(1, 3) = 4 * gaussPoint.x();
		//out_dn(1, 4) = 4 * (gaussPoint.z() - gaussPoint.y());
		//out_dn(1, 5) = -4 * gaussPoint.x();


		//Vec3 V_xi(0, 0, 0), V_eta(0, 0, 0);
		//Vector3 tmp;
		//for (int i = 0; i < 6; ++i)
		//{
		//	tmp = pos[tri.phyxels[i]] * out_dn(0, i);
		//	V_xi += Vec3(tmp.x, tmp.y, tmp.z);

		//	tmp = pos[tri.phyxels[i]] * out_dn(1, i);
		//	V_eta += Vec3(tmp.x, tmp.y, tmp.z);
		//}

		Vector3 inv_light_dir = Vector3(0.5f, 1.0f, -0.8f);
		inv_light_dir.Normalise();

		const Vector4 color = Vector4(0.f, 0.f, 1.f, 1.f);
		auto draw_line = [&](const Vector3& v1, const Vector3& v2)
		{
#if 1
			NCLDebug::DrawHairLine(v1, v2, color);
#else
			Vector3 norm = v2 - v1;
			norm.Normalise();

			float d = norm.Dot(inv_light_dir);
			d = min(max(0, d), 1);
			//d *= 0.5f;
			//NCLDebug::DrawHairLine(v1, v2, Vector4(color.x * d, color.y * d, color.z * d, color.w));

			NCLDebug::DrawHairLine(v1, v2, Vector4(color.x + d, color.y + d, color.z + d, color.w));
#endif
		};

	
		Vector3 xOldP;
		for (int ix = 0; ix <= nsteps; ++ix)
		{

			Vector3 oldp;
			for (int iy = 0; iy <= (nsteps - ix); ++iy)
			{
				Vector3 p = get_point(ix, iy);

				//NCLDebug::DrawPoint(p, 0.02f, Vector4(dir.x, dir.y, dir.z, 1.f));
				//NCLDebug::DrawHairLine(p, p + dir * 0.01f , Vector4(dir.x, dir.y, dir.z, 1.f));
				if (iy == nsteps - ix)
				{
					if (ix > 0)
						draw_line(xOldP, p);

					xOldP = p;
				}

				if (ix > 0)
				{
					Vector3 p3 = get_point(ix - 1, iy);

					draw_line(p, p3);
					if (iy > 0)draw_line(oldp, p3);
					
				}
				
				if (iy > 0)
				{
					draw_line(p, oldp);
				}

				oldp = p;
			}

			
		}



	}
}
