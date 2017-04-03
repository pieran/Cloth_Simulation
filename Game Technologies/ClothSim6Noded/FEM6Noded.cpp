#include "FEM6Noded.h"
#include <ncltech\NCLDebug.h>
#include "utils.h"

FEM6Noded::FEM6Noded()
{
	m_TimeStep = max_timestep;
	m_StepCounter = 0;
	m_StepsBeforeIncrease = 10;
	m_RenderMode = 0;

	InitGaussWeights();

	const float Y = 2500.0f;	//Youngs Modulus
	const float v = 0.3f;		//Poisson coefficient
	E.setZero();
	E(0, 0) = 1.0f;
	E(1, 0) = v;
	E(0, 1) = v;
	E(1, 1) = 1.0f;
	E(2, 2) = (1.0f - v) * 0.5f;
	E *= Y / (1.0f - v * v);
}

FEM6Noded::~FEM6Noded()
{
}

void FEM6Noded::simulation_OnClothDesignChanged(uint nTris, FETriangle* tris, uint nVerts, FEVertDescriptor* verts, uint nTangents, Vector3* tangents)
{


	//Gen Vertices
	m_NumTotal = nVerts;
	m_NumTangents = nTangents;

	m_Solver.AllocateMemory(m_NumTotal + m_NumTangents);
	m_Solver.m_A.resize(m_NumTotal + m_NumTangents);

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


	m_PhyxelsTangent.resize(m_NumTangents);
	m_PhyxelsTangentInitial.resize(m_NumTangents);
	memcpy(&m_PhyxelsTangent[0], tangents, m_NumTangents * sizeof(Vector3));
	memcpy(&m_PhyxelsTangentInitial[0], tangents, m_NumTangents * sizeof(Vector3));

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
			m_PhyxelForces[i] = Vector3(0.f, 0.f, 0.f);// GRAVITY / float(m_NumTotal);
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
		memcpy(m_Triangles[i].phyxels, tris[i].phyxels, sizeof(FETriangle));
	
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
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTangents; ++i)
	{
		m_Solver.m_Constraints[m_NumTotal + i] = Matrix3::Identity;
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

				//m_PhyxelsVel[i].y = -cos(angle) * spd;// *(static_itr < 33 ? 1.0f : -1.0f);
				//m_PhyxelsVel[i].z = sin(angle) * spd;// *(static_itr++ < 33 ? 1.0f : -1.0f);

				m_Solver.m_X[i] = m_PhyxelsVel[i];
			}
		}
		for (uint i = 0; i < m_NumTangents; ++i)
		{
			m_Solver.m_X[m_NumTotal + i] = Vector3(0.f, 0.f, 0.f);
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

			for (int i = 0; i < (int)m_NumTangents; ++i)
			{
				m_PhyxelsTangent[i] = m_PhyxelsTangent[i] + m_Solver.m_X[m_NumTotal + i];
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


void FEM6Noded::CalcBMatrix(const FETriangle& tri, const std::vector<Vector3>& pos, const std::vector<Vector3>& tans, const Vec3& gaussPoint, const float warp_angle, const Eigen::Matrix<float, 45, 1>& displacements, BMatrix_C1& out_b, JaMatrix& out_ja, Eigen::Matrix<float, 6, 45>& out_g)
{
	BaseMatrix DN, a;
	Mat33 T;
	BMatrix_C1 B_l, B_nl;
	Eigen::Matrix<float, 3, 6> derivatives;
	Eigen::Matrix<float, 6, 1> delta;

	//BuildTransformationMatrix(tri, pos, gaussPoint, warp_angle, T, out_ja, DN);
	CalcRotationC1(tri, pos, tans, Vector3(gaussPoint.x(), gaussPoint.y(), gaussPoint.z()), warp_angle, T, out_ja, DN);

	Mat22 Ja_inv = out_ja.inverse();

	a = Ja_inv * DN;



	for (int i = 0; i < 15; ++i)
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
	for (int n = 0; n < 15; ++n)
	{
		for (int i = 0; i < 3; ++i)
		{
			out_g(i, n * 3) = T(i, 0) * a(0, n);
			out_g(i, n * 3 + 1) = T(i, 1) * a(0, n);
			out_g(i, n * 3 + 2) = T(i, 2) * a(0, n);
		}

		for (int i = 3; i < 6; ++i)
		{
			out_g(i, n * 3) = T(i - 3, 0) * a(1, n);
			out_g(i, n * 3 + 1) = T(i - 3, 1) * a(1, n);
			out_g(i, n * 3 + 2) = T(i - 3, 2) * a(1, n);
		}
	}

	//Find displacement dependant terms of local Bmatrix
	delta = out_g * displacements;

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

	out_b = B_l +B_nl;
}

void FEM6Noded::SimpleCorotatedBuildAMatrix(float dt)
{
	const float dt2 = dt * dt;
	const float dt2_viscos = dt2 + V_SCALAR * dt;


	//#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTotal; ++i)
	{
		m_Solver.m_A(i, i) = Matrix3::Identity * m_PhyxelsMass[i] * 1.001f;
		m_Solver.m_B[i] = m_PhyxelsVel[i] * m_PhyxelsMass[i] +m_PhyxelForces[i] * dt;
	}

	for (int i = 0; i < (int)m_NumTangents; ++i)
	{
		m_Solver.m_A(m_NumTotal + i, m_NumTotal + i) = Matrix3::Identity;
		m_Solver.m_B[m_NumTotal + i] = Vector3(0.f, 0.f, 0.f);
	}

	Matrix3 RK, RKR, ReT;
	Vector3 posAB, posIAB;




	BMatrix_C1 B_nl, B_0;
	JaMatrix Ja;
	Eigen::Matrix<float, 6, 45> G;

	Eigen::Matrix<float, 45, 1> d_g, d_0, force;
	d_0.setZero();

	Eigen::Matrix<float, 45, 45> K_E, K_S, K_T;
	Eigen::Matrix<float, 6, 6> M; M.setZero();
	
	for (uint i = 0; i < m_NumTriangles; ++i)
	{
		FETriangle& tri = m_Triangles[i];

		K_E.setZero();
		K_S.setZero();

		for (int j = 0; j < 6; ++j)
		{
			Vector3 pos = (m_PhyxelsPos[tri.phyxels[j]] - m_PhyxelsPosInitial[tri.phyxels[j]]);

			d_g(j * 3 + 0, 0) = pos.x; d_g(j * 3 + 1, 0) = pos.y; d_g(j * 3 + 2, 0) = pos.z;
		}

		for (int j = 0; j < 9; ++j)
		{
			Vector3 tan = (m_PhyxelsTangent[tri.tangents[j]] - m_PhyxelsTangentInitial[tri.tangents[j]]);// *tri.tan_multipliers[j];

			d_g(j * 3 + 18, 0) = tan.x; d_g(j * 3 + 19, 0) = tan.y; d_g(j * 3 + 20, 0) = tan.z;
		}

		 force.setZero();
		for (uint j = 0; j < 12; ++j)
		{
			CalcBMatrix(tri, m_PhyxelsPosInitial, m_PhyxelsTangentInitial, GaussPoint12.row(j), 0.0f, d_g, B_nl, Ja, G);

			Vec3 strain = B_nl * d_g;
			Vec3 stress = E * strain;

			CalcBMatrix(tri, m_PhyxelsPos, m_PhyxelsTangent, GaussPoint12.row(j), 0.0f, d_0, B_0, Ja, G);
		
			float area = Ja.determinant() * 0.5f;
			if (area < 0)
			{
				printf("ERROR:: Element %d:%d has a negative area!!!\n", i, j);
			}
			
			float tfactor = GaussWeight12(j) * area;


			
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
		for (uint j = 0; j < 15; ++j)
		{
			uint idxJ = (j < 6) ? tri.phyxels[j] : m_NumTotal + tri.tangents[j-6];
			if (j < 6)
				m_Solver.m_B[idxJ] -= Vector3(force(j * 3), force(j * 3 + 1), force(j * 3 + 2)) * dt;
			if (j >= 6)
				m_Solver.m_B[idxJ] -= Vector3(force(j * 3), force(j * 3 + 1), force(j * 3 + 2)) * dt;// *tri.tan_multipliers[j - 6];

		//	m_PhyxelsPos[tri.phyxels[j]] += Vector3(force(j * 3), force(j * 3 + 1), force(j * 3 + 2)) * dt * 0.001f;

			for (uint k = 0; k < 15; ++k)
			{
				uint idxK = (k < 6) ? tri.phyxels[k] : m_NumTotal + tri.tangents[k-6];
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
}


void FEM6Noded::BuildNaturalCoordinateBasisVectors(const FETriangle& tri, const Vector3& gaussPoint, BaseMatrix& out_dn)
{
	Vector3 gp = gaussPoint;
	Vector3 gp2 = gp * gp;
	Vector3 gp3 = gp * gp * gp;
	Vector3 gp4 = gp2 * gp2;
	Vector3 gp5 = gp3 * gp2;

	float gxy = gp.x * gp.y;
	float gx2y = gp2.x * gp.y;
	float gxy2 = gp.x * gp2.y;
	float gx3y = gp3.x * gp.y;
	float gx2y2 = gp2.x * gp2.y;
	float gxy3 = gp.x * gp3.y;

	//Build natural coordinate basis vectors
	float coeffs_x[]
	{
		-10 * gp.x + 42 * gp2.x - 32 * gp3.x + 6 * (2 * gxy - 3 * gx2y - 2 * gxy2),
		6 * (gp2.y - gp3.y - 2 * gxy2),
		0,
		-16 * gp.y + 48 * (2 * gxy + gp2.y) - 32 * (3 * gx2y + gp3.y) - 96 * gxy2,
		0,
		-16 * gp.z + 48 * (2 * gp.x*gp.z + gp2.z) - 32 * (gp3.z + 3 * gp.z*gp2.x) - 96 * gp.x*gp2.z,

		0.5 * (-gp.y + 2 * gxy + 3 * gp2.y) + 3 * gx2y - 4 * gxy2 - gp3.y,
		0.5 * (-gp.z + 2 * gp.x*gp.z + 3 * gp2.z) + 3 * gp2.x*gp.z - 4 * gp.x*gp2.z - gp3.z,
		0.5 * (-gp.y + 6 * gxy + gp2.y) - 3 * gx2y - 4 * gxy2 + gp3.y,
		0,
		0.5 * (gp.z - 6 * gp.x*gp.z - gp2.z) + 3 * gp2.x*gp.z - gp3.z + 4 * gp.x*gp2.z,
		0,

		-4 * gp.y + 12 * (2 * gxy + gp2.y) - 8 * (3 * gx2y + 4 * gxy2 + gp3.y),
		0,
		-4 * gp.z + 12 * (2 * gp.x * gp.z + gp2.z) - 8 * (3 * gp2.x * gp.z + 4 * gp.x * gp2.z + gp3.z)
	};
	float coeffs_y[]
	{
		6 * (gp2.x - gp3.x - 2 * gx2y),
		-10 * gp.y + 42 * gp2.y - 32 * gp3.y + 6 * (2 * gxy - 3 * gxy2 - 2 * gx2y),
		6 * (gp2.z - gp3.z - 2 * gp.y*gp2.z),
		-16 * gp.x + 48 * (gp2.x + 2 * gxy) - 32 * (gp3.x + 3 * gxy2) - 96 * gx2y,
		-16 * gp.z + 48 * (2 * gp.y*gp.z + gp2.z) - 32 * (gp3.z + 3 * gp.z*gp2.y) - 96 * gp.y*gp2.z,
		0,

		0.5 * (-gp.x + gp2.x + 6 * gxy) + gp3.x - 4 * gx2y - 3 * gxy2,
		0,
		0.5 * (-gp.x + 3 * gp2.x + 2 * gxy) - gp3.x - 4 * gx2y + 3 * gxy2,
		0.5 * (-gp.z + 2 * gp.y*gp.z + 3 * gp2.z) + 3 * gp2.y*gp.z - 4 * gp.y*gp2.z - gp3.z,
		0,
		0.5 * (gp.z - 6 * gp.y*gp.z - gp2.z) + 3 * gp2.y*gp.z - gp3.z + 4 * gp.y*gp2.z,

		-4 * gp.x + 12 * (gp2.x + 2 * gxy) - 8 * (gp3.x + 4 * gx2y + 3 * gxy2),
		-4 * gp.z + 12 * (2 * gp.y * gp.z + gp2.z) - 8 * (3 * gp2.y * gp.z + 4 * gp.y * gp2.z + gp3.z),
		0
	};
	float coeffs_z[]
	{
		0,
		0,
		-10 * gp.z + 42 * gp2.z - 32 * gp3.z + 6 * (2 * gp.y*gp.z - 3 * gp.y*gp2.z - 2 * gp2.y*gp.z),
		0,
		-16 * gp.y + 48 * (gp2.y + 2 * gp.z * gp.y) - 32 * (3 * gp2.z*gp.y + gp3.y) - 96 * gp2.y*gp.z,
		-16 * gp.x + 48 * (gp2.x + 2 * gp.z * gp.x) - 32 * (3 * gp2.z*gp.x + gp3.x) - 96 * gp2.x*gp.z,

		0,
		0.5 * (-gp.x + gp2.x + 6 * gp.x*gp.z) + gp3.x - 4 * gp2.x*gp.z - 3 * gp.x*gp2.z,
		0,
		0.5 * (-gp.y + gp2.y + 6 * gp.y*gp.z) + gp3.y - 4 * gp2.y*gp.z - 3 * gp.y*gp2.z,
		0.5 * (gp.x - 3 * gp2.x - 2 * gp.x*gp.z) + gp3.x - 3 * gp.x*gp2.z + 4 * gp2.x*gp.z,
		0.5 * (gp.y - 3 * gp2.y - 2 * gp.y*gp.z) + gp3.y - 3 * gp.y*gp2.z + 4 * gp2.y*gp.z,

		0,
		-4 * gp.y + 12 * (gp2.y + 2 * gp.y * gp.z) - 8 * (gp3.y + 4 * gp2.y * gp.z + 3 * gp.y * gp2.z),
		-4 * gp.x + 12 * (gp2.x + 2 * gp.x * gp.z) - 8 * (gp3.x + 4 * gp2.x * gp.z + 3 * gp.x * gp2.z)
	};

	for (int i = 0; i < 6; ++i)
	{
		out_dn(0, i) = coeffs_x[i] - coeffs_z[i];
		out_dn(1, i) = coeffs_y[i] - coeffs_z[i];
	}

	for (int i = 0; i < 9; ++i)
	{
		out_dn(0, 6 + i) = (coeffs_x[6 + i] - coeffs_z[6 + i]) *tri.tan_multipliers[i];
		out_dn(1, 6 + i) = (coeffs_y[6 + i] - coeffs_z[6 + i]) *tri.tan_multipliers[i];
	}
}

void FEM6Noded::CalcRotationC1(const FETriangle& tri, const std::vector<Vector3>& pos, const std::vector<Vector3>& tans, const Vector3& gaussPoint, float warp_angle, Eigen::Matrix3f& out_rot, JaMatrix& out_ja, BaseMatrix& out_dn)
{
	BuildNaturalCoordinateBasisVectors(tri, gaussPoint, out_dn);
	
	//Build local fill direction (local warp direction)
	Vector3 vFill(-sin(warp_angle), cos(warp_angle), 0.0f);

	//Compute pure deformation in XZ and YZ axis
	Vector3 dir_xz, dir_yz;

	for (int i = 0; i < 6; ++i)
	{
		dir_xz += pos[tri.phyxels[i]] * out_dn(0, i);
		dir_yz += pos[tri.phyxels[i]] * out_dn(1, i);
	}
	for (int i = 0; i < 9; ++i)
	{
		dir_xz += tans[tri.tangents[i]] * out_dn(0, 6 + i);
		dir_yz += tans[tri.tangents[i]] * out_dn(1, 6 + i);
	}

	dir_xz.Normalise();
	dir_yz.Normalise();

	//Build the Rotation Matrix
	Vector3 V_z = Vector3::Cross(dir_xz, dir_yz);  V_z.Normalise();
	Vector3 V_x = Vector3::Cross(vFill, V_z); V_x.Normalise();
	Vector3 V_y = Vector3::Cross(V_z, V_x); V_y.Normalise();

	out_rot(0, 0) = V_x.x; out_rot(0, 1) = V_x.y; out_rot(0, 2) = V_x.z;
	out_rot(1, 0) = V_y.x; out_rot(1, 1) = V_y.y; out_rot(1, 2) = V_y.z;
	out_rot(2, 0) = V_z.x; out_rot(2, 1) = V_z.y; out_rot(2, 2) = V_z.z;


	//	//Build the Jacobian Matrix
	out_ja(0, 0) = Vector3::Dot(dir_xz, V_x);
	out_ja(0, 1) = Vector3::Dot(dir_xz, V_y);
	out_ja(1, 0) = Vector3::Dot(dir_yz, V_x);
	out_ja(1, 1) = Vector3::Dot(dir_yz, V_y);
}

void FEM6Noded::Render_DrawingToVisualiser()
{

	unsigned int i;
	for (i = 0; i < m_NumTotal; ++i)
	{
		if (m_PhyxelIsStatic[i])
			NCLDebug::DrawPoint(m_PhyxelsPos[i], 0.02f, Vector4(1.0f, 0.0f, 0.0f, 0.5f));
		else
			NCLDebug::DrawPoint(m_PhyxelsPos[i], 0.02f, Vector4(0.0f, 0.0f, 0.0f, 0.3f));

	}

	for (i = 0; i < m_NumTriangles; ++i)
	{
		const auto& t = m_Triangles[i];

		Vector3 a = m_PhyxelsPos[t.phyxels[0]];
		Vector3 b = m_PhyxelsPos[t.phyxels[1]];
		Vector3 c = m_PhyxelsPos[t.phyxels[2]];

		Vector3 ab = m_PhyxelsPos[t.phyxels[3]];
		Vector3 bc = m_PhyxelsPos[t.phyxels[4]];
		Vector3 ac = m_PhyxelsPos[t.phyxels[5]];

		int nsteps = 16;
		float step = 1.f / float(nsteps);
		Vector3 inv_light_dir = Vector3(0.5f, 1.0f, -0.8f);
		inv_light_dir.Normalise();
		const Vector4 color = Vector4(0.f, 0.f, 1.f, 1.f);


		Vector3 p1 = a;
		Vector3 p2 = ab;
		Vector3 p3 = b;

		Vector3 ta1 = m_PhyxelsTangent[t.t12];
		Vector3 ta2 = m_PhyxelsTangent[t.t13];
		Vector3 tb1 = m_PhyxelsTangent[t.t21];
		Vector3 tb2 = m_PhyxelsTangent[t.t23];
		Vector3 tc1 = m_PhyxelsTangent[t.t31];
		Vector3 tc2 = m_PhyxelsTangent[t.t32];

		Vector3 tab = m_PhyxelsTangent[t.tab];
		Vector3 tbc = m_PhyxelsTangent[t.tbc];
		Vector3 tac = m_PhyxelsTangent[t.tca];


#if 0
		const int width = 9;
		auto calc_offset = [&](int pidx, int xoffset, int yoffset, float* divisor)
		{
			int x = pidx % width + xoffset;
			int y = pidx / width + yoffset;

			if (x >= 0 && x < width
				&& y >= 0 && y < width)
			{
				int nidx = y * width + x;
				*divisor += 1.f;
				return m_PhyxelsPos[pidx] - m_PhyxelsPos[nidx];
			}

			return Vector3(0, 0, 0);
		};

		uint idxA = t.phyxels[0];
		uint idxB = t.phyxels[1];
		uint idxC = t.phyxels[2];
		uint idxAB = t.phyxels[3];
		uint idxBC = t.phyxels[4];
		uint idxAC = t.phyxels[5];

		float divisor = 0.f;
		if (i % 2 == 0)
		{
			ta1 = (calc_offset(idxA, -1, 1, &divisor) - calc_offset(idxA, 1, -1, &divisor)); ta1 = ta1 / divisor * 2.f; divisor = 0.f;
			ta2 = (calc_offset(idxA, 0, 1, &divisor) - calc_offset(idxA, 0, -1, &divisor)); ta2 = ta2 / divisor * 2.f; divisor = 0.f;
			tb1 = (calc_offset(idxB, 1, -1, &divisor) - calc_offset(idxB, -1, 1, &divisor)); tb1 = tb1 / divisor * 2.f; divisor = 0.f;
			tb2 = (calc_offset(idxB, 1, 0, &divisor) - calc_offset(idxB, -1, 0, &divisor)); tb2 = tb2 / divisor * 2.f; divisor = 0.f;
			tc1 = (calc_offset(idxC, 0, 1, &divisor) - calc_offset(idxC, 0, -1, &divisor)); tc1 = -tc1 / divisor * 2.f; divisor = 0.f;
			tc2 = (calc_offset(idxC, 1, 0, &divisor) - calc_offset(idxC, -1, 0, &divisor)); tc2 = -tc2 / divisor * 2.f; divisor = 0.f;

			tab = (calc_offset(idxAB, 1, 1, &divisor) + calc_offset(idxAB, 0, 1, &divisor) + calc_offset(idxAB, 1, 0, &divisor) - calc_offset(idxAB, -1, 0, &divisor) - calc_offset(idxAB, 0, -1, &divisor) - calc_offset(idxAB, -1, -1, &divisor)); tab = tab / divisor * 2.0f; divisor = 0.f;
			tab = (calc_offset(idxAB, 1, 1, &divisor) - calc_offset(idxAB, -1, -1, &divisor)); tab = tab / divisor; divisor = 0.f;
			tbc = (calc_offset(idxBC, 1, -2, &divisor) - calc_offset(idxBC, -1, 2, &divisor)); tbc = tbc / divisor;	divisor = 0.f;
			tac = (calc_offset(idxAC, -2, 1, &divisor) - calc_offset(idxAC, 2, -1, &divisor)); tac = tac / divisor; divisor = 0.f;
		}
		else
		{
			ta1 = (calc_offset(idxA, -1, 0, &divisor) - calc_offset(idxA, 1, 0, &divisor)); ta1 = ta1 / divisor * 2.f; divisor = 0.f;
			ta2 = (calc_offset(idxA, -1, 1, &divisor) - calc_offset(idxA, 1, -1, &divisor)); ta2 = ta2 / divisor * 2.f; divisor = 0.f;
			tb1 = (calc_offset(idxB, 1, 0, &divisor) - calc_offset(idxB, -1, 0, &divisor)); tb1 = tb1 / divisor * 2.f; divisor = 0.f;
			tb2 = (calc_offset(idxB, 0, 1, &divisor) - calc_offset(idxB, 0, -1, &divisor)); tb2 = tb2 / divisor * 2.f; divisor = 0.f;
			tc1 = (calc_offset(idxC, -1, 1, &divisor) - calc_offset(idxC, 1, -1, &divisor)); tc1 = -tc1 / divisor * 2.f; divisor = 0.f;
			tc2 = (calc_offset(idxC, 0, 1, &divisor) - calc_offset(idxC, 0, -1, &divisor)); tc2 = -tc2 / divisor * 2.f; divisor = 0.f;

			tab = (calc_offset(idxAB, -1, 2, &divisor) - calc_offset(idxAB, 1, -2, &divisor)); tab = tab / divisor; divisor = 0.f;
			tbc = (calc_offset(idxBC, 2, -1, &divisor) - calc_offset(idxBC, -2, 1, &divisor)); tbc = tbc / divisor;	divisor = 0.f;
			tac = (calc_offset(idxAC, -1, -1, &divisor) - calc_offset(idxAC, 1, 1, &divisor)); tac = tac / divisor; divisor = 0.f;
		}
#endif

		auto get_point = [&](const Vector3& gp)
		{
			Vector3 gp2 = gp * gp;
			Vector3 gp3 = gp * gp * gp;
			Vector3 gp4 = gp2 * gp2;
			Vector3 gp5 = gp3 * gp2;

			float gxy = gp.x * gp.y;
			float gx2y = gp2.x * gp.y;
			float gxy2 = gp.x * gp2.y;
			float gx3y = gp3.x * gp.y;
			float gx2y2 = gp2.x * gp2.y;
			float gxy3 = gp.x * gp3.y;

			float FormFunction[]{
				-5 * gp2.x + 14 * gp3.x - 8 * gp4.x + 6 * (gx2y - gx3y - gx2y2),												//PA
				-5 * gp2.y + 14 * gp3.y - 8 * gp4.y + 6 * (gxy2 - gxy3 - gx2y2),												//PB
				-5 * gp2.z + 14 * gp3.z - 8 * gp4.z + 6 * (gp.y*gp2.z - gp.y*gp3.z - gp2.y*gp2.z),								//PC - Note: Y here is mute and can be replaced by 'X'
				-16 * gxy + 48 * (gx2y + gxy2) - 32 * (gx3y + gxy3) - 48 * gx2y2,												//PAB
				-16 * gp.y * gp.z + 48 * (gp2.y*gp.z + gp2.z * gp.y) - 32 * (gp3.z*gp.y + gp.z*gp3.y) - 48 * gp2.y*gp2.z,		//PBC
				-16 * gp.x * gp.z + 48 * (gp2.x*gp.z + gp2.z * gp.x) - 32 * (gp3.z*gp.x + gp.z*gp3.x) - 48 * gp2.x*gp2.z,		//PAC


				0.5 * (-gxy + gx2y + 3 * gxy2) + gx3y - 2 * gx2y2 - gxy3,														//T12
				0.5 * (-gp.x*gp.z + gp2.x*gp.z + 3 * gp.x*gp2.z) + gp3.x*gp.z - 2 * gp2.x*gp2.z - gp.x*gp3.z,					//T13

				0.5 * (-gxy + 3 * gx2y + gxy2) - gx3y - 2 * gx2y2 + gxy3,														//T21
				0.5 * (-gp.y*gp.z + gp2.y*gp.z + 3 * gp.y*gp2.z) + gp3.y*gp.z - 2 * gp2.y*gp2.z - gp.y*gp3.z,					//T23

				0.5 * (gp.x*gp.z - 3 * gp2.x*gp.z - gp.x*gp2.z) + gp3.x*gp.z - gp.x*gp3.z + 2 * gp2.x*gp2.z,					//T31
				0.5 * (gp.y*gp.z - 3 * gp2.y*gp.z - gp.y*gp2.z) + gp3.y*gp.z - gp.y*gp3.z + 2 * gp2.y*gp2.z,					//T32

				-4 * gxy + 12 * (gx2y + gxy2) - 8 * (gx3y + 2 * gx2y2 + gxy3),													//TAB
				-4 * gp.y * gp.z + 12 * (gp2.y * gp.z + gp.y * gp2.z) - 8 * (gp3.y * gp.z + 2 * gp2.y * gp2.z + gp.y * gp3.z),	//TBC
				-4 * gp.x * gp.z + 12 * (gp2.x * gp.z + gp.x * gp2.z) - 8 * (gp3.x * gp.z + 2 * gp2.x * gp2.z + gp.x * gp3.z),	//TAC
			};

			Vector3 p;
			for (int i = 0; i < 6; ++i)
			{
				p += m_PhyxelsPos[t.phyxels[i]] * FormFunction[i];
			}
			for (int i = 0; i < 9; ++i)
			{
				p += m_PhyxelsTangent[t.tangents[i]] * (FormFunction[6 + i] * t.tan_multipliers[i]);
			}

			return p;
		};

	



		float colFactor = 1.f;
		if (Window::GetKeyboard()->KeyDown(KEYBOARD_F))
		{
			colFactor = 0.4f;

			Eigen::Matrix3f rot;
			JaMatrix Ja;
			BaseMatrix Dn;

			const float scalar = 0.03f;
			for (int j = 0; j < 12; ++j)
			{
				Vector3 gp = Vector3(GaussPoint12(j,0), GaussPoint12(j,1), GaussPoint12(j,2));

				Vector3 p = get_point(gp);

				CalcRotationC1(t, m_PhyxelsPos, m_PhyxelsTangent, gp, 0.f, rot, Ja, Dn);

				Vector3 x = Vector3(rot(0, 0), rot(0, 1), rot(0, 2));
				Vector3 y = Vector3(rot(1, 0), rot(1, 1), rot(1, 2));
				Vector3 z = Vector3(rot(2, 0), rot(2, 1), rot(2, 2));

				NCLDebug::DrawHairLine(p, p + x * scalar, Vector4(0.0f, 0.0f, 1.f, 1.f));
				NCLDebug::DrawHairLine(p, p + y * scalar, Vector4(0.0f, 1.0f, 0.f, 1.f));
				NCLDebug::DrawHairLine(p, p + z * scalar, Vector4(1.0f, 0.0f, 0.f, 1.f));
			}

		}
	
		std::function<void(const Vector3& v1, const Vector3& v2, const Vector3& gp)> draw_line = [&](const Vector3& v1, const Vector3& v2, const Vector3& gp)
		{
			NCLDebug::DrawHairLine(v1, v2, Vector4(Vector3(gp.x, gp.y, gp.z) * colFactor + Vector3(0.6, 0.6, 0.6) * (1.f - colFactor), 1.f));
		};

#define TMP_USE_NEW_TAN_BMATRIX 1
#if !TMP_USE_NEW_TAN_BMATRIX
		Eigen::Matrix<float, 18, 1> d_g;
		BMatrix B_nl;
		JaMatrix Ja;
		GMatrix G;
#else
		Eigen::Matrix<float, 45, 1> d_g;
		BMatrix_C1 B_nl;
		JaMatrix Ja;
		Eigen::Matrix<float, 6, 45> G;
#endif


		if (m_RenderMode < 2)
		{
			draw_line = [&](const Vector3& v1, const Vector3& v2, const Vector3& gp)
			{
				for (int j = 0; j < 6; ++j)
				{
					Vector3 pos = (m_PhyxelsPos[t.phyxels[j]] - m_PhyxelsPosInitial[t.phyxels[j]]);

					d_g(j * 3 + 0, 0) = pos.x; d_g(j * 3 + 1, 0) = pos.y; d_g(j * 3 + 2, 0) = pos.z;
				}

				for (int j = 0; j < 9; ++j)
				{
					Vector3 tan = (m_PhyxelsTangent[t.tangents[j]] - m_PhyxelsTangentInitial[t.tangents[j]])* t.tan_multipliers[j];

					d_g(j * 3 + 18, 0) = tan.x; d_g(j * 3 + 19, 0) = tan.y; d_g(j * 3 + 20, 0) = tan.z;
				}

				CalcBMatrix(t, m_PhyxelsPosInitial, m_PhyxelsTangentInitial, Vec3(gp.x, gp.y, gp.z), 0.0f, d_g, B_nl, Ja, G);

				Vec3 strain = B_nl * d_g;
				Vec3 stress = E * strain;
				
				float strain_scalar = strain.norm() * 1.f;
				float stress_scalar = stress.norm() * 0.001f;
			
				Vector3 col = hsv2rgb(Vector3(stress_scalar, 1.f, 1.f));

				if (m_RenderMode == 1)
					col = hsv2rgb(Vector3(strain_scalar, 1.f, 1.f));

				NCLDebug::DrawHairLine(v1, v2, Vector4(col.x, col.y, col.z, 1.f));
			};
		}


		Vector3 gp, p;
		Vector3 xOldP;
		for (int ix = 0; ix <= nsteps; ++ix)
		{

			Vector3 oldp;
			for (int iy = 0; iy <= (nsteps - ix); ++iy)
			{
				gp.x = ix * step;
				gp.y = iy * step;
				gp.z = 1.f - (gp.x + gp.y);

				p = get_point(gp);

				if (iy == nsteps - ix)
				{
					if (ix > 0)
						draw_line(xOldP, p, gp);
					
					xOldP = p;
				}

				if (ix > 0)
				{
					Vector3 tgp = Vector3((ix - 1) * step, iy * step, 0.f);
					tgp.z = 1.f - (tgp.x + tgp.y);
					Vector3 p3 = get_point(tgp);

					draw_line(p, p3, tgp);
					if (iy > 0)draw_line(oldp, p3, tgp);

				}
				
				if (iy > 0)
				{
					draw_line(p, oldp, gp);
				}

				oldp = p;
			}
		}
		


	}
}
