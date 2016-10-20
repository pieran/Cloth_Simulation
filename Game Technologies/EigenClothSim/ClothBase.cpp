#include "ClothBase.h"

ClothBase::ClothBase(ClothIntegrationScheme integrator) : m_NumPhyxels(0)
{
	SetIntegrator(integrator);
}

ClothBase::~ClothBase()
{
	m_NumPhyxels = 0;
	m_NumAllocatedPhyxels = 0;
	m_NumAllocatedTris = 0;
}



uint ClothBase::AddPoint(const uint& triIdx, const Vector3& barycentric_coords)
{
	//Get next available index
	if (m_NumPhyxels < m_NumAllocatedPhyxels)
	{
		uint idxA = m_TriIndicies[triIdx];
		uint idxB = m_TriIndicies[triIdx+1];
		uint idxC = m_TriIndicies[triIdx+2];

		uint nIdx = m_NumPhyxels;
		m_NumPhyxels++;

		m_PhyxelIsStatic[nIdx] = false;

		//Get Deformed Position/Velocity
		m_ClothStatesPos[nIdx] = barycentric_interpolate<Vector3>(m_ClothStatesPos[idxA], m_ClothStatesPos[idxB], m_ClothStatesPos[idxC], barycentric_coords);
		m_ClothStatesVel[nIdx] = barycentric_interpolate<Vector3>(m_ClothStatesVel[idxA], m_ClothStatesVel[idxB], m_ClothStatesVel[idxC], barycentric_coords);

		//Get Initial Position
		m_PhyxelExtForces[nIdx] = barycentric_interpolate<Vector3>(m_PhyxelExtForces[idxA], m_PhyxelExtForces[idxB], m_PhyxelExtForces[idxC], barycentric_coords);

		//Get Mass (WRONG!!!)
		m_PhyxelMass[nIdx] = barycentric_interpolate<float>(m_PhyxelMass[idxA], m_PhyxelMass[idxB], m_PhyxelMass[idxC], barycentric_coords);

		//Get Texture Coordinates
		m_PhyxelTexCoords[nIdx] = barycentric_interpolate<Vector2>(m_PhyxelTexCoords[idxA], m_PhyxelTexCoords[idxB], m_PhyxelTexCoords[idxC], barycentric_coords);
		return nIdx;
	}

	throw new std::exception("OUT_OF_MEMORY");
	return 0;
}

uint ClothBase::AddTriangle(const uint& idxA, const uint& idxB, const uint& idxC)
{
	//Get next available index
	if (m_NumTris < m_NumAllocatedTris)
	{
		uint nIdx = m_NumTris * 3;
		m_NumTris++;

		m_TriIndicies[nIdx] = idxA;
		m_TriIndicies[nIdx + 1] = idxB;
		m_TriIndicies[nIdx + 2] = idxC;
		return nIdx;
	}

	throw new std::exception("OUT_OF_MEMORY");
	return 0;
}

void ClothBase::UpdateTriangleMaterial(const uint& triIdx)
{
	//Stub.
}

void ClothBase::InitializeClothDesign(const ClothDesignEntity<3>* design)
{
	//Reset simulation variables
	m_AccumTime = 0.0f;
	m_TimeStep = 0.0066f;// 0.0016f;

	uint i;


	auto& verts = design->Vertices();
	m_NumPhyxels = verts.size();
	m_NumAllocatedPhyxels = m_NumPhyxels + MAX_ADDITIONAL_PHYXELS;
	ClothBase::AllocateDataStructures(m_NumPhyxels);

	auto& tris = design->Triangles();
	m_NumTris = tris.size();
	m_NumAllocatedTris = m_NumTris + MAX_ADDITIONAL_TRIANGLES;

	auto& quads = design->Quads();
	m_NumQuads = quads.size();



	//Compute the default surface area of the cloth (for mass distribution)
	double totalArea = 0.0;
	for (i = 0; i < m_NumTris; ++i)
	{
		Vec3 a = verts[tris[i].verts[0]].position;
		Vec3 b = verts[tris[i].verts[1]].position;
		Vec3 c = verts[tris[i].verts[2]].position;

		Vec3 e1 = a - c;
		Vec3 e2 = b - c;

		//Area + Mass of Triangle/Phyxels
		float area = e1.cross(e2).norm() * 0.5f;
		totalArea += area;
	}
	m_InitialArea = totalArea;



	//Phyxel Memory
	ClothState clothstate;

	for (uint i = 0; i < m_NumClothStates; ++i)
	{
		ClothBase::GetClothState(0, &clothstate);
		memset(&clothstate.phyxelVel[0].x, 0, m_NumPhyxels * sizeof(Vector3));
		for (i = 0; i < m_NumPhyxels; ++i)
		{
			const Vertex& v = verts[i];
			clothstate.phyxelPos[i] = Vector3(v.position.x(), v.position.y(), v.position.z());
		}
	}

	m_PhyxelForces.resize(m_NumAllocatedPhyxels);
	m_PhyxelExtForces.resize(m_NumAllocatedPhyxels);
	m_PhyxelNormals.resize(m_NumAllocatedPhyxels);
	m_PhyxelIsStatic.resize(m_NumAllocatedPhyxels);
	m_PhyxelMass.resize(m_NumAllocatedPhyxels);

	memset(&m_PhyxelExtForces[0].x, 0, m_NumAllocatedPhyxels * sizeof(Vector3));
	memset(&m_PhyxelNormals[0].x, 0, m_NumAllocatedPhyxels * sizeof(Vector3));

	float mass = (float)((totalArea * mass_density) / (double)(m_NumPhyxels));
	for (i = 0; i < m_NumPhyxels; ++i)
	{
		const Vertex& v = verts[i];
		m_PhyxelIsStatic[i] = v.flags & VERTEXFLAGS_IS_STATIC;
		m_PhyxelMass[i] = mass;

		if (!m_PhyxelIsStatic[i])
		{
			m_PhyxelForces[i] = GRAVITY * float(m_PhyxelMass[i]);
		}
		else
		{
			m_PhyxelForces[i] = Vector3(0.0f, 0.0f, 0.0f);
		}
	}


	//Render-Vertex Allocation
	m_RenderValue.resize(m_NumAllocatedPhyxels);
	m_PhyxelTexCoords.resize(m_NumAllocatedPhyxels);

	m_RenderValueMinMax.x = 0.0f;
	m_RenderValueMinMax.y = 40.0f;

	for (i = 0; i < m_NumPhyxels; ++i)
	{
		const Vertex& v = verts[i];

		m_RenderValue[i] = (float(i * 5.0f) / float(m_NumPhyxels));// ((rand() % 100) / 100.0f) * 5.0f;
		m_PhyxelTexCoords[i] = Vector2(v.texCoord.x(), v.texCoord.y());
	}

	//Tri Memory
	m_TriIndicies.resize(m_NumAllocatedTris * 3);
	memcpy(&m_TriIndicies[0], &tris[0].verts[0], m_NumTris * 3 * sizeof(uint));

	//Quad Memory
	m_QuadIndicies.resize(m_NumQuads * 4);
	memcpy(&m_QuadIndicies[0], &quads[0].vert_a, m_NumQuads * 4 * sizeof(uint));


	GenNormals(&m_ClothStatesPos[0]); //Compute Normals again after final iteration for smooth rendering

}

void ClothBase::UpdateSimulation(float global_timestep)
{
	m_Profiler.ResetAllTiming();
	m_Profiler.BeginTiming(PROFILERID_CLOTH_TOTALTIME);

	m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
	if (m_Callback_OnOuterUpdateScene != nullptr) m_Callback_OnOuterUpdateScene(global_timestep);
	if (m_Callback_OnOuterCollisions != nullptr)  m_Callback_OnOuterCollisions();
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

	m_AccumTime += global_timestep;
	
	switch (m_ClothIntegrationScheme)
	{
	case	CLOTH_INTEGRATION_LEAPFROG:
		for (; m_AccumTime >= m_TimeStep; m_AccumTime -= m_TimeStep)
		{
			StepSimulationLeapFrog();
			if (m_Callback_OnInnerFullTimestep != nullptr) m_Callback_OnInnerFullTimestep(m_TimeStep);
		}
		break;

	case	CLOTH_INTEGRATION_RK4:
		for (; m_AccumTime >= m_TimeStep; m_AccumTime -= m_TimeStep)
		{
			StepSimulationRK4();
			if (m_Callback_OnInnerFullTimestep != nullptr) m_Callback_OnInnerFullTimestep(m_TimeStep);
		}
		break;
	default:
	case 	CLOTH_INTEGRATION_EXPLICIT:
		for (; m_AccumTime >= m_TimeStep; m_AccumTime -= m_TimeStep)
		{
			m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
			if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep, m_NumPhyxels, &m_ClothStatesPos[0], &m_ClothStatesVel[0]);
			if (m_Callback_OnInnerCollisions != nullptr)  m_Callback_OnInnerCollisions();
			m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

			StepSimulationExplicit(m_TimeStep, &m_ClothStatesPos[0], &m_ClothStatesVel[0], &m_ClothStatesPos[0], &m_ClothStatesVel[0]);
			if (m_Callback_OnInnerFullTimestep != nullptr) m_Callback_OnInnerFullTimestep(m_TimeStep);
		}
		break;
	}
	
	if (m_Callback_OnInnerCollisions != nullptr)  m_Callback_OnInnerCollisions();
	GenNormals(&m_ClothStatesPos[0]); //Compute Normals again after final iteration for smooth rendering

	m_Profiler.EndTiming(PROFILERID_CLOTH_TOTALTIME);
}

void ClothBase::StepSimulationExplicit(float timestep, const Vector3* in_clothstatepos, const Vector3* in_clothstatevel, Vector3* out_clothstatepos, Vector3* out_clothstatevel)
{
	m_Profiler.BeginTiming(PROFILERID_CLOTH_BUILDROTATIONANDNORMALS);
	GenNormals(in_clothstatepos);
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_BUILDROTATIONANDNORMALS);

	m_Profiler.BeginTiming(PROFILERID_CLOTH_COMPUTEFORCES);
	ApplyAirResistance(in_clothstatepos, in_clothstatevel);
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_COMPUTEFORCES);

	ComputeVelocity(timestep, in_clothstatepos, in_clothstatevel, out_clothstatevel);

	memset(&m_PhyxelExtForces[0].x, 0, m_NumPhyxels * sizeof(Vector3));

	if (out_clothstatepos != NULL)
	{
#pragma omp parallel for
		for (int i = 0; i < (int)m_NumPhyxels; ++i)
		{
			out_clothstatepos[i] = in_clothstatepos[i] + out_clothstatevel[i] * timestep;
		}
	}
}

void ClothBase::AdvanceVelocityHalfTimestep()
{
	if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep * 0.5f, m_NumPhyxels, &m_ClothStatesPos[0], &m_ClothStatesVel[0]);
	
	StepSimulationExplicit(m_TimeStep, &m_ClothStatesPos[0], &m_ClothStatesVel[0], NULL, &m_ClothStatesVel[0]);
	
	if (m_Callback_OnInnerFullTimestep != nullptr) m_Callback_OnInnerFullTimestep(m_TimeStep * 0.5f);
}

void ClothBase::StepSimulationLeapFrog()
{
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumPhyxels; ++i)
	{
		m_ClothStatesPos[i] += m_ClothStatesVel[i] * m_TimeStep;
	}


	m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
	if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep * 0.5f, m_NumPhyxels, &m_ClothStatesPos[0], &m_ClothStatesVel[0]);
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

	StepSimulationExplicit(m_TimeStep, &m_ClothStatesPos[0], &m_ClothStatesVel[0], NULL, &m_ClothStatesVel[0]);
}

void ClothBase::StepSimulationRK4()
{
	//memcpy(&m_PhyxelPos[m_NumPhyxels], &m_PhyxelPos[0], m_NumPhyxels * sizeof(Vector3));
	Vector3* av = &m_ClothStatesVel[0];
	Vector3* bv = &m_ClothStatesVel[m_NumPhyxels];
	Vector3* cv = &m_ClothStatesVel[m_NumPhyxels * 2];
	Vector3* dv = &m_ClothStatesVel[m_NumPhyxels * 3];

	Vector3* ap = &m_ClothStatesPos[0];
	Vector3* bp = &m_ClothStatesPos[m_NumPhyxels];
	Vector3* cp = &m_ClothStatesPos[m_NumPhyxels * 2];
	Vector3* dp = &m_ClothStatesPos[m_NumPhyxels * 3];

	//a = evaluate(state, t, 0.0f, Derivative());
	//a = last velocity	
		
	memcpy(av, bv, m_NumPhyxels * sizeof(Vector3));
	memcpy(av, cv, m_NumPhyxels * sizeof(Vector3));
	memcpy(av, dv, m_NumPhyxels * sizeof(Vector3));


/*	double	k1 = f(x, y),
		k2 = f(x + m_TimeStep / 2, y + k1 * m_TimeStep / 2),
		k3 = f(x + m_TimeStep / 2, y + k2 * m_TimeStep / 2),
		k4 = f(x + m_TimeStep, y + k3 * m_TimeStep);
	return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6 * m_TimeStep;


	av = calcVel(oldTime, oldPos);
	bv = calcVel(oldTime + m_TimeStep * 0.5f, oldPos + av * m_TimeStep * 0.5f);
	cv = calcVel(oldTime + m_TimeStep * 0.5f, oldPos + bv * m_TimeStep * 0.5f);
	dv = calcVel(oldTime + m_TimeStep, oldPos + cv * m_TimeStep);*/


	//av = calcVel(oldTime, oldPos);
	//StepSimulationExplicit(0.0f, ap, av, NULL, av);
	
	//bv = calcVel(oldTime + m_TimeStep * 0.5f, oldPos + av * m_TimeStep * 0.5f);
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumPhyxels; ++i)
	{
		bp[i] = ap[i]  +av[i] * (m_TimeStep * 0.5f);
	}

	m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
	if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep * 0.5f, m_NumPhyxels, bp, bv);
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

	StepSimulationExplicit(m_TimeStep * 0.5f, bp, bv, NULL, bv);



	//cv = calcVel(oldTime + m_TimeStep * 0.5f, oldPos + bv * m_TimeStep * 0.5f);
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumPhyxels; ++i)
	{
		cp[i] = ap[i] +bv[i] * (m_TimeStep * 0.5f);
	}

	m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
	if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep * 0.5f, m_NumPhyxels, cp, cv);
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

	StepSimulationExplicit(m_TimeStep * 0.5f, cp, cv, NULL, cv);



	//dv = calcVel(oldTime + m_TimeStep, oldPos + cv * m_TimeStep);
#pragma omp parallel for
	for (int i = 0; i < (int)m_NumPhyxels; ++i)
	{
		dp[i] = ap[i] +cv[i] * m_TimeStep;
	}

	m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
	if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep, m_NumPhyxels, dp, dv);
	m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

	StepSimulationExplicit(m_TimeStep, dp, dv, NULL, dv);



	/*
			//b = evaluate(state, t, dt*0.5f, a);
	#pragma omp parallel for
			for (int i = 0; i < (int)m_NumPhyxels; ++i)
			{
				bp[i] = ap[i] + av[i] * (m_TimeStep * 0.5f);
			}
			memcpy(bv, av, m_NumPhyxels * sizeof(Vector3));

			m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
			if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep * 0.5f, m_NumPhyxels, bp, bv);
			if (m_Callback_OnInnerCollisions != nullptr)  m_Callback_OnInnerCollisions();
			m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

			StepSimulationExplicit(m_TimeStep * 0.5f, bp, bv, NULL, bv);

	
	#pragma omp parallel for
			for (int i = 0; i < (int)m_NumPhyxels; ++i)
			{
				cp[i] = ap[i] + bv[i] * (m_TimeStep * 0.5f);
			}

			memcpy(cv, bv, m_NumPhyxels * sizeof(Vector3));
			m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
			if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep * 0.5f, m_NumPhyxels, cp, cv);
			if (m_Callback_OnInnerCollisions != nullptr)  m_Callback_OnInnerCollisions();
			m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

			//c = evaluate(state, t, dt*0.5f, b);
			StepSimulationExplicit(m_TimeStep * 0.5f, cp, cv, NULL, cv);
	

	#pragma omp parallel for
			for (int i = 0; i < (int)m_NumPhyxels; ++i)
			{
				dp[i] = cp[i] - cv[i] * m_TimeStep;
			}
			memcpy(dv, cv, m_NumPhyxels * sizeof(Vector3));
			m_Profiler.BeginTiming(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
			if (m_Callback_OnInnerApplySubTimestep != nullptr) m_Callback_OnInnerApplySubTimestep(m_TimeStep, m_NumPhyxels, dp, dv);
			if (m_Callback_OnInnerCollisions != nullptr)  m_Callback_OnInnerCollisions();
			m_Profiler.EndTimingAccumulative(PROFILERID_CLOTH_EXTERNALCOLLISIONS);
	
			//d = evaluate(state, t, dt, c);
			StepSimulationExplicit(m_TimeStep, dp, dv, NULL, dv);
	*/

	#pragma omp parallel for
			for (int i = 0; i < (int)m_NumPhyxels; ++i)
			{
				Vector3 dxdt = (av[i] + (bv[i] + cv[i]) * 2.0f + dv[i]) / 6.0f;
	
				av[i] = dxdt;
				ap[i] = ap[i] + dxdt * m_TimeStep;
			}
}

void ClothBase::SetIntegrator(ClothIntegrationScheme integrator)
{
	if (integrator != m_ClothIntegrationScheme && m_NumPhyxels > 0)
	{
		m_ClothIntegrationScheme = integrator;
		AllocateDataStructures(m_NumPhyxels);
	}

	m_ClothIntegrationScheme = integrator;
}

bool ClothBase::GetClothState(uint idx, ClothState* out_clothstate)
{
	if (idx >= m_NumClothStates)
	{
		return false;
	}

	out_clothstate->numPhyxels = m_NumPhyxels;
	out_clothstate->phyxelPos = &m_ClothStatesPos[m_NumPhyxels * idx];
	out_clothstate->phyxelVel = &m_ClothStatesVel[m_NumPhyxels * idx];

	return true;
}

void ClothBase::AllocateDataStructures(uint numPhyxels)
{
	m_NumPhyxels = numPhyxels;
	m_NumAllocatedPhyxels = m_NumPhyxels + MAX_ADDITIONAL_PHYXELS;
	switch (m_ClothIntegrationScheme)
	{
	case	CLOTH_INTEGRATION_LEAPFROG:
		m_NumClothStates = 2;
		break;

	case	CLOTH_INTEGRATION_RK4:
		m_NumClothStates = 4;
		break;
	default:
	case 	CLOTH_INTEGRATION_EXPLICIT:
		m_NumClothStates = 1;
		break;
	}

	uint totalLength = m_NumAllocatedPhyxels * m_NumClothStates;
	m_ClothStatesPos.resize(totalLength);
	m_ClothStatesVel.resize(totalLength);
}


void ClothBase::ApplyAirResistance(const Vector3* in_clothpos, const Vector3* in_clothvel)
{
//#pragma omp parallel for
	for (int i = 0; i < (int)m_NumTris; ++i)
	{
		const int i3 = i * 3;
		uint idxA = m_TriIndicies[i3];
		uint idxB = m_TriIndicies[i3 + 1];
		uint idxC = m_TriIndicies[i3 + 2];


		Vector3 triNormals[3] = { m_PhyxelNormals[idxA], m_PhyxelNormals[idxB], m_PhyxelNormals[idxC] };
		Vector3 triPositions[3] = { in_clothpos[idxA], in_clothpos[idxB], in_clothpos[idxC] };

		Vector3 face_normal = (triNormals[0] + triNormals[1] + triNormals[2]) / 3.0f;
		face_normal.Normalise();

		float face_area = Vector3::Cross(triPositions[0] - triPositions[2], triPositions[1] - triPositions[2]).Length() * 0.5f;




		const float Cd = 0.55f; //Drag Coefficient
		const float Cl = 0.0f; //Lift Coefficient
		const float p = 1.293f; //Air Density


								//Vrel = Vi - u (Where u is velocity field of the wind)
								//F = 0.5f * Cd * p * |Vrel|^2 * A * (n . Vrel) . (-Vrel)
		Vector3 vrel = (in_clothvel[idxA] + in_clothvel[idxB] + in_clothvel[idxC]) / 3.0f - WIND;
		Vector3 vrel_norm = vrel; vrel_norm.Normalise();

		if (Vector3::Dot(face_normal, vrel) < 0) face_normal = -face_normal;

		Vector3 lift_dir = Vector3::Cross(Vector3::Cross(face_normal, vrel_norm), vrel_norm);

		float vrelsq = Vector3::Dot(vrel, vrel);
		float drag_coef = 0.5f * Cd * p * vrelsq * face_area * Vector3::Dot(face_normal, vrel_norm);
		float lift_coef = 0.5f * Cl * p * vrelsq * face_area * Vector3::Dot(lift_dir, vrel_norm);

		Vector3 Fid = (-vrel_norm * drag_coef) / 3.0f;		//Drag
		Vector3 Fil = (lift_dir * lift_coef) / 3.0f;		//Lift

		Vector3 totalForce = (Fid + Fil);

		m_PhyxelExtForces[idxA] += totalForce;
		m_PhyxelExtForces[idxB] += totalForce;
		m_PhyxelExtForces[idxC] += totalForce;

		/*Vector3* B_i = &m_TriB[i3];
		B_i[0] = totalForce;
		B_i[1] = totalForce;
		B_i[2] = totalForce;*/
	}
}

void ClothBase::GenNormals(const Vector3* in_clothpos)
{
	memset(&m_PhyxelNormals[0], 0, m_NumPhyxels * sizeof(Vector3));
#if USE_QUAD_NORMALS
	for (int i = 0; i < (int)m_NumQuads; ++i)
	{
		uint vert_a = m_QuadIndicies[i * 4 + 0];
		uint vert_b = m_QuadIndicies[i * 4 + 1];
		uint vert_c = m_QuadIndicies[i * 4 + 2];
		uint vert_d = m_QuadIndicies[i * 4 + 3];

		Vector3 a = in_clothpos[vert_a];
		Vector3 b = in_clothpos[vert_b];
		Vector3 c = in_clothpos[vert_c];
		Vector3 d = in_clothpos[vert_d];

		Vector3 normalA = Vector3::Cross(a - c, b - c);
		Vector3 normalB = Vector3::Cross(d - c, a - c);
		normalA.Normalise();
		normalB.Normalise();

		Vector3 normal = normalA + normalB;
		normal.Normalise();

		m_PhyxelNormals[vert_a] += normal;
		m_PhyxelNormals[vert_b] += normal;
		m_PhyxelNormals[vert_c] += normal;
		m_PhyxelNormals[vert_d] += normal;
	}
#else
	for (int i = 0; i < (int)m_NumTris; ++i)
	{
		uint vert_a = m_TriIndicies[i * 3 + 0];
		uint vert_b = m_TriIndicies[i * 3 + 1];
		uint vert_c = m_TriIndicies[i * 3 + 2];

		Vector3 a = in_clothpos[vert_a];
		Vector3 b = in_clothpos[vert_b];
		Vector3 c = in_clothpos[vert_c];

		Vector3 normal = Vector3::Cross(a - c, b - c);
		normal.Normalise();

		m_PhyxelNormals[vert_a] += normal;
		m_PhyxelNormals[vert_b] += normal;
		m_PhyxelNormals[vert_c] += normal;
	}
#endif

#pragma omp parallel for
	for (int i = 0; i < (int)m_NumPhyxels; ++i)
	{
		m_PhyxelNormals[i].Normalise();
	}
}

bool ClothBase::GetClosestPointsForCutting(const Ray& start, const Ray& end, int* out_idx_start, int* out_idx_end)
{
	std::vector<int> possible_start_points, possible_end_points;

	const float limit_distance = 0.2f;

	/*float dist_start = FLT_MAX;
	float dist_end = FLT_MAX;
	int idx_start = -1;
	int idx_end = -1;

	float t, cdist_start, cdist_end;
	Vector3 pos, proj;
	for (int i = 0; i < (int)m_NumPhyxels; ++i) //TODO: Implement multi-threaded variation with reduction (Or use Spatial subdivision schemes)
	{
		//Calculate Distance to Ray
		pos = m_ClothStatesPos[i];

		t = Vector3::Dot(pos - start.wsPos, start.wsDir);
		cdist_start = (pos - start.wsPos - start.wsDir * t).Length();
		if (cdist_start < limit_distance)
		{
			possible_start_points.push_back(i);
		}

		t = Vector3::Dot(pos - end.wsPos, end.wsDir);
		cdist_end = (pos - end.wsPos - end.wsDir * t).Length();
		if (cdist_end < limit_distance)
		{
			possible_end_points.push_back(i);
		}
	}


	for (int i = 0; i < (int)possible_start_points.size(); ++i)
	{
		pos = m_ClothStatesPos[possible_start_points[i]];
		cdist_start = (start.wsPos - pos).LengthSquared();

		if (cdist_start < dist_start)
		{
			dist_start = cdist_start;
			idx_start = possible_start_points[i];
		}
	}

	for (int i = 0; i < (int)possible_end_points.size(); ++i)
	{
		pos = m_ClothStatesPos[possible_end_points[i]];
		cdist_end = (end.wsPos - pos).LengthSquared();

		if (cdist_end < dist_start)
		{
			dist_start = cdist_end;
			idx_end = possible_end_points[i];
		}
	}


	if (out_idx_start) *out_idx_start = idx_start;
	if (out_idx_end) *out_idx_end = idx_end;

	return (idx_start > -1 && idx_end > -1);*/


	float dist_start = FLT_MAX;
	float dist_end = FLT_MAX;
	int idx_start = -1;
	int idx_end = -1;

	float t, cdist_start, cdist_end;
	Vector3 pos, proj;
	for (int i = 0; i < (int)m_NumPhyxels; ++i) //TODO: Implement multi-threaded variation with reduction (Or use Spatial subdivision schemes)
	{
		//Calculate Distance to Ray
		pos = m_ClothStatesPos[i];
	
		t = Vector3::Dot(pos - start.wsPos, start.wsDir);
		cdist_start = (pos - start.wsPos - start.wsDir * t).LengthSquared();
		cdist_start += (start.wsPos - pos).LengthSquared() * 0.02f;
		if (cdist_start < dist_start)
		{
			dist_start = cdist_start;
			idx_start = i;
		}
		

		t = Vector3::Dot(pos - end.wsPos, end.wsDir);
		cdist_end = (pos - end.wsPos - end.wsDir * t).LengthSquared();
		cdist_end += (end.wsPos - pos).LengthSquared() * 0.02f;
		if (cdist_end < dist_end)
		{
			dist_end = cdist_end;
			idx_end = i;
		}
	}

	if (out_idx_start) *out_idx_start = idx_start;
	if (out_idx_end) *out_idx_end = idx_end;

	return (idx_start > -1 && idx_end > -1);
}