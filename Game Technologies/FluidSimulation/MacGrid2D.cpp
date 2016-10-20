#include "MacGrid2D.h"
#include <nclgl\Matrix4.h>

MacGrid2D::MacGrid2D(const BoundingBox2D& region, uint div_x, uint div_y)
	: m_BoundingRegion(region)
	, m_DivX(div_x)
	, m_DivY(div_y)
	, m_ExternalPressure(0.0f)
	, m_ExternalVelocityX(0.0f)
	, m_ExternalVelocityY(0.0f)
	, m_Solver()
{
	m_Solver.AllocateMemory(m_DivX * m_DivY);
	m_Pressure.resize(m_DivX * m_DivY);
	m_VelocityX.resize((m_DivX + 1) * m_DivY);
	m_VelocityY.resize(m_DivX * (m_DivY + 1));


	memset(&m_Pressure[0], 0, m_Pressure.size() * sizeof(float));
	memset(&m_VelocityX[0], 0, m_VelocityX.size() * sizeof(float));
	memset(&m_VelocityY[0], 0, m_VelocityY.size() * sizeof(float));

	m_CellDimensions.x = m_BoundingRegion.scale.x / (float)m_DivX;
	m_CellDimensions.y = m_BoundingRegion.scale.y / (float)m_DivY;



	m_IntermediatePressure.resize(m_DivX * m_DivY);
	m_IntermediateVelocityX.resize((m_DivX+1) * m_DivY);
	m_IntermediateVelocityY.resize(m_DivX * (m_DivY+1));

	Vector2 center;
	center.x = m_BoundingRegion.offset_min.x + m_BoundingRegion.scale.x * 0.5f;
	center.y = m_BoundingRegion.offset_min.y + m_BoundingRegion.scale.y * 0.5f;


	const float speed = 0.1f;
	const float max_len = m_BoundingRegion.scale.Length() * 0.5f;
#pragma omp parallel for

	for (int iy = 0; iy < (int)(m_DivY + 1); ++iy)
	{
		for (int ix = 0; ix < (int)(m_DivX + 1); ++ix)
		{
			Vector2 wsPos = GetCellWorldSpace(ix, iy);
			Vector2 relPos = wsPos - center;
			float dist = relPos.Length();
			
			//SetPressure(ix, iy, 1.0f - dist / max_len);
			Vector3 rot = Vector3(-relPos.y, relPos.x, 0.0f);

			//rot.Normalise();
			rot *= speed;

			SetVelocityX(ix, iy, rot.x);
			//SetVelocityY(ix, iy, rot.y);
		}
	}
	

}

MacGrid2D::~MacGrid2D()
{

}



float MacGrid2D::CalcMaxVelocity()
{
	const int len = m_DivX * m_DivY;
	float velSq_max = 0.0f;

	//#pragma omp parallel for reduction(max : velSq_max)

	for (int iy = 0; iy < (int)m_DivY; ++iy)
	{
		for (int ix = 0; ix < (int)m_DivX; ++ix)
		{
			//Trilinear Interpolation
			Vector2 vel = GetVelocityMidPoint(ix, iy);

			float vel_squared = vel.LengthSquared();
			if (vel_squared > velSq_max)
				velSq_max = vel_squared;
		}
	}
	

	return sqrt(velSq_max);
}

void MacGrid2D::ApplyGravity(float sub_timestep)
{
	const int size = m_DivX * (m_DivY);

#pragma omp parallel for
	for (int iy = 0; iy < (int)m_DivY; ++iy)
	{
		for (int ix = 0; ix < (int)m_DivX; ++ix)
		{
			if (m_Pressure[iy * m_DivX + ix] > 0.0f)
				m_VelocityY[iy * m_DivX + ix] -= 9.81f * sub_timestep;
			else
				m_VelocityY[iy * m_DivX + ix] = 0.0f;
		}
	}
}

void MacGrid2D::AdvectFluid(float sub_timestep)
{
#define ADVECT_INDIVIDUAL_COMPONENTS FALSE

	#pragma omp parallel for
	for (int iy = 0; iy < (int)m_DivY; ++iy)
	{
		for (int ix = 0; ix < (int)m_DivX; ++ix)
		{
			uint cellIdx = iy * m_DivX + ix;

			//Note <var>0 corresponds to the cell at current time, and <var>1 corresponds to cell at current time - delta time (old cell from which to advect)
			Vector2 vel0 = GetVelocityMidPoint(ix, iy);

			//float pressure = m_Pressure[(iz * m_DivY + iy) * m_DivX + ix];

			Vector2 wsPos0 = GetCellWorldSpace(ix, iy);
			Vector2 wsPos1 = wsPos0 - vel0 * sub_timestep;

			m_IntermediatePressure[cellIdx] = GetPressure(wsPos1);

#if (ADVECT_INDIVIDUAL_COMPONENTS == FALSE)
			Vector2 vel1 = GetVelocityTrilinear(wsPos1);
			m_IntermediateVelocityX[cellIdx] = vel1.x;
			m_IntermediateVelocityY[cellIdx] = vel1.y;
#endif
		}
	}


#if (ADVECT_INDIVIDUAL_COMPONENTS == TRUE)
	//Advect Velocity X
#pragma omp parallel for
	for (int iy = 0; iy < (int)m_DivY; ++iy)
	{
		for (int ix = 0; ix < (int)m_DivX + 1; ++ix)
		{
			Vector2 vel0;
			vel0.x = GetVelocityX(ix, iy);
			vel0.y = (GetVelocityY(ix, iy - 1) + GetVelocityY(ix, iy-1) + GetVelocityY(ix+1, iy) + GetVelocityY(ix+1, iy-1)) * 0.25f;
		
			Vector2 wsPos0;
			wsPos0.x = m_CellDimensions.x * ((float)ix) + m_BoundingRegion.offset_min.x;
			wsPos0.y = m_CellDimensions.y * ((float)iy + 0.5f) + m_BoundingRegion.offset_min.y;

			Vector2 wsPos1 = wsPos0 - vel0 * sub_timestep;

			Vector2 vel1 = GetVelocityTrilinear(wsPos1);
			m_IntermediateVelocityX[iy * (m_DivX + 1) + ix] = vel1.x;
		}
	}

	//Advect Velocity Y
#pragma omp parallel for
	for (int iy = 0; iy < (int)m_DivY + 1; ++iy)
	{
		for (int ix = 0; ix < (int)m_DivX; ++ix)
		{
			Vector2 vel0;
			vel0.x = (GetVelocityX(ix - 1, iy) + GetVelocityX(ix - 1, iy) + GetVelocityX(ix, iy + 1) + GetVelocityX(ix - 1, iy + 1)) * 0.25f;
			vel0.y = GetVelocityY(ix, iy);

			Vector2 wsPos0;
			wsPos0.x = m_CellDimensions.x * ((float)ix + 0.5f) + m_BoundingRegion.offset_min.x;
			wsPos0.y = m_CellDimensions.y * ((float)iy) + m_BoundingRegion.offset_min.y;

			Vector2 wsPos1 = wsPos0 - vel0 * sub_timestep;

			Vector2 vel1 = GetVelocityTrilinear(wsPos1);
			m_IntermediateVelocityY[iy * m_DivX + ix] = vel1.y;
		}
	}
#endif







	const int padx = 0;
	const int pady =0;

	//Set Pressure
#pragma omp parallel for
	for (int iy = pady; iy < (int)m_DivY - pady; ++iy)
	{
		for (int ix = padx; ix < (int)m_DivX - padx; ++ix)
		{
			uint cellIdx = iy * m_DivX + ix;
			m_Pressure[cellIdx] = m_IntermediatePressure[cellIdx];
		}
	}

	
	//Set Velocity X
	#pragma omp parallel for
	for (int iy = pady; iy < (int)m_DivY - pady; ++iy)
	{
		for (int ix = padx; ix < (int)m_DivX + 1 - padx; ++ix)
		{
#if (ADVECT_INDIVIDUAL_COMPONENTS == TRUE)
			uint cellIdx = iy * (m_DivX+1) + ix;
			m_VelocityX[cellIdx] = m_IntermediateVelocityX[cellIdx];
#else
			float vel0 = CellExists(ix - 1, iy) ? m_IntermediateVelocityX[iy * m_DivX + ix - 1] : m_ExternalVelocityX;
			float vel1 = CellExists(ix, iy) ? m_IntermediateVelocityX[iy * m_DivX + ix] : m_ExternalVelocityX;

			SetVelocityX(ix, iy, (vel0 + vel1) * 0.5f);
#endif
		}
	}
	

	//Set Velocity Y
	#pragma omp parallel for
	for (int iy = pady; iy < (int)m_DivY + 1 - pady; ++iy)
	{
		for (int ix = padx; ix < (int)m_DivX - padx; ++ix)
		{
#if (ADVECT_INDIVIDUAL_COMPONENTS == TRUE)
			uint cellIdx = iy * m_DivX + ix;
			m_VelocityY[cellIdx] = m_IntermediateVelocityY[cellIdx];
#else
			float vel0 = CellExists(ix, iy - 1) ? m_IntermediateVelocityY[(iy - 1) * m_DivX + ix] : m_ExternalVelocityY;
			float vel1 = CellExists(ix, iy) ? m_IntermediateVelocityY[iy * m_DivX + ix] : m_ExternalVelocityY;

			SetVelocityY(ix, iy, (vel0 + vel1) * 0.5f);
#endif
		}
	}

}

void MacGrid2D::ProjectPressure(float sub_timestep)
{
	const float total_pressure = 1.0f;	//Change in fluid volume??
	const float dtp = total_pressure / sub_timestep; 
	const float dx = 1.0f; //Change in global velocity??
	const float invdx = 1.0f / dx;

	const float default_a = -dtp * invdx * invdx;

	//Build Matrix
#pragma omp parallel for
	for (int iy = 0; iy < (int)m_DivY - 0; ++iy)
	{
		for (int ix = 0; ix < (int)m_DivX - 0; ++ix)
		{
			uint idx = iy * m_DivX + ix;
			m_Solver.m_A[idx][idx] = 4.0f * invdx * invdx * dtp;
			m_Solver.m_A[idx][iy * m_DivX + ix + 1]   = default_a;
			m_Solver.m_A[idx][iy * m_DivX + ix - 1]   = default_a;
			m_Solver.m_A[idx][(iy + 1) * m_DivX + ix] = default_a;
			m_Solver.m_A[idx][(iy - 1) * m_DivX + ix] = default_a;


			m_Solver.m_B[idx] = -(GetVelocityX(ix + 1, iy) - GetVelocityX(ix, iy) + GetVelocityY(ix, iy + 1) - GetVelocityY(ix, iy)) * invdx;
		}
	}

	m_Solver.SetTolerance(1E-3f);
	m_Solver.SetMaxIterations(100);
	m_Solver.Solve();

	int num_total = (int)(m_DivX * m_DivY);
#pragma omp parallel for
	for (int i = 0; i < num_total; ++i)
	{
		int ix = i % m_DivX;
		int iy = (i) / m_DivX;
		//if (ix > 0 && ix < m_DivX - 1 && iy > 0 && iy < m_DivY -1)
			//m_Pressure[i] -= m_Solver.m_B[i] * sub_timestep;

		//m_VelocityX[iy * (m_DivX+1) + ix]	-= m_Solver.m_B[i] / dtp;
		//m_VelocityY[iy * m_DivX + ix]		-= m_Solver.m_B[i] / dtp;
	}



	//Set Velocity X
#pragma omp parallel for
	for (int iy = 0; iy < (int)m_DivY; ++iy)
	{
		for (int ix = 0; ix < (int)m_DivX + 1; ++ix)
		{
			float pd0 = (ix > 0) ? m_Solver.m_B[iy * m_DivX + ix - 1] : m_ExternalPressure;
			float pd1 = m_Solver.m_B[iy * m_DivX + ix];

			float combined = (pd0 + pd1) * 0.5f;

			m_VelocityX[iy * (m_DivX + 1) + ix] -= combined / dtp;
		}
	}


	//Set Velocity Y
#pragma omp parallel for
	for (int iy = 0; iy < (int)m_DivY + 1; ++iy)
	{
		for (int ix = 0; ix < (int)m_DivX; ++ix)
		{
			float pd0 = (iy > 0) ? m_Solver.m_B[(iy-1) * m_DivX + ix] : m_ExternalPressure;
			float pd1 = m_Solver.m_B[iy * m_DivX + ix];

			float combined = (pd0 + pd1) * 0.5f;

			m_VelocityY[iy * m_DivX + ix] -= combined / dtp;
		}
	}

}