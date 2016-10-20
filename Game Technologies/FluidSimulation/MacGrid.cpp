#include "MacGrid.h"
#include <nclgl\Matrix4.h>

MacGrid::MacGrid(const BoundingBox& region, uint div_x, uint div_y, uint div_z)
	: m_BoundingRegion(region)
	, m_DivX(div_x)
	, m_DivY(div_y)
	, m_DivZ(div_z)
	, m_ExternalPressure(0.0f)
	, m_ExternalVelocityX(0.0f)
	, m_ExternalVelocityY(0.0f)
	, m_ExternalVelocityZ(0.0f)
{

	m_Pressure.resize(m_DivX * m_DivY * m_DivZ);
	m_VelocityX.resize((m_DivX + 1) * m_DivY * m_DivZ);
	m_VelocityY.resize(m_DivX * (m_DivY + 1) * m_DivZ);
	m_VelocityZ.resize(m_DivX * m_DivY * (m_DivZ + 1));

	memset(&m_Pressure[0], 0, m_Pressure.size() * sizeof(float));
	memset(&m_VelocityX[0], 0, m_VelocityX.size() * sizeof(float));
	memset(&m_VelocityY[0], 0, m_VelocityY.size() * sizeof(float));
	memset(&m_VelocityZ[0], 0, m_VelocityZ.size() * sizeof(float));

	m_CellDimensions.x = m_BoundingRegion.scale.x / (float)m_DivX;
	m_CellDimensions.y = m_BoundingRegion.scale.y / (float)m_DivY;
	m_CellDimensions.z = m_BoundingRegion.scale.z / (float)m_DivZ;



	m_IntermediatePressure.resize(m_DivX * m_DivY * m_DivZ);
	m_IntermediateVelocityX.resize(m_DivX * m_DivY * m_DivZ);
	m_IntermediateVelocityY.resize(m_DivX * m_DivY * m_DivZ);
	m_IntermediateVelocityZ.resize(m_DivX * m_DivY * m_DivZ);

	Vector3 center;
	center.x = m_BoundingRegion.offset_min.x + m_BoundingRegion.scale.x * 0.5f;
	center.y = m_BoundingRegion.offset_min.y + m_BoundingRegion.scale.y * 0.5f;
	center.z = m_BoundingRegion.offset_min.z + m_BoundingRegion.scale.z * 0.5f;

	Matrix4 rot90 = Matrix4::Rotation(90.0f, Vector3(1, 0, 0));

	const float speed = 1.1f;
#pragma omp parallel for
	for (int iz = 0; iz < (int)(m_DivZ+1); ++iz)
	{
		for (int iy = 0; iy < (int)(m_DivY+1); ++iy)
		{
			for (int ix = 0; ix < (int)(m_DivX+1); ++ix)
			{
				Vector3 wsPos = GetCellWorldSpace(ix, iy, iz);
				Vector3 relPos = wsPos - center;

				Vector3 rotX = Vector3(0.0f, -relPos.z, relPos.y);
				Vector3 rotY = Vector3(-relPos.z, 0.0f, relPos.x);
				Vector3 rotZ = Vector3(-relPos.y, relPos.x, 0.0f);

				Vector3 tangent = rotZ;

				tangent.Normalise();
				tangent *= speed;

				SetVelocityX(ix, iy, iz, tangent.x);
				SetVelocityY(ix, iy, iz, tangent.y);
				SetVelocityZ(ix, iy, iz, tangent.z);


				/*SetVelocityX(ix, iy, iz, ix / (float)(m_DivX));
				SetVelocityY(ix, iy, iz, iy / (float)(m_DivY));
				SetVelocityZ(ix, iy, iz, 0.0f);*/
			}
		}
	}

}

MacGrid::~MacGrid()
{

}



float MacGrid::CalcMaxVelocity()
{
	const int len = m_DivX * m_DivY * m_DivZ;
	float velSq_max = 0.0f;

//#pragma omp parallel for reduction(max : velSq_max)
	for (int iz = 0; iz < (int)m_DivZ; ++iz)
	{
		for (int iy = 0; iy < (int)m_DivY; ++iy)
		{
			for (int ix = 0; ix < (int)m_DivX; ++ix)
			{
				//Trilinear Interpolation
				Vector3 vel = GetVelocityMidPoint(ix, iy, iz);

				float vel_squared = vel.LengthSquared();
				if (vel_squared > velSq_max)
					velSq_max = vel_squared;
			}
		}
	}

	return sqrt(velSq_max);
}

void MacGrid::ApplyGravity(float sub_timestep)
{
	/*const int size = m_DivX * (m_DivY+1) * m_DivZ;

#pragma omp parallel for
	for (int i = 0; i < size; ++i)
	{
		m_VelocityY[i] -= 9.81f * sub_timestep;
	}*/
}

void MacGrid::AdvectFluid(float sub_timestep)
{
//#pragma omp parallel for
	for (int iz = 0; iz < (int)m_DivZ; ++iz)
	{
		for (int iy = 0; iy < (int)m_DivY; ++iy)
		{
			for (int ix = 0; ix < (int)m_DivX; ++ix)
			{
				uint cellIdx = (iz * m_DivY + iy) * m_DivX + ix;

				//Note <var>0 corresponds to the cell at current time, and <var>1 corresponds to cell at current time - delta time (old cell from which to advect)
				Vector3 vel0 = GetVelocityMidPoint(ix, iy, iz);

				//float pressure = m_Pressure[(iz * m_DivY + iy) * m_DivX + ix];
				
				Vector3 wsPos0= GetCellWorldSpace(ix, iy, iz);

				Vector3 wsPos1 = wsPos0 - vel0 * sub_timestep;

				Vector3 relCell;
				GetCellLocalSpace(wsPos1, relCell);

				int ix1 = (int)floorf(relCell.x);
				int iy1 = (int)floorf(relCell.y);
				int iz1 = (int)floorf(relCell.z);

#if USE_INTERPOLATED_ADVECTION_CELL==TRUE
				//Pressure
				float alphax = relCell.x - (float)ix1;
				float alphay = relCell.y - (float)iy1;
				float alphaz = relCell.z - (float)iz1;

				Vector3 vel1 = GetVelocityTrilinear(ix1, iy1, iz1, alphax, alphay, alphaz);

				m_IntermediatePressure[cellIdx] = GetPressureTrilinear(ix1, iy1, iz1, alphax, alphay, alphaz);
				m_IntermediateVelocityX[cellIdx] = vel1.x;
				m_IntermediateVelocityY[cellIdx] = vel1.y;
				m_IntermediateVelocityZ[cellIdx] = vel1.z;
#else


				Vector3 vel1 = GetVelocityMidPoint(ix1, iy1, iz1);
				float pressure1 = GetPressure(ix1, iy1, iz1);

				m_IntermediatePressure[cellIdx] = pressure1;
				m_IntermediateVelocityX[cellIdx] = vel1.x;
				m_IntermediateVelocityY[cellIdx] = vel1.y;
				m_IntermediateVelocityZ[cellIdx] = vel1.z;
#endif
			}
		}
	}	
	
	const int padx = 2;
	const int pady = 2;
	const int padz = 0;

	//Set Pressure
#pragma omp parallel for
	for (int iz = padz; iz < (int)m_DivZ - padz; ++iz)
	{
		for (int iy = pady; iy < (int)m_DivY - pady; ++iy)
		{
			for (int ix = padx; ix < (int)m_DivX - padx; ++ix)
			{
				uint cellIdx = (iz * m_DivY + iy) * m_DivX + ix;

				m_Pressure[cellIdx] = m_IntermediatePressure[cellIdx];
			}
		}
	}

/*
	//Set Velocity X
#pragma omp parallel for
	for (int iz = padz; iz < (int)m_DivZ- padz; ++iz)
	{
		for (int iy = pady; iy < (int)m_DivY - pady; ++iy)
		{
			for (int ix = padx; ix < (int)m_DivX+1 - padx; ++ix)
			{
				float vel0 = CellExists(ix - 1, iy, iz) ? m_IntermediateVelocityX[(iz * m_DivY + iy) * m_DivX + ix - 1] : m_ExternalVelocityX;
				float vel1 = CellExists(ix, iy, iz) ? m_IntermediateVelocityX[(iz * m_DivY + iy) * m_DivX + ix] : m_ExternalVelocityX;

				SetVelocityX(ix, iy, iz, (vel0 + vel1) * 0.5f);
			}
		}
	}

	//Set Velocity Y
#pragma omp parallel for
	for (int iz = padz; iz < (int)m_DivZ - padz; ++iz)
	{
		for (int iy = pady; iy < (int)m_DivY + 1 - pady; ++iy)
		{
			for (int ix = padx; ix < (int)m_DivX - padx; ++ix)
			{
				float vel0 = CellExists(ix, iy - 1, iz) ? m_IntermediateVelocityY[(iz * m_DivY + iy-1) * m_DivX + ix] : m_ExternalVelocityY;
				float vel1 = CellExists(ix, iy, iz) ? m_IntermediateVelocityY[(iz * m_DivY + iy) * m_DivX + ix] : m_ExternalVelocityY;

				SetVelocityY(ix, iy, iz, (vel0 + vel1) * 0.5f);
			}
		}
	}

	//Set Velocity Z
#pragma omp parallel for
	for (int iz = padz; iz < (int)m_DivZ + 1 - padz; ++iz)
	{
		for (int iy = pady; iy < (int)m_DivY - pady; ++iy)
		{
			for (int ix = padx; ix < (int)m_DivX - padx; ++ix)
			{
				float vel0 = CellExists(ix, iy, iz - 1) ? m_IntermediateVelocityZ[((iz - 1) * m_DivY + iy) * m_DivX + ix] : m_ExternalVelocityZ;
				float vel1 = CellExists(ix, iy, iz) ? m_IntermediateVelocityZ[(iz * m_DivY + iy) * m_DivX + ix] : m_ExternalVelocityZ;

				SetVelocityZ(ix, iy, iz, (vel0 + vel1) * 0.5f);
			}
		}
	}*/

}