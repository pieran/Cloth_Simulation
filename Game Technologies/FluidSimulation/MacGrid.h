#pragma once

#include <nclgl\Vector3.h>
#include <nclcore\common.h>
#include <nclcore\PArray.h>

#define USE_INTERPOLATED_ADVECTION_CELL TRUE

struct BoundingBox
{
	Vector3 offset_min;
	Vector3 scale;
};

class MacGrid
{
	friend class FluidRenderer;
public:
	MacGrid(const BoundingBox& region, uint div_x, uint div_y, uint div_z);
	~MacGrid();

	//inline BoundingBox& BoundingArea() { return m_BoundingRegion; } ERR: Does not update m_CellDimensions!!
	inline const BoundingBox& BoundingRegion() const { return m_BoundingRegion; }
	
	float CalcMaxVelocity();
	inline float MacGrid::GetMinGridDimension() { return min(m_CellDimensions.x, min(m_CellDimensions.y, m_CellDimensions.z)); }

	void ApplyGravity(float sub_timestep);
	void AdvectFluid(float sub_timestep);

	void SetCellPressure(float pressure, int start_x, int start_y, int start_z, int end_x, int end_y, int end_z)
	{
#pragma omp parallel for
		for (int z = start_z; z <= end_z; ++z)
		{
			for (int y = start_y; y <= end_y; ++y)
			{
				for (int x = start_x; x <= end_x; ++x)
				{
					m_Pressure[(z * m_DivY + y) * m_DivX + x] = pressure;
				}
			}
		}
	}

	inline bool CellExists(int ix, int iy, int iz) { return (ix >= 0 && ix < m_DivX) && (iy >= 0 && iy < m_DivY) && (iz >= 0 && iz < m_DivZ); }
	inline bool VelocityXExists(int ix, int iy, int iz) { return (ix >= 0 && ix < m_DivX + 1) && (iy >= 0 && iy < m_DivY) && (iz >= 0 && iz < m_DivZ); }
	inline bool VelocityYExists(int ix, int iy, int iz) { return (ix >= 0 && ix < m_DivX) && (iy >= 0 && iy < m_DivY + 1) && (iz >= 0 && iz < m_DivZ); }
	inline bool VelocityZExists(int ix, int iy, int iz) { return (ix >= 0 && ix < m_DivX) && (iy >= 0 && iy < m_DivY) && (iz >= 0 && iz < m_DivZ + 1); }

	inline float GetPressure(int ix, int iy, int iz) { return CellExists(ix, iy, iz) ? m_Pressure[(iz * m_DivY + iy) * m_DivX + ix] : m_ExternalPressure; }
	inline float GetVelocityX(int ix, int iy, int iz) { return VelocityXExists(ix, iy, iz) ? m_VelocityX[(iz * m_DivY + iy) * (m_DivX + 1) + ix] : m_ExternalVelocityX; }
	inline float GetVelocityY(int ix, int iy, int iz) { return VelocityYExists(ix, iy, iz) ? m_VelocityY[(iz * (m_DivY+1) + iy) * m_DivX + ix] : m_ExternalVelocityY; }
	inline float GetVelocityZ(int ix, int iy, int iz) { return VelocityZExists(ix, iy, iz) ? m_VelocityZ[(iz * m_DivY + iy) * m_DivX + ix] : m_ExternalVelocityZ; }

	inline void SetVelocityX(int ix, int iy, int iz, float vx)
	{
		if (VelocityXExists(ix, iy, iz)) 
			m_VelocityX[(iz * m_DivY + iy) * (m_DivX + 1) + ix] = vx;
	}
	inline void SetVelocityY(int ix, int iy, int iz, float vy)
	{
		if (VelocityYExists(ix, iy, iz))
			m_VelocityY[(iz * (m_DivY + 1) + iy) * m_DivX + ix] = vy;
	}
	inline void SetVelocityZ(int ix, int iy, int iz, float vz)
	{
		if (VelocityZExists(ix, iy, iz))
			m_VelocityZ[(iz * m_DivY + iy) * m_DivX + ix] = vz;
	}

	inline float GetPressure(const Vector3& ws_pos)
	{
		int ix, iy, iz;
		GetCellLocalSpaceFloored(ws_pos, ix, iy, iz);

		return GetPressure(ix, iy, iz);
		/*return Interpolate(m_Pressure,
		(int)floorf(celx), (int)floorf(cely), (int)floorf(celz),
		celx - floorf(celx),
		cely - floorf(cely),
		celz - floorf(celz));*/
	}

	inline Vector3 GetVelocityTrilinear(const Vector3& ws_pos)
	{
		Vector3 cel_space;
		GetCellLocalSpace(ws_pos, cel_space);
		cel_space = cel_space - 0.5f;
		int ix = (int)floorf(cel_space.x);
		int iy = (int)floorf(cel_space.y);
		int iz = (int)floorf(cel_space.z);

		float alphax = cel_space.x - (float)ix;
		float alphay = cel_space.y - (float)iy;
		float alphaz = cel_space.z - (float)iz;

#if 1
		return GetVelocityTrilinear(ix, iy, iz, alphax, alphay, alphaz);
		//return GetVelocityMidPoint(ix, iy, iz);
#else


		Vector3 p_blb = GetVelocityMidPoint(ix, iy, iz);
		Vector3 p_brb = GetVelocityMidPoint(ix + 1, iy, iz);
		Vector3 p_blt = GetVelocityMidPoint(ix, iy + 1, iz);
		Vector3 p_brt = GetVelocityMidPoint(ix + 1, iy + 1, iz);

		Vector3 p_flb = GetVelocityMidPoint(ix, iy, iz + 1);
		Vector3 p_frb = GetVelocityMidPoint(ix + 1, iy, iz + 1);
		Vector3 p_flt = GetVelocityMidPoint(ix, iy + 1, iz + 1);
		Vector3 p_frt = GetVelocityMidPoint(ix + 1, iy + 1, iz + 1);

		Vector3 p_bl = p_blb * (1.0f - alphay) + p_blt * alphay;
		Vector3 p_br = p_brb * (1.0f - alphay) + p_brt * alphay;
		Vector3 p_fl = p_flb * (1.0f - alphay) + p_flt * alphay;
		Vector3 p_fr = p_frb * (1.0f - alphay) + p_frt * alphay;

		Vector3 p_b = p_bl * (1.0f - alphax) + p_br * alphax;
		Vector3 p_f = p_fl * (1.0f - alphax) + p_fr * alphax;

		return p_b * (1.0f - alphaz) + p_f * alphaz;
#endif
	}

	inline Vector3 GetVelocityTrilinear(int ix, int iy, int iz, float alphax, float alphay, float alphaz)
	{
		Vector3 out_vel;
		out_vel.x = (1.0f - alphax) * GetVelocityX(ix, iy, iz) + alphax * GetVelocityX(ix + 1, iy, iz);
		out_vel.y = (1.0f - alphay) * GetVelocityY(ix, iy, iz) + alphay * GetVelocityY(ix, iy + 1, iz);
		out_vel.z = (1.0f - alphaz) * GetVelocityZ(ix, iy, iz) + alphaz * GetVelocityZ(ix, iy, iz + 1);
		return out_vel;
	}

	inline Vector3 GetVelocityMidPoint(int ix, int iy, int iz) //GetVelocityMidPoint(ix, iy, iz, 0.5f, 0.5f, 0.5f)
	{
		Vector3 out_vel;
		out_vel.x = (GetVelocityX(ix, iy, iz) + GetVelocityX(ix + 1, iy, iz)) * 0.5f;
		out_vel.y = (GetVelocityY(ix, iy, iz) + GetVelocityY(ix, iy + 1, iz)) * 0.5f;
		out_vel.z = (GetVelocityZ(ix, iy, iz) + GetVelocityZ(ix, iy, iz + 1)) * 0.5f;

		return out_vel;
	}
	
	inline float GetPressureTrilinear(int ix, int iy, int iz, float alphax, float alphay, float alphaz)
	{
		float p_blb = GetPressure(ix, iy, iz);
		float p_brb = GetPressure(ix + 1, iy, iz);
		float p_blt = GetPressure(ix, iy + 1, iz);
		float p_brt = GetPressure(ix + 1, iy + 1, iz);

		float p_flb = GetPressure(ix, iy, iz + 1);
		float p_frb = GetPressure(ix + 1, iy, iz + 1);
		float p_flt = GetPressure(ix, iy + 1, iz + 1);
		float p_frt = GetPressure(ix + 1, iy + 1, iz + 1);

		float p_bl = (1.0f - alphay) * p_blb + alphay * p_blt;
		float p_br = (1.0f - alphay) * p_brb + alphay * p_brt;
		float p_fl = (1.0f - alphay) * p_flb + alphay * p_flt;
		float p_fr = (1.0f - alphay) * p_frb + alphay * p_frt;

		float p_b = (1.0f - alphax) * p_bl + alphax * p_br;
		float p_f = (1.0f - alphax) * p_fl + alphax * p_fr;

		return (1.0f - alphaz) * p_b + alphaz * p_f;
	}



	inline Vector3 GetCellWorldSpace(int ix, int iy, int iz) //Given from centre of cell (where pressure/heat etc is stored)
	{
		Vector3 out_pos;
		out_pos.x = m_CellDimensions.x * ((float)ix + 0.5f) + m_BoundingRegion.offset_min.x;
		out_pos.y = m_CellDimensions.y * ((float)iy + 0.5f) + m_BoundingRegion.offset_min.y;
		out_pos.z = m_CellDimensions.z * ((float)iz + 0.5f) + m_BoundingRegion.offset_min.z;
		return out_pos;
	}

	inline void GetCellLocalSpace(const Vector3& wsPos, Vector3& ls)
	{
		GetCellLocalSpace(wsPos, ls.x, ls.y, ls.z);
	}

	inline void GetCellLocalSpace(const Vector3& wsPos, float& fx, float& fy, float& fz)
	{
		fx = (wsPos.x - m_BoundingRegion.offset_min.x) / m_CellDimensions.x - 0.5f;
		fy = (wsPos.y - m_BoundingRegion.offset_min.y) / m_CellDimensions.y - 0.5f;
		fz = (wsPos.z - m_BoundingRegion.offset_min.z) / m_CellDimensions.z - 0.5f;
	}

	inline void GetCellLocalSpaceFloored(const Vector3& wsPos, int& ix, int& iy, int& iz)
	{
		ix = (int)floorf((wsPos.x - m_BoundingRegion.offset_min.x) / m_CellDimensions.x - 0.5f);
		iy = (int)floorf((wsPos.y - m_BoundingRegion.offset_min.y) / m_CellDimensions.y - 0.5f);
		iz = (int)floorf((wsPos.z - m_BoundingRegion.offset_min.z) / m_CellDimensions.z - 0.5f);
	}


protected:
	//External variables for 'outside' the MAC grid
	float m_ExternalPressure;
	float m_ExternalVelocityX;
	float m_ExternalVelocityY;
	float m_ExternalVelocityZ;


	int m_DivX, m_DivY, m_DivZ; //For Centre of Voxel (Velocity components have +1 dimensions on there half axis)
	BoundingBox m_BoundingRegion;
	Vector3 m_CellDimensions;

	//MAC GRID
	PArray<float> m_Pressure;
	PArray<float> m_VelocityX;
	PArray<float> m_VelocityY;
	PArray<float> m_VelocityZ;

	//SUB_MAC_GRID (used to store advected cell data at all centres (intemediate data structure used by advection algorithm)
	PArray<float> m_IntermediatePressure;
	PArray<float> m_IntermediateVelocityX;
	PArray<float> m_IntermediateVelocityY;
	PArray<float> m_IntermediateVelocityZ;

};