#pragma once

#include <nclgl\Vector3.h>
#include <nclcore\common.h>
#include <nclcore\PArray.h>
#include <nclgl\Vector2.h>
#include "mpcg.h"

#define USE_INTERPOLATED_ADVECTION_CELL TRUE

struct BoundingBox2D
{
	Vector2 offset_min;
	Vector2 scale;
};

class MacGrid2D
{
	friend class FluidRenderer;
public:
	MacGrid2D(const BoundingBox2D& region, uint div_x, uint div_y);
	~MacGrid2D();

	//inline BoundingBox& BoundingArea() { return m_BoundingRegion; } ERR: Does not update m_CellDimensions!!
	inline const BoundingBox2D& BoundingRegion() const { return m_BoundingRegion; }

	float CalcMaxVelocity();
	inline float GetMinGridDimension() { return min(m_CellDimensions.x, m_CellDimensions.y); }

	void ApplyGravity(float sub_timestep);
	void AdvectFluid(float sub_timestep);
	void ProjectPressure(float sub_timestep);


	void SetCellPressure(float pressure, int start_x, int start_y, int end_x, int end_y)
	{
#pragma omp parallel for
		for (int y = start_y; y <= end_y; ++y)
		{
			for (int x = start_x; x <= end_x; ++x)
			{
				m_Pressure[y * m_DivX + x] = pressure;
			}
		}
	}

	inline bool CellExists(int ix, int iy) { return (ix >= 0 && ix < m_DivX) && (iy >= 0 && iy < m_DivY); }
	inline bool VelocityXExists(int ix, int iy) { return (ix >= 0 && ix < m_DivX + 1) && (iy >= 0 && iy < m_DivY); }
	inline bool VelocityYExists(int ix, int iy) { return (ix >= 0 && ix < m_DivX) && (iy >= 0 && iy < m_DivY + 1); }

	inline float GetPressure(int ix, int iy) { return CellExists(ix, iy) ? m_Pressure[iy * m_DivX + ix] : m_ExternalPressure; }
	inline float GetVelocityX(int ix, int iy) { return VelocityXExists(ix, iy) ? m_VelocityX[iy * (m_DivX + 1) + ix] : m_ExternalVelocityX; }
	inline float GetVelocityY(int ix, int iy) { return VelocityYExists(ix, iy) ? m_VelocityY[iy * m_DivX + ix] : m_ExternalVelocityY; }

	inline void SetVelocityX(int ix, int iy, float vx)
	{
		if (VelocityXExists(ix, iy))
			m_VelocityX[iy * (m_DivX + 1) + ix] = vx;
	}
	inline void SetVelocityY(int ix, int iy, float vy)
	{
		if (VelocityYExists(ix, iy))
			m_VelocityY[iy * m_DivX + ix] = vy;
	}
	inline void SetPressure(int ix, int iy, float p)
	{
		if (CellExists(ix, iy))
			m_Pressure[iy * m_DivX + ix] = p;
	}
		

	inline float GetPressure(const Vector2& ws_pos)
	{
		Vector2 cel_space;
		GetCellLocalSpace(ws_pos, cel_space);

		int ix = (int)floorf(cel_space.x);
		int iy = (int)floorf(cel_space.y);

		float alphax = cel_space.x - (float)ix;
		float alphay = cel_space.y - (float)iy;

		return GetPressureTrilinear(ix, iy, alphax, alphay);
	}

	inline Vector2 GetVelocityTrilinear(const Vector2& ws_pos)
	{
#define USE_FULL_TRILINEAR_VELOCITY TRUE

		Vector2 cel_space;
		GetCellLocalSpace(ws_pos, cel_space);
#if (USE_FULL_TRILINEAR_VELOCITY == FALSE)
		cel_space.x += 0.5f;
		cel_space.y += 0.5f;
#endif
		int ix = (int)floorf(cel_space.x);
		int iy = (int)floorf(cel_space.y);

		float alphax = cel_space.x - (float)ix;
		float alphay = cel_space.y - (float)iy;

#if (USE_FULL_TRILINEAR_VELOCITY == FALSE)
		return GetVelocityTrilinear(ix, iy, alphax, alphay);
		//return GetVelocityMidPoint(ix, iy, iz);
#else
		Vector2 p_blb = GetVelocityMidPoint(ix, iy);
		Vector2 p_brb = GetVelocityMidPoint(ix + 1, iy);
		Vector2 p_blt = GetVelocityMidPoint(ix, iy + 1);
		Vector2 p_brt = GetVelocityMidPoint(ix + 1, iy + 1);

		Vector2 p_bl = p_blb * (1.0f - alphay) + p_blt * alphay;
		Vector2 p_br = p_brb * (1.0f - alphay) + p_brt * alphay;

		return p_bl * (1.0f - alphax) + p_br * alphax;
#endif
	}

	inline Vector2 GetVelocityTrilinear(int ix, int iy, float alphax, float alphay)
	{
		Vector2 out_vel;
		out_vel.x = (1.0f - alphax) * GetVelocityX(ix, iy) + alphax * GetVelocityX(ix + 1, iy);
		out_vel.y = (1.0f - alphay) * GetVelocityY(ix, iy) + alphay * GetVelocityY(ix, iy + 1);
		
		return out_vel;
	}

	inline Vector2 GetVelocityMidPoint(int ix, int iy) //GetVelocityMidPoint(ix, iy, iz, 0.5f, 0.5f, 0.5f)
	{
		Vector2 out_vel;
		out_vel.x = (GetVelocityX(ix, iy) + GetVelocityX(ix + 1, iy)) * 0.5f;
		out_vel.y = (GetVelocityY(ix, iy) + GetVelocityY(ix, iy + 1)) * 0.5f;
		

		return out_vel;
	}

	inline float GetPressureTrilinear(int ix, int iy, float alphax, float alphay)
	{
		float p_blb = GetPressure(ix, iy);
		float p_brb = GetPressure(ix + 1, iy);
		float p_blt = GetPressure(ix, iy + 1);
		float p_brt = GetPressure(ix + 1, iy + 1);

		float p_bl = (1.0f - alphay) * p_blb + alphay * p_blt;
		float p_br = (1.0f - alphay) * p_brb + alphay * p_brt;

		return (1.0f - alphax) * p_bl + alphax * p_br;
	}



	inline Vector2 GetCellWorldSpace(int ix, int iy) //Given from centre of cell (where pressure/heat etc is stored)
	{
		Vector2 out_pos;
		out_pos.x = m_CellDimensions.x * ((float)ix + 0.5f) + m_BoundingRegion.offset_min.x;
		out_pos.y = m_CellDimensions.y * ((float)iy + 0.5f) + m_BoundingRegion.offset_min.y;
		
		return out_pos;
	}

	inline void GetCellLocalSpace(const Vector2& wsPos, Vector2& ls)
	{
		GetCellLocalSpace(wsPos, ls.x, ls.y);
	}

	inline void GetCellLocalSpace(const Vector2& wsPos, float& fx, float& fy)
	{
		fx = (wsPos.x - m_BoundingRegion.offset_min.x) / m_CellDimensions.x - 0.5f;
		fy = (wsPos.y - m_BoundingRegion.offset_min.y) / m_CellDimensions.y - 0.5f;
	}

	inline void GetCellLocalSpaceFloored(const Vector2& wsPos, int& ix, int& iy)
	{
		ix = (int)floorf((wsPos.x - m_BoundingRegion.offset_min.x) / m_CellDimensions.x - 0.5f);
		iy = (int)floorf((wsPos.y - m_BoundingRegion.offset_min.y) / m_CellDimensions.y - 0.5f);
	}


protected:
	//External variables for 'outside' the MAC grid
	float m_ExternalPressure;
	float m_ExternalVelocityX;
	float m_ExternalVelocityY;

	uint m_DivX, m_DivY; //For Centre of Voxel (Velocity components have +1 dimensions on there half axis)
	BoundingBox2D m_BoundingRegion;
	Vector2 m_CellDimensions;

	//MAC GRID
	PArray<float> m_Pressure;
	PArray<float> m_VelocityX;
	PArray<float> m_VelocityY;

	//SUB_MAC_GRID (used to store advected cell data at all centres (intemediate data structure used by advection algorithm)
	PArray<float> m_IntermediatePressure;
	PArray<float> m_IntermediateVelocityX;
	PArray<float> m_IntermediateVelocityY;


	//Solver for pressure
	MPCG m_Solver;
};