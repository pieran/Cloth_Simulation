#pragma once

#include "MacGrid2D.h"

class FluidSim2D
{
	friend class MyScene;
public:
	FluidSim2D(const BoundingBox2D& region, uint div_x, uint div_y);
	~FluidSim2D();

	void InitializeFluid();
	void StepSimulation(float global_timestep);

protected:

	void AdvectFluid(float sub_timestep);
	void SolvePressure(float sub_timestep);
	void AdvectBoundaries(float sub_timestep);


protected:
	float		m_AccumTime;
	MacGrid2D	m_MACGrid;

};