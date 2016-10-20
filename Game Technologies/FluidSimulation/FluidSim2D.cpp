#include "FluidSim2D.h"

FluidSim2D::FluidSim2D(const BoundingBox2D& region, uint div_x, uint div_y)
	: m_MACGrid(region, div_x, div_y)
	, m_AccumTime(0.0f)
{

}

FluidSim2D::~FluidSim2D()
{

}

void FluidSim2D::InitializeFluid()
{

}

void FluidSim2D::StepSimulation(float global_timestep)
{
	m_AccumTime += global_timestep;

	float sub_timestep = min(global_timestep, m_MACGrid.GetMinGridDimension() / max(1.0f, m_MACGrid.CalcMaxVelocity()));
	//Better Solution =>
	// sub_timestep += sqrt(dh * CalcMaxForce())

	while (m_AccumTime >= sub_timestep)
	{
		printf("SubTimeStep: %5.2f\n", sub_timestep);
		//Update Velocity Field
		m_MACGrid.ApplyGravity(sub_timestep);

		//Advect Fluid
		AdvectFluid(sub_timestep);
		
		//Pressure Projection
		SolvePressure(sub_timestep);

		//Advect Boundaries


		m_AccumTime -= sub_timestep;
		sub_timestep = m_MACGrid.CalcMaxVelocity();
	}
}


void FluidSim2D::AdvectFluid(float sub_timestep)
{
	m_MACGrid.AdvectFluid(sub_timestep);
}

void FluidSim2D::SolvePressure(float sub_timestep)
{
	//Build A matrix and B vector in order to solve linear equation Ax = b
	
	//Iterate over each grid point
	m_MACGrid.ProjectPressure(sub_timestep);
}

void FluidSim2D::AdvectBoundaries(float sub_timestep)
{

}