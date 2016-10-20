#pragma once

#include "PArray.h"
#include "SimulationProfiler.h"
#include "ClothDesignEntity.h"

#include <nclgl\Vector3.h>
#include <nclgl\Vector2.h>

#include <functional>


//Maximum number of user-generated new points TODO: Make dynamic without large memory reallocation hit
#define MAX_ADDITIONAL_PHYXELS 1000
#define MAX_ADDITIONAL_TRIANGLES 500


//Set to false to directly calculate normals from individual triangles (probably best for cutting purposes.. but im not sure)
#define USE_QUAD_NORMALS TRUE


template <class T>
T barycentric_interpolate(const T& a, const T& b, const T& c, const Vector3& b_coords)
{
	return a * b_coords.x
		+ b * b_coords.y
		+ c * b_coords.z;
}

struct ClothState
{
	uint numPhyxels;
	Vector3* phyxelPos;
	Vector3* phyxelVel;	
};

enum ClothIntegrationScheme
{
	CLOTH_INTEGRATION_EXPLICIT,
	CLOTH_INTEGRATION_LEAPFROG,
	CLOTH_INTEGRATION_RK4
};

struct Ray //For Cutting Ray-Intersection Tests
{
	Vector3 wsPos;
	Vector3 wsDir;
};

class ClothBase
{
	friend class ClothScenarioViewer;
	friend class MyScene;
public:
	ClothBase(ClothIntegrationScheme integrator = CLOTH_INTEGRATION_EXPLICIT);
	virtual ~ClothBase();


	virtual void InitializeClothDesign(const ClothDesignEntity<3>* design);
	void UpdateSimulation(float global_timestep);

	//Triangle index given as exact offset in 'm_TriIndicies' to the start of the triangle.
	//Barycentric coordinates given as per wikipedia
	//returns index of new point inside phyxel array
	//throws OUT_OF_MEMORY exception if there is no more memory available for new points to be added (WILL FIX LATER)
	virtual uint AddPoint(const uint& triIdx, const Vector3& barycentric_coords);

	virtual uint AddTriangle(const uint& idxA, const uint& idxB, const uint& idxC);
	virtual void UpdateTriangleMaterial(const uint& triIdx); //MUST be called after changing any triangle indicies or creating a new triangle, this allows the springs/fem to readjust the properties to match the new triangle.

	//GETTERS & SETTERS
	//------------------
	bool GetClothState(uint idx, ClothState* out_clothstate);

	inline void SetCallback_OnOuterUpdateScene(const std::function<void(float)>& callback) { m_Callback_OnOuterUpdateScene = callback; }
	inline void SetCallback_OnOuterCollisions(const std::function<void()>& callback) { m_Callback_OnOuterCollisions = callback; }

	inline void SetCallback_OnInnerFullTimestep(const std::function<void(float)>& callback) { m_Callback_OnInnerFullTimestep = callback; }
	inline void SetCallback_OnInnerApplySubTimestep(const std::function<void(float, uint, Vector3*, Vector3*)>& callback) { m_Callback_OnInnerApplySubTimestep = callback; }
	inline void SetCallback_OnInnerCollisions(const std::function<void()>& callback) { m_Callback_OnInnerCollisions = callback; }


	void SetIntegrator(ClothIntegrationScheme integrator);
	inline ClothIntegrationScheme GetIntegrator() const { return m_ClothIntegrationScheme; }

	void AdvanceVelocityHalfTimestep();
	bool GetClosestPointsForCutting(const Ray& start, const Ray& end, int* out_idx_start, int* out_idx_end);


	inline SimulationProfiler& Profiler() { return m_Profiler; }
	inline const SimulationProfiler& Profiler() const { return m_Profiler; }

	inline float& Timestep() { return m_TimeStep; }
	inline const float& Timestep() const { return m_TimeStep; }

	inline bool IsPhyxelStatic(uint idx) const { return m_PhyxelIsStatic[idx]; }



protected:
	virtual void ComputeVelocity(float subtimestep, const Vector3* in_clothpos, const Vector3* in_clothvel, Vector3* out_clothvel) = 0;

	void StepSimulationExplicit(float timestep, const Vector3* in_clothstatepos, const Vector3* in_clothstatevel, Vector3* out_clothstatepos, Vector3* out_clothstatevel);
	void StepSimulationLeapFrog();
	void StepSimulationRK4();

	void AllocateDataStructures(uint numPhyxels);


	
	void GenNormals(const Vector3* in_clothpos);
	void ApplyAirResistance(const Vector3* in_clothpos, const Vector3* in_clothvel);

public:
	SimulationProfiler m_Profiler;

	ClothIntegrationScheme m_ClothIntegrationScheme;

	//CLOTH STATE
	//----------------------
	//	INTEGRATION SPECIFIC
	uint	m_NumPhyxels;
	uint	m_NumAllocatedPhyxels;	//Number of possible phyxels!
	uint	m_NumTris, m_NumAllocatedTris, m_NumQuads;


	uint	m_NumClothStates; //Used by integrator schemes to store sub-steps
	PArray<Vector3> m_ClothStatesPos;
	PArray<Vector3> m_ClothStatesVel;

	//	GENERIC
	PArray<Vector3> m_PhyxelForces;					//Internal forces generated by simulation (E.G. Gravity)
	PArray<Vector3> m_PhyxelExtForces;				//External Forces from external collisions/constraints
	PArray<Vector3> m_PhyxelNormals;				//Normals of each phyxel
	PArray<bool>	m_PhyxelIsStatic;				//Is the phyxel static? (Purely for testing of constraints and material properties)
	PArray<float>	m_PhyxelMass;					//Mass of given phyxel (Currently uniform but planned to make more generic)



	//VISUAL OUTPUT
	//-----------------------
	PArray<Vector2> m_PhyxelTexCoords;				//Purely Visual!
	PArray<float>	m_RenderValue;					//Output render values (TODO: MAKE MORE GENERIC!)
	Vector2			m_RenderValueMinMax;


	//POINT->ELEMENT DESCRIPTORS
	//-----------------------
	PArray<uint>	m_TriIndicies;
	PArray<uint>	m_QuadIndicies;					//Allows for a slightly better normal generation/interpolation method


	//CALL-BACKS
	//-----------------------
	std::function<void(float)> m_Callback_OnOuterUpdateScene;
	std::function<void(void)>  m_Callback_OnOuterCollisions;
	
	std::function<void(float)>  m_Callback_OnInnerFullTimestep;
	std::function<void(float, uint, Vector3*, Vector3*)>  m_Callback_OnInnerApplySubTimestep;
	std::function<void(void)>   m_Callback_OnInnerCollisions;

private:

	float	m_TimeStep;
	float	m_InitialArea;
	float	m_AccumTime;
};