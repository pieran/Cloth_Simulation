#pragma once

#include <ncltech\NCLDebug.h>
#include "ClothDesignEntity.h"

typedef unsigned int uint;

#include <nclgl\Matrix3.h>
#include <map>
#include <mutex>

#define GRAVITY Vector3(0, -9.81f, 0)
#define PBDSOLVER_ITERATIONS 20

struct PBDDistanceConstraint { int p1, p2;	float rest_length, k; float k_prime; };
struct PBDBendingConstraint { int p1, p2, p3; float rest_length, w, k; float k_prime; };

struct Ellipsoid
{
	Matrix4 Transform;
	Matrix4 invTransform;
	float radius;
};

class ClothRenderObject;
class ClothScene;

class PositionBasedDynamicsMS
{
	friend class ClothRenderObject;
	friend class ClothScene;

public:
	PositionBasedDynamicsMS();
	~PositionBasedDynamicsMS();

	//Simulation
	void simulation_OnClothDesignChanged(ClothDesignEntity* design);
	void Simulation_StepSimulation(float dt);

	inline void AddCollisionEllipsoid(Ellipsoid& e) { m_colEllipsoids.push_back(e); }

	//Render
	void Render_DrawingToVisualiser();

protected:
	void AddDistanceConstraint(int a, int b, float k);
	void AddBendingConstraint(int pa, int pb, int pc, float k);
	void BuildNormals();
	void UpdateInternalConstraints(float dt);
	void UpdateDistanceConstraints();
	void UpdateBendingConstraints();
	void GroundCollision();
	void EllipsoidCollision(const Ellipsoid& e);

protected:
	unsigned int m_NumWidth, m_NumHeight;
	unsigned int m_NumTotal, m_NumSprings;
	float m_TotalArea;

	unsigned int m_NumTriangles;
	std::vector<Triangle> m_Triangles; //Used only for normal calculation

	//Phyxel Data
	std::vector<Vector3>	m_PhyxelsPos;
	std::vector<Vector2>	m_PhyxelTexCoords;
	std::vector<Vector3>	m_PhyxelsPosNew;		//New position, based on Verlet integration
	std::vector<Vector3>	m_PhyxelsVel;
	std::vector<Vector3>	m_PhyxelForces;
	std::vector<bool>		m_PhyxelIsStatic;
	std::vector<float>		m_PhyxelsMass;
	std::vector<float>		m_PhyxelsInvMass;
	std::vector<Vector3>	m_PhyxelNormals;
	std::vector<Vector3>	m_PhyxelRi;	//Relative offset from centre of mass, for use when applying intertial damping to the cloth

	//Spring Data
	std::vector<PBDDistanceConstraint> d_constraints;
	std::vector<PBDBendingConstraint> b_constraints;

	//Collidable Objects
	std::vector<Ellipsoid> m_colEllipsoids;
};