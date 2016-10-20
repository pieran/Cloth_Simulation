#pragma once

#include "ClothScenario.h"
#include "FEDefault.h"
#include "ClothDesignEntity.h"
#include <ncltech\SimpleMeshObject.h>
#include <nclgl\OBJMesh.h>
#include <nclgl\Vector3.h>
#include <nclgl\Vector2.h>

#define WEIGHTED_MODE 
#define WEIGHT_RADIUS 0.02f

class ScenarioCollisionYusef : public ClothScenario
{
public:
	ScenarioCollisionYusef(ClothBase* clothsim, int numDivsX = 10, int numDivsY = 10, const Vector2& scale = Vector2(1.0f, 1.0f), const Vector3& pos = Vector3(0, 1.0f, -5.0f))
		: ClothScenario(clothsim), m_ClothDivsX(numDivsX), m_ClothDivsY(numDivsY), m_ClothScale(scale), m_ClothPosition(pos), m_AccumTime(0.0f)
	{
		using namespace std::placeholders;
		clothsim->SetCallback_OnInnerFullTimestep(std::bind(&ScenarioCollisionYusef::OnInnerFullTimestep, this, _1));
		clothsim->SetCallback_OnInnerApplySubTimestep(std::bind(&ScenarioCollisionYusef::OnInnerApplySubTimestep, this, _1, _2, _3, _4));
		clothsim->SetCallback_OnInnerCollisions(std::bind(&ScenarioCollisionYusef::OnInnerCollisions, this));
		m_AccumTime = 0.0f;

		m_ColRadius = 0.205f;
		m_ColCentre = pos - Vector3(0.0f, 0.15f + m_ColRadius, 0.1f);
		


		SimpleMeshObject* obj = new SimpleMeshObject("ColBall");
		obj->SetMesh(new OBJMesh(MESHDIR"sphere.obj"), true);
		obj->SetLocalTransform(Matrix4::Translation(m_ColCentre) * Matrix4::Scale(Vector3(m_ColRadius-0.01f, m_ColRadius - 0.01f, m_ColRadius - 0.01f)));
		obj->SetBoundingRadius(1.2f);
		obj->SetColour(Vector4(1.0f, 0.0f, 0.0f, 1.0f));
		
		m_WhiteTex = SOIL_load_OGL_texture(TEXTUREDIR"white.bmp", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);
		glBindTexture(GL_TEXTURE_2D, m_WhiteTex);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST); //No linear interpolation to get crisp checkerboard no matter the scalling
		glBindTexture(GL_TEXTURE_2D, 0);
		obj->SetTexture(m_WhiteTex, true);

		this->AddChildObject(obj);



		m_PointA = new SimpleMeshObject("ColBall");
		m_PointA->SetMesh(new OBJMesh(MESHDIR"sphere.obj"), true);
		m_PointA->SetBoundingRadius(1.2f);
		m_PointA->SetColour(Vector4(1.0f, 0.0f, 1.0f, 1.0f));
		m_PointA->SetTexture(m_WhiteTex, false);
		m_PointA->SetLocalTransform(Matrix4::Scale(Vector3(0.0f, 0.0f, 0.0f)));
		this->AddChildObject(m_PointA);

		m_PointB = new SimpleMeshObject("ColBall");
		m_PointB->SetMesh(new OBJMesh(MESHDIR"sphere.obj"), true);
		m_PointB->SetBoundingRadius(1.2f);
		m_PointB->SetColour(Vector4(1.0f, 0.0f, 1.0f, 1.0f));
		m_PointB->SetTexture(m_WhiteTex, false);
		m_PointB->SetLocalTransform(Matrix4::Scale(Vector3(0.0f, 0.0f, 0.0f)));
		this->AddChildObject(m_PointB);

	}

	virtual ~ScenarioCollisionYusef() {}

	//Cleanup any non-cloth data structures
	virtual void OnDeleteScene() override
	{

	}

	//Set All data in cloth class and scenario class
	virtual void OnInitializeScene() override
	{
		m_AccumTime = 0.5f;


		Matrix4 transform = Matrix4::Translation(m_ClothPosition)*Matrix4::Rotation(90.0f, Vector3(1, 0, 0));
		Mat44 tmp;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				tmp(i, j) = transform[j * 4 + i];
			}
		}


		//Create vertical cloth with fixed points at the top and bottom rows
		ClothDesignEntity<3> design;
		design.GenerateDefaultClothVertices(m_ClothDivsX, m_ClothScale.x, m_ClothScale.y, &tmp);
		design.GenerateDefaultClothTriangles(m_ClothDivsX);

		//Set Static Vertices
		for (Vertex& v : design.Vertices())
		{
			v.flags = 0;
		}
		//design.Vertices()[0].flags = VERTEXFLAGS_IS_STATIC;
		//design.Vertices()[m_ClothDivsX - 1].flags = VERTEXFLAGS_IS_STATIC;

		design.Vertices()[(m_ClothDivsY - 1) * m_ClothDivsX + 0].flags = VERTEXFLAGS_IS_STATIC;
		design.Vertices()[(m_ClothDivsY - 1) * m_ClothDivsX + m_ClothDivsX - 1].flags = VERTEXFLAGS_IS_STATIC;

		m_ClothSim->InitializeClothDesign(&design);

		auto add_mass = [&](uint idx, float mass)
		{
			m_ClothSim->m_PhyxelMass[idx] += mass;
			m_ClothSim->m_PhyxelForces[idx] = GRAVITY * m_ClothSim->m_PhyxelMass[idx];
			m_WeightedVerts.push_back(idx);

			SimpleMeshObject* obj = new SimpleMeshObject("ColBall");
			obj->SetMesh(new OBJMesh(MESHDIR"sphere.obj"), true);
			obj->SetBoundingRadius(1.2f);
			obj->SetColour(Vector4(0.5f, 0.4f, 0.1f, 0.5f));
			obj->SetTexture(m_WhiteTex, false);
			obj->SetLocalTransform(Matrix4::Scale(Vector3(WEIGHT_RADIUS, WEIGHT_RADIUS, WEIGHT_RADIUS)));
			obj->Physics()->SetPosition(m_ClothSim->m_ClothStatesPos[idx]);
			this->AddChildObject(obj);

			m_WeightedSpheres.push_back(obj);
		};
		
		//add_mass(0, 1.00f);
		//add_mass(m_ClothDivsX - 1, 1.00f);
		//add_mass((m_ClothDivsY - 1) * m_ClothDivsX + 0, 1.00f);
		//add_mass((m_ClothDivsY - 1) * m_ClothDivsX + m_ClothDivsX - 1, 1.00f);

		PBD* pbdsim = (PBD*)dynamic_cast<PBD*>(m_ClothSim);
		if (pbdsim)
		{
			pbdsim->InitializePBDConstraints();
		}


		XPBD* xpbd = dynamic_cast<XPBD*>(m_ClothSim);
		if (xpbd != NULL)
		{
			//xpbd->InitializePBDConstraints();
			XPBDSphereConstraint sc;
			sc.centre = m_ColCentre;
			sc.radius = m_ColRadius;
			xpbd->m_SphereConstraints.push_back(sc);
		}

	}

	//Called once per cloth update
	void OnOuterUpdateScene(float cloth_timestep)
	{

	}

	//Called once per timestep (possibly multiple cloth steps to each scene timestep) NOTE: Sum of sub_timesteps each frame will equal cloth_timestep so no need to update on both callbacks!
	void OnInnerFullTimestep(float timestep)
	{
		m_AccumTime += timestep;
	}

	void OnInnerApplySubTimestep(float sub_timestep, uint numPhyxels, Vector3* clothPos, Vector3* clothVel)
	{

		/*for (uint i = 0; i < numPhyxels; ++i)
		{
			if (m_ClothSim->IsPhyxelStatic(i))
			{
				float spd = 3.85f;

				clothVel[i].y = sin((m_AccumTime + sub_timestep) * 4.0f) * spd;
			}
		}*/
	}

	//Called once per timestep (possibly multiple cloth steps to each scene timestep)
	void OnOuterCollisions()
	{

	}


	//Called once per cloth update
	void OnInnerCollisions()
	{
		XPBD* xpbd = dynamic_cast<XPBD*>(m_ClothSim);
		if (xpbd != NULL)
		{
			return;
		}

		//Iterate through each Phyxel and handle collisions with the ball

		const float radiusSq = m_ColRadius * m_ColRadius;

		FEBase* fesim = (FEBase*)dynamic_cast<FEBase*>(m_ClothSim);


#pragma omp parallel for
		for (int i = 0; i < m_ClothSim->m_NumPhyxels; ++i)
		{
			Vector3 pos = m_ClothSim->m_ClothStatesPos[i];

			Vector3 axis = pos - m_ColCentre;
			float distSquared = Vector3::Dot(axis, axis);
			
			if (distSquared < radiusSq)
			{
				float dist = sqrt(distSquared);

				float excess = 1.0f - dist / m_ColRadius;
				pos += axis * excess;

				m_ClothSim->m_ClothStatesPos[i] = pos;

				Vector3 axisNrm = axis;
				axisNrm.Normalise();

				Vector3 vel = m_ClothSim->m_ClothStatesVel[i];



				// If the objects are moving away from each other we don't need to apply an impulse
				float relativeMovement = -Vector3::Dot(vel, axisNrm);
				//if (relativeMovement > 0.01f) {

					/*float e = 0.1f; //Bouncy-ness
					float jn = -1 * (1 + e) * Vector3::Dot(vel, axisNrm);


					Vector3 tangent = vel - (axisNrm * Vector3::Dot(vel, axisNrm)); tangent.Normalise();
					float jt = -0.5f * Vector3::Dot(vel, tangent);

					vel += axisNrm * jn + tangent * jt;
					m_ClothSim->m_ClothStatesVel[i] = vel;*/
				//}
			




					if (fesim)
					{
						if (relativeMovement > 0.0001f)
							fesim->SetSolverConstraint(i, Matrix3::ZeroMatrix);// (Matrix3::Identity - Matrix3::OuterProduct(axisNrm, axisNrm)));
						else
							fesim->SetSolverConstraint(i, Matrix3::Identity);
					}
					else
					{
						//m_ClothSim->m_PhyxelIsStatic[i] = true;
					}
		
				
				m_ClothSim->m_ClothStatesVel[i] = Vector3(0.0f, 0.0f, 0.0f);

				//m_ClothSim->m_ClothStatesVel[i] = Vector3(0.0f, 0.0f, 0.0f);
				
			}
			else if (fesim)
			{
				fesim->SetSolverConstraint(i, Matrix3::Identity);
			}
			else
			{
				//m_ClothSim->m_PhyxelIsStatic[i] = false;

				
			}

		}

		const float max_sq_dist = WEIGHT_RADIUS * WEIGHT_RADIUS;
		if (m_WeightedSpheres.size() > 0)
		{
			for (size_t i = 0; i < m_WeightedSpheres.size() - 1; ++i)
			{
				uint idxI = m_WeightedVerts[i];
				Vector3 posI = m_ClothSim->m_ClothStatesPos[idxI];

				for (size_t j = i + 1; j < m_WeightedSpheres.size(); ++j)
				{
					uint idxJ = m_WeightedVerts[j];
					Vector3 posJ = m_ClothSim->m_ClothStatesPos[idxJ];

					Vector3 axis = posI - posJ;
					float distSquared = Vector3::Dot(axis, axis);
					if (distSquared < max_sq_dist)
					{
						float dist = sqrt(distSquared);

						float excess = 1.0f - dist / WEIGHT_RADIUS;
						posI += (axis * excess) * 0.5f;
						posJ -= (axis * excess) * 0.5f;

						m_ClothSim->m_ClothStatesPos[idxI] = posI;
						m_ClothSim->m_ClothStatesPos[idxJ] = posJ;

						m_ClothSim->m_ClothStatesVel[idxI] = Vector3(0.0f, 0.0f, 0.0f);
						m_ClothSim->m_ClothStatesVel[idxJ] = Vector3(0.0f, 0.0f, 0.0f);

						Vector3 axisNrm = axis;
						axisNrm.Normalise();

						Vector3 vel1 = m_ClothSim->m_ClothStatesVel[idxI];
						Vector3 vel2 = m_ClothSim->m_ClothStatesVel[idxJ];

						// If the objects are moving away from each other we don't need to apply an impulse
						Vector3 rv = vel1 - vel2;
						float relativeMovement = -Vector3::Dot(rv, axisNrm);
						if (relativeMovement > 0.01f) {
							float e = 0.5f; //Bouncy-ness
							float jn = -1 * (1 + e) * Vector3::Dot(rv, axisNrm);

							Vector3 tangent = rv - (axisNrm * Vector3::Dot(rv, axisNrm)); tangent.Normalise();
							float jt = -1.0f * Vector3::Dot(rv, tangent);

							vel1 += axisNrm * jn + tangent * jt;
							vel2 -= axisNrm * jn + tangent * jt;
							m_ClothSim->m_ClothStatesVel[idxI] = vel1;
							m_ClothSim->m_ClothStatesVel[idxJ] = vel2;
						}
					}
				}
			}
		
			for (size_t i = 0; i < m_WeightedSpheres.size(); ++i)
			{
				m_WeightedSpheres[i]->Physics()->SetPosition(m_ClothSim->m_ClothStatesPos[m_WeightedVerts[i]]);
			}
		}
	}


	virtual void OnClothCut(int start_idx, int end_idx, const Ray& start_ray, const Ray& end_ray, const Plane& cut_plane) override
	{
		Vector3 posA = m_ClothSim->m_ClothStatesPos[start_idx];
		Vector3 posB = m_ClothSim->m_ClothStatesPos[end_idx];

		const Vector3 color = Vector3(1.0f, 0.7f, 0.85f);
		NCLDebug::Log(color, "Measurement A-B");
		NCLDebug::Log(color, "------------------------------");
		NCLDebug::Log(color, "    Point A: (%5.2f, %5.2f, %5.2f)", posA.x, posA.y, posA.z);
		NCLDebug::Log(color, "    Point B: (%5.2f, %5.2f, %5.2f)", posB.x, posB.y, posB.z);
		NCLDebug::Log(color, "    Distance: %fm", (posB-posA).Length());
		NCLDebug::Log(color, "");

		m_PointA->SetLocalTransform(Matrix4::Translation(posA) * Matrix4::Scale(Vector3(0.01f, 0.01f, 0.01f)));
		m_PointB->SetLocalTransform(Matrix4::Translation(posB) * Matrix4::Scale(Vector3(0.01f, 0.01f, 0.01f)));
	}


protected:
	SimpleMeshObject* m_PointA, *m_PointB;

	std::vector<uint> m_WeightedVerts;
	std::vector<SimpleMeshObject*> m_WeightedSpheres;

	GLuint m_WhiteTex;

	int		m_ClothDivsX, m_ClothDivsY;
	Vector2 m_ClothScale;
	Vector3 m_ClothPosition;
	float	m_AccumTime;
	Vector3 m_ColCentre;
	float	m_ColRadius;
};
