#pragma once
#include <ncltech\GameObject.h>
#include <nclgl\Shader.h>
#include "SPHCommon.h"

class FluidScreenSpaceRenderer : public GameObject
{
	friend class MyScene;
public:
	FluidScreenSpaceRenderer(const std::string& name = "");
	virtual ~FluidScreenSpaceRenderer();

	void SetBackgroundSceneCubemap(GLuint texIdx) { m_BackgroundSceneCubeMap = texIdx; }
	void SetWorldSpace(const Vector3& world_origin, const Vector3& sim_ratio);
	void BuildBuffers(int max_particles);
	void UpdateBuffers(int num_particles, Particle* particles);
	void UpdateBuffers(int num_particles, float4* particles);
protected:
	virtual void OnRenderObject();				//Handles OpenGL calls to Render the object
	virtual void OnUpdateObject(float dt);		//Override to handle things like AI etc on update loop
	virtual void OnAttachedToScene();			//Called when attached to a scene object

	void BuildFBO();
	void ReloadShaders();
protected:
	uint m_MaxParticles, m_NumParticles;
	GLuint m_glArray;
	GLuint m_glVerts;
	Vector3* m_Vertices;
	GLuint m_BackgroundSceneCubeMap;

	Matrix4 m_ModelMatrix;

	Shader* m_shdrParticles;
	Shader* m_shdrComputeNormals;
	Shader* m_shdrBilateralBlur;
	Shader* m_shdrFinalDisplay;

	GLuint m_fboParticle;
	GLuint m_texParticleNormals;
	GLuint m_texParticleAbsorbtion;
	GLuint m_texParticleDepth;

	GLuint m_texParticleAbsorbtion2;//For Blur Ping-Pong Blur
	GLuint m_texParticleDepth2; 
};