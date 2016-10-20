#include "FluidScreenSpaceRenderer.h"
#include <ncltech\Scene.h>
#include <assert.h>

FluidScreenSpaceRenderer::FluidScreenSpaceRenderer(const std::string& name)
	: GameObject(name)
	, m_glArray(NULL)
	, m_glVerts(NULL)
	, m_shdrParticles(NULL)
	, m_shdrFinalDisplay(NULL)
	, m_shdrComputeNormals(NULL)
	, m_shdrBilateralBlur(NULL)
	, m_MaxParticles(0)
	, m_NumParticles(0)
	, m_Vertices(NULL)
	, m_fboParticle(0)
	, m_texParticleNormals(0)
	, m_texParticleAbsorbtion(0)
	, m_texParticleAbsorbtion2(0)
	, m_texParticleDepth(0)
	, m_texParticleDepth2(0)
	, m_BackgroundSceneCubeMap(0)
{
	this->SetBoundingRadius(FLT_MAX);
}

FluidScreenSpaceRenderer::~FluidScreenSpaceRenderer()
{
	m_MaxParticles = 0;
	m_NumParticles = 0;

	m_BackgroundSceneCubeMap = 0;

	if (m_Vertices)
	{
		delete[] m_Vertices;
		m_Vertices = NULL;
	}

	if (m_glArray)
	{
		glDeleteBuffers(1, &m_glVerts);
		glDeleteVertexArrays(1, &m_glArray);
		m_glArray = NULL;
	}

	if (m_shdrParticles)
	{
		delete m_shdrParticles;
		m_shdrParticles = NULL;

		delete m_shdrFinalDisplay;
		m_shdrFinalDisplay = NULL;

		delete m_shdrComputeNormals;
		m_shdrComputeNormals = NULL;

		delete m_shdrBilateralBlur;
		m_shdrBilateralBlur = NULL;
	}

	if (m_fboParticle)
	{
		glDeleteFramebuffers(1, &m_fboParticle);
		glDeleteTextures(1, &m_texParticleNormals);
		glDeleteTextures(1, &m_texParticleDepth);
		glDeleteTextures(1, &m_texParticleDepth2);
		glDeleteTextures(1, &m_texParticleAbsorbtion);
		glDeleteTextures(1, &m_texParticleAbsorbtion2);
	}
}

void FluidScreenSpaceRenderer::OnAttachedToScene()
{
	ReloadShaders();
	BuildFBO();
}

void FluidScreenSpaceRenderer::SetWorldSpace(const Vector3& world_origin, const Vector3& sim_ratio)
{
	m_ModelMatrix = Matrix4::Scale(sim_ratio);
	m_ModelMatrix.SetPositionVector(world_origin);
}

void FluidScreenSpaceRenderer::BuildBuffers(int max_particles)
{
	if (max_particles == m_MaxParticles)
		return;

	m_MaxParticles = max_particles;

	m_Vertices = new Vector3[m_MaxParticles];

	if (!m_glArray) glGenVertexArrays(1, &m_glArray);
	glBindVertexArray(m_glArray);

	//Buffer vertex data	
	if (!m_glVerts) glGenBuffers(1, &m_glVerts);
	glBindBuffer(GL_ARRAY_BUFFER, m_glVerts);
	glBufferData(GL_ARRAY_BUFFER, m_MaxParticles*sizeof(Vector3), NULL, GL_STREAM_DRAW);
	glVertexAttribPointer(VERTEX_BUFFER, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(VERTEX_BUFFER);

	glBindVertexArray(0);
}

void FluidScreenSpaceRenderer::UpdateBuffers(int num_particles, Particle* particles)
{
	m_NumParticles = num_particles;

#pragma omp parallel for
	for (int i = 0; i < m_NumParticles; ++i)
	{
		m_Vertices[i] = particles[i].pos;
	}

	glBindBuffer(GL_ARRAY_BUFFER, m_glVerts);
	glBufferSubData(GL_ARRAY_BUFFER, 0, m_NumParticles*sizeof(Vector3), m_Vertices);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void FluidScreenSpaceRenderer::UpdateBuffers(int num_particles, float4* particles)
{
	m_NumParticles = num_particles;

#pragma omp parallel for
	for (int i = 0; i < m_NumParticles; ++i)
	{
		m_Vertices[i].x = particles[i].x;
		m_Vertices[i].y = particles[i].y;
		m_Vertices[i].z = particles[i].z;
	}

	glBindBuffer(GL_ARRAY_BUFFER, m_glVerts);
	glBufferSubData(GL_ARRAY_BUFFER, 0, m_NumParticles*sizeof(Vector3), m_Vertices);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void FluidScreenSpaceRenderer::OnRenderObject()
{
	if (m_NumParticles == 0)
		return;
#define BLUR_ONE_PASS TRUE

	float particle_radius = 0.235f;// sqrt(0.4f * 0.4f + 0.4f * 0.4f + 0.4f * 0.4f);
	Matrix4 projMatrix = m_Scene->GetProjMatrix();
	Matrix4 viewMatrix = m_Scene->GetViewMatrix();
	Matrix4 invProjMatrix = Matrix4::Inverse(projMatrix);
	Matrix4 invProjViewMatrix = Matrix4::Inverse(projMatrix * viewMatrix);

	glBindVertexArray(m_glArray);

	glBindFramebuffer(GL_FRAMEBUFFER, m_fboParticle);
	glViewport(0, 0, m_Scene->GetWidth(), m_Scene->GetHeight());
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	
	//First draw the particles again, building up an eye space fluid depth map
#if BLUR_ONE_PASS == TRUE
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texParticleAbsorbtion2, 0);
#else
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texParticleAbsorbtion, 0);
#endif
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_Scene->m_ScreenDTex, 0);
	glClear(GL_COLOR_BUFFER_BIT);

	glUseProgram(m_shdrParticles->GetProgram());
	glUniformMatrix4fv(glGetUniformLocation(m_shdrParticles->GetProgram(), "mdlMatrix"), 1, false, (float*)&m_ModelMatrix);
	glUniformMatrix4fv(glGetUniformLocation(m_shdrParticles->GetProgram(), "viewMatrix"), 1, false, (float*)&viewMatrix);
	glUniformMatrix4fv(glGetUniformLocation(m_shdrParticles->GetProgram(), "projMatrix"), 1, false, (float*)&projMatrix);
	glUniform1f(glGetUniformLocation(m_shdrParticles->GetProgram(), "sphereRadius"), particle_radius);
	glUniform1f(glGetUniformLocation(m_shdrParticles->GetProgram(), "sphereAlpha"), 0.01f);

	glDepthMask(GL_FALSE);
	glBlendFunc(GL_ONE, GL_ONE);
	glDrawArrays(GL_POINTS, 0, m_NumParticles);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthMask(GL_TRUE);
	
	//Then Draw The Particles to the depth buffer
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texParticleNormals, 0);
#if BLUR_ONE_PASS == TRUE
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_texParticleDepth2, 0);
#else
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_texParticleDepth, 0);
#endif

	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	//glDepthFunc(GL_ALWAYS);
	glDrawBuffers(0, GL_NONE);
	glDrawArrays(GL_POINTS, 0, m_NumParticles);
	//glDepthFunc(GL_LEQUAL);
	


	//Blur the depth buffer
	//glDepthMask(GL_FALSE);
	GLenum buffers[] = { GL_COLOR_ATTACHMENT0 }; glDrawBuffers(1, buffers);

	Vector2 imgScale = Vector2(1.0f / (float)m_Scene->GetWidth(), 1.0f / (float)m_Scene->GetHeight());
	
	

	glUseProgram(m_shdrBilateralBlur->GetProgram());
	glUniform1i(glGetUniformLocation(m_shdrBilateralBlur->GetProgram(), "texDepth"), 0);
	glUniform1i(glGetUniformLocation(m_shdrBilateralBlur->GetProgram(), "texAbsorb"), 1);
	glUniformMatrix4fv(glGetUniformLocation(m_shdrBilateralBlur->GetProgram(), "projMatrix"), 1, false, (float*)&projMatrix);
	glUniformMatrix4fv(glGetUniformLocation(m_shdrBilateralBlur->GetProgram(), "invProjMatrix"), 1, false, (float*)&invProjMatrix);

#if BLUR_ONE_PASS == FALSE
	//Bind Depth Texture #2
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_texParticleDepth2, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texParticleAbsorbtion2, 0);
	//glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	glDepthFunc(GL_ALWAYS);

	glUniform2f(glGetUniformLocation(m_shdrBilateralBlur->GetProgram(), "tsBlurDir"), imgScale.x, 0.0f);
	glUniform2f(glGetUniformLocation(m_shdrBilateralBlur->GetProgram(), "esBlurDir"), 1.0f, 0.0f);
	glActiveTexture(GL_TEXTURE0);  glBindTexture(GL_TEXTURE_2D, m_texParticleDepth);
	glActiveTexture(GL_TEXTURE1);  glBindTexture(GL_TEXTURE_2D, m_texParticleAbsorbtion);
	glDrawArrays(GL_POINTS, 0, 1);
	glDepthFunc(GL_LEQUAL);
#endif


	//Bind Depth Texture #1 (Original)
	//GLenum buffers[] = { GL_COLOR_ATTACHMENT0 };
	//glDrawBuffers(1, buffers);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_texParticleDepth, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texParticleAbsorbtion, 0);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	glUniform2f(glGetUniformLocation(m_shdrBilateralBlur->GetProgram(), "tsBlurDir"), 0.0f, imgScale.y);
	glUniform2f(glGetUniformLocation(m_shdrBilateralBlur->GetProgram(), "esBlurDir"), 0.0f, 1.0f);
	glActiveTexture(GL_TEXTURE0);  glBindTexture(GL_TEXTURE_2D, m_texParticleDepth2);
	glActiveTexture(GL_TEXTURE1);  glBindTexture(GL_TEXTURE_2D, m_texParticleAbsorbtion2);
	//glDepthFunc(GL_ALWAYS);
	glDrawArrays(GL_POINTS, 0, 1);
	//glDepthFunc(GL_LEQUAL);

	//Compute surface normals
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texParticleNormals, 0);

	
	glUseProgram(m_shdrComputeNormals->GetProgram());
	glUniform2fv(glGetUniformLocation(m_shdrComputeNormals->GetProgram(), "texelSize"), 1, &imgScale.x);
	glUniform1i(glGetUniformLocation(m_shdrComputeNormals->GetProgram(), "texDepth"), 1);
	glUniformMatrix4fv(glGetUniformLocation(m_shdrComputeNormals->GetProgram(), "invProjMatrix"), 1, false, (float*)&invProjViewMatrix);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, m_texParticleDepth);

	glDepthMask(GL_FALSE);
	glDrawArrays(GL_POINTS, 0, 1);
	glDepthMask(GL_TRUE);


	//Render Fullscreen Quad
	glBindFramebuffer(GL_FRAMEBUFFER, m_Scene->m_ScreenFBO);
	//glDisable(GL_CULL_FACE);
	Vector3 camPos = m_Scene->m_Camera->GetPosition();

	glUseProgram(m_shdrFinalDisplay->GetProgram());
	glUniformMatrix4fv(glGetUniformLocation(m_shdrFinalDisplay->GetProgram(), "invProjViewMatrix"), 1, false, (float*)&invProjViewMatrix);
	glUniform3fv(glGetUniformLocation(m_shdrFinalDisplay->GetProgram(), "ambientColour"), 1, &m_Scene->m_AmbientColour.x);
	glUniform3fv(glGetUniformLocation(m_shdrFinalDisplay->GetProgram(), "invLightDir"), 1, &m_Scene->m_InvLightDirection.x);
	glUniform3fv(glGetUniformLocation(m_shdrFinalDisplay->GetProgram(), "cameraPos"), 1, &camPos.x);
	glUniform1i(glGetUniformLocation(m_shdrFinalDisplay->GetProgram(), "texNormals"), 0);
	glUniform1i(glGetUniformLocation(m_shdrFinalDisplay->GetProgram(), "texAbsorbtion"), 1);
	glUniform1i(glGetUniformLocation(m_shdrFinalDisplay->GetProgram(), "texDepth"), 2);
	glUniform1i(glGetUniformLocation(m_shdrFinalDisplay->GetProgram(), "texReflect"), 3);
	
	glActiveTexture(GL_TEXTURE3); glBindTexture(GL_TEXTURE_CUBE_MAP, m_BackgroundSceneCubeMap);
	glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, m_texParticleDepth);
	glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, m_texParticleAbsorbtion);	
	glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D,  m_texParticleNormals);



	glDrawArrays(GL_POINTS, 0, 1);
	glBindVertexArray(0);

	




	glClearColor(0.6f, 0.6f, 0.6f, 1.0f);

	if (m_Scene->GetCurrentShader() != NULL)
		glUseProgram(m_Scene->GetCurrentShader()->GetProgram());
	else
		glUseProgram(0);
}

void FluidScreenSpaceRenderer::OnUpdateObject(float dt)
{
	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_X))
	{
		ReloadShaders();
		BuildFBO();
	}
}

void BuildTextureComponent(GLuint* texture, GLint fbo_component, GLint internalFormat, bool linear, bool repeat, int width, int height) {
	if (*texture == NULL)
		glGenTextures(1, texture);

	glBindTexture(GL_TEXTURE_2D, *texture);
	glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, width, height, 0,
		(fbo_component == GL_DEPTH_ATTACHMENT) ? GL_DEPTH_COMPONENT : GL_RGB,
		(fbo_component == GL_DEPTH_ATTACHMENT) ? GL_FLOAT : GL_UNSIGNED_BYTE, NULL);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, linear ? GL_LINEAR : GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, linear ? GL_LINEAR : GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, repeat ? GL_REPEAT : GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, repeat ? GL_REPEAT : GL_CLAMP_TO_EDGE);

	if (fbo_component == GL_DEPTH_ATTACHMENT) {

		//Has to be GL_NONE in order to correctly sample from glsl shader, 
		//	enable again to allow for hardware optimisation (but without the chance to correctly sample manually)
		//glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);


		//Handle depth as a normal texture/sampler2D
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	}

	glBindTexture(GL_TEXTURE_2D, 0);
}

void FluidScreenSpaceRenderer::BuildFBO()
{
	int width = m_Scene->GetWidth();
	int height = m_Scene->GetHeight();

	BuildTextureComponent(&m_texParticleNormals, GL_COLOR_ATTACHMENT0, GL_RGBA8, false, false, width, height);
	BuildTextureComponent(&m_texParticleAbsorbtion, GL_COLOR_ATTACHMENT0, GL_R32F, false, false, width, height);
	BuildTextureComponent(&m_texParticleAbsorbtion2, GL_COLOR_ATTACHMENT0, GL_R32F, false, false, width, height);
	BuildTextureComponent(&m_texParticleDepth, GL_DEPTH_ATTACHMENT, GL_DEPTH_COMPONENT32, false, false, width, height);
	BuildTextureComponent(&m_texParticleDepth2, GL_DEPTH_ATTACHMENT, GL_DEPTH_COMPONENT32, false, false, width, height);


	if (m_fboParticle == NULL) {
		glGenFramebuffers(1, &m_fboParticle);
	}

	glBindFramebuffer(GL_FRAMEBUFFER, m_fboParticle);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_texParticleDepth, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texParticleNormals, 0);

	GLenum buffers[] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, buffers);


	GLenum status;
	status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status != GL_FRAMEBUFFER_COMPLETE) {
		printf("FRAMEBUFFER ERROR: %d\n", status);
		assert(false);
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void FluidScreenSpaceRenderer::ReloadShaders()
{
	if (m_shdrParticles) delete m_shdrParticles;
	if (m_shdrComputeNormals) delete m_shdrComputeNormals;
	if (m_shdrBilateralBlur) delete m_shdrBilateralBlur;
	if (m_shdrFinalDisplay) delete m_shdrFinalDisplay;


	m_shdrParticles = new Shader(SHADERDIR"fluid/Particle.vert", SHADERDIR"fluid/Particle.frag", SHADERDIR"fluid/Particle.geom");
	m_shdrParticles->LinkProgram();

	m_shdrComputeNormals = new Shader(SHADERDIR"fluid/FullScreen.vert", SHADERDIR"fluid/ComputeNormals.frag", SHADERDIR"fluid/FullScreen.geom");
	m_shdrComputeNormals->LinkProgram();

	m_shdrBilateralBlur = new Shader(SHADERDIR"fluid/FullScreen.vert", SHADERDIR"fluid/BilateralFilter.frag", SHADERDIR"fluid/FullScreen.geom");
	m_shdrBilateralBlur->LinkProgram();

	m_shdrFinalDisplay = new Shader(SHADERDIR"fluid/FullScreen.vert", SHADERDIR"fluid/FullScreen.frag", SHADERDIR"fluid/FullScreen.geom");
	m_shdrFinalDisplay->LinkProgram();
}