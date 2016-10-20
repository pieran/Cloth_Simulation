#include "ClothScenarioViewer.h"
#include "FEBase.h"
#include <ncltech\Scene.h>
#include <ncltech\NCLDebug.h>
#include <assert.h>
#include <sstream>

ClothScenarioViewer::ClothScenarioViewer(int numSimultaneousViews)
{
	m_Paused = true;
	
	m_CutState = CUTSTATE_IDLE;

	m_CurrentScenario = -1;
	m_ScenarioLookups.reserve(10);
	m_Scenarios.reserve(10);
	m_NumScenariosInView = numSimultaneousViews;
	m_ArrayObjInvalidated = false;

	m_BoundingRadius = FLT_MAX;
	m_DrawElements = false;
	m_DrawStrain = false;

	glGenVertexArrays(1, &m_glArrayObj);
	glBindVertexArray(m_glArrayObj);
	glGenBuffers(1, &m_glBufferVerts);
	glGenBuffers(1, &m_glBufferColours);
	glGenBuffers(1, &m_glBufferNormals);
	glGenBuffers(1, &m_glBufferIndices);
	glGenBuffers(1, &m_glBufferTexCoords);
	glBindVertexArray(0);

	//m_glClothTexture = SOIL_load_OGL_texture(TEXTUREDIR"checkerboard.tga", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);
	m_glClothTexture = SOIL_load_OGL_texture(TEXTUREDIR"cloth_seamless_yellow.tga", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);

	glBindTexture(GL_TEXTURE_2D, m_glClothTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST); //No linear interpolation to get crisp checkerboard no matter the scalling
	glGenerateMipmap(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);

	BuildGradientMap();

	m_RenderShader = new Shader(SHADERDIR"cloth/cloth_vertex.glsl", SHADERDIR"cloth/cloth_fragment_spec.glsl");
	if (!m_RenderShader->LinkProgram())
	{
		assert(false);
	}

	glUseProgram(m_RenderShader->GetProgram());
	glUniform1i(glGetUniformLocation(m_RenderShader->GetProgram(), "colorGradient"), 4);
}

ClothScenarioViewer::~ClothScenarioViewer()
{
	if (m_glArrayObj)
	{
		glDeleteVertexArrays(1, &m_glArrayObj);
		glDeleteBuffers(1, &m_glBufferVerts);
		glDeleteBuffers(1, &m_glBufferColours);
		glDeleteBuffers(1, &m_glBufferNormals);
		glDeleteBuffers(1, &m_glBufferIndices);
		glDeleteBuffers(1, &m_glBufferTexCoords);

		glDeleteTextures(1, &m_glClothTexture);
		glDeleteTextures(1, &m_glClothOutGradientMap);

		delete m_RenderShader;
	}

	for (uint i = 0; i < m_Scenarios.size(); ++i)
	{
		m_Scenarios[i]->OnDeleteScene();
		//delete m_Scenarios[i];
	}
	m_Scenarios.clear();
	m_ScenarioLookups.clear();
}


void ClothScenarioViewer::AddScenario(const std::string& name, ClothScenario* scenario)
{
	assert(GetScenarioIndex(name) == -1);

	m_ScenarioLookups.push_back(name);
	m_Scenarios.push_back(scenario);

	scenario->OnInitializeScene();

	this->AddChildObject(scenario);

	if (m_CurrentScenario == -1)
		GotoScenario(0);
}

void ClothScenarioViewer::GotoNextScenario()
{
	GotoScenario((m_CurrentScenario + 1) % m_Scenarios.size());
}

void ClothScenarioViewer::GotoPreviousScenario()
{
	GotoScenario((m_CurrentScenario -1) < 0 ? m_Scenarios.size() - 1 : m_CurrentScenario - 1);
}

void ClothScenarioViewer::GotoScenario(const std::string& name)
{
	int idx = GetScenarioIndex(name);
	if (idx > -1)
	{
		GotoScenario(idx);
	}
}

void ClothScenarioViewer::GotoScenario(int idx)
{
	assert(idx >= 0 && idx < m_Scenarios.size());
	//Cleanup old Scenario
	/*if (m_CurrentScenario > -1)
	{
		for (int i = 0; i < m_NumScenariosInView; ++i)
		{
			int idx = (m_CurrentScenario + i) % m_Scenarios.size();

			if (m_Scenarios[idx] != NULL)
				m_Scenarios[idx]->OnDeleteScene();
		}
	}*/
	
	//Initialize New Scenario
	m_CurrentScenario = idx;
	//for (int i = 0; i < m_NumScenariosInView; ++i)
	//{
	//	int idx = (m_CurrentScenario + i) % m_Scenarios.size();

	//	if (m_Scenarios[idx] != NULL)
	//		m_Scenarios[idx]->OnInitializeScene();
	//}
}

const std::string ClothScenarioViewer::GetScenarioName() const
{
	if (m_CurrentScenario == -1)
		return undefined_string;
	else
	{
		std::stringstream name;
		for (int i = 0; i < m_NumScenariosInView; ++i)
		{
			if (i > 0)
				name << ", ";

			int idx = (m_CurrentScenario + i) % m_Scenarios.size();

			if (m_Scenarios[idx] != NULL)
				name << "(" << m_ScenarioLookups[idx] << ")";
			else
				name << "(" << undefined_string << ")";
		}
		return name.str();
	}
}

const std::string& ClothScenarioViewer::GetScenarioName(int idx) const
{
	return (m_CurrentScenario > -1) ? m_ScenarioLookups[idx] : undefined_string;
}

int ClothScenarioViewer::GetScenarioIndex() const
{
	return m_CurrentScenario;
}

int ClothScenarioViewer::GetScenarioIndex(const std::string& name) const
{
	for (uint i = 0; i < m_ScenarioLookups.size(); ++i)
	{
		if (m_ScenarioLookups[i].compare(name) == 0)
			return i;
	}

	return -1;
}

ClothScenario* ClothScenarioViewer::GetScenario()
{
	return (m_CurrentScenario == -1) ? NULL : m_Scenarios[m_CurrentScenario];
}

ClothScenario* ClothScenarioViewer::GetScenario(const std::string& name)
{
	int idx = GetScenarioIndex(name);
	
	return (idx == -1)
		? NULL
		: m_Scenarios[idx];
}

ClothScenario* ClothScenarioViewer::GetScenario(int idx)
{
	assert(idx >= 0 && idx < m_Scenarios.size());
	return m_Scenarios[idx];
}

const ClothScenario* ClothScenarioViewer::GetScenario() const
{
	return (m_CurrentScenario == -1) ? NULL : m_Scenarios[m_CurrentScenario];
}

const ClothScenario* ClothScenarioViewer::GetScenario(const std::string& name) const
{
	int idx = GetScenarioIndex(name);

	return (idx == -1)
		? NULL
		: m_Scenarios[idx];
}

const ClothScenario* ClothScenarioViewer::GetScenario(int idx) const
{
	assert(idx >= 0 && idx < m_Scenarios.size());
	return m_Scenarios[idx];
}

bool ClothScenarioViewer::DeleteScenario(const std::string& name)
{
	int idx = GetScenarioIndex(name);
	if (idx != -1)
	{
		return DeleteScenario(idx);
	}
	return false;
}

bool ClothScenarioViewer::DeleteScenario(int idx)
{
	if (idx < 0 || idx >= m_Scenarios.size())
		return false; //INVALID ARRAY INDEX

	delete m_Scenarios[idx];

	m_ScenarioLookups.erase(m_ScenarioLookups.begin() + idx);
	m_Scenarios.erase(m_Scenarios.begin() + idx);

	if (m_CurrentScenario > idx || (m_CurrentScenario == idx == m_Scenarios.size() == 0))
	{
		m_CurrentScenario--;
	}
}


void ClothScenarioViewer::LoadShaderConfiguration()
{
	glDisable(GL_CULL_FACE);
	glUseProgram(m_RenderShader->GetProgram());

	Vector3 camPos = m_Scene->m_Camera->GetPosition();
	glUniformMatrix4fv(glGetUniformLocation(m_RenderShader->GetProgram(), "modelMatrix"), 1, false, (float*)&m_Scene->modelMatrix);
	glUniformMatrix4fv(glGetUniformLocation(m_RenderShader->GetProgram(), "viewMatrix"), 1, false, (float*)&m_Scene->viewMatrix);
	glUniformMatrix4fv(glGetUniformLocation(m_RenderShader->GetProgram(), "projMatrix"), 1, false, (float*)&m_Scene->projMatrix);
	glUniformMatrix4fv(glGetUniformLocation(m_RenderShader->GetProgram(), "textureMatrix"), 1, false, (float*)&m_Scene->textureMatrix);
	glUniform1i(glGetUniformLocation(m_RenderShader->GetProgram(), "diffuseTex"), 0);
	glUniform3fv(glGetUniformLocation(m_RenderShader->GetProgram(), "ambientColour"), 1, &m_Scene->m_AmbientColour.x);
	glUniform3fv(glGetUniformLocation(m_RenderShader->GetProgram(), "invLightDir"), 1, &m_Scene->m_InvLightDirection.x);
	glUniform3fv(glGetUniformLocation(m_RenderShader->GetProgram(), "cameraPos"), 1, &camPos.x);
	glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "specularIntensity"), m_Scene->m_SpecularIntensity);

	glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "isTextured"), m_DrawStrain ? 0.0f : 1.0f);
	glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "minVal"), 0.0f);// m_Cloth->m_RenderValueMinMax.x);
	glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "maxVal"), 5.0f);// m_Cloth->m_RenderValueMinMax.y);


																					//NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Min/Max: %5.2f / %5.2f", m_Cloth->m_RenderValueMinMax.x, m_Cloth->m_RenderValueMinMax.y);
	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_1D, m_glClothOutGradientMap);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, m_glClothTexture);
}


void ClothScenarioViewer::OnRenderObject()
{
	GameObject::OnRenderObject();

	//Cutting
	if (m_CutState == CUTSTATE_INPROGRESS)
	{
		
		Ray sRay, eRay;
		GetMouseRayCasts(sRay, eRay);

		NCLDebug::DrawHairLine(sRay.wsPos, eRay.wsPos, Vector4(1.0f, 0.8f, 1.2f, 1.0f));
	}
	
	ClothScenario* scenario = GetScenario();
	if (scenario != NULL)
	{

		GLboolean cullface;
		glGetBooleanv(GL_CULL_FACE, &cullface);

		GLint pid;
		glGetIntegerv(GL_CURRENT_PROGRAM, &pid);

		bool useCustomShader = m_Scene->GetCurrentRenderStage() != RENDER_TYPE_SHADOW_EXTRUDE;
		if (useCustomShader) LoadShaderConfiguration();

		glBindVertexArray(m_glArrayObj);
		if (useCustomShader) glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "isTextured"), 1.0f);
		glDrawElements(GL_TRIANGLES, m_NumIndicies, GL_UNSIGNED_INT, 0);

		if (m_DrawElements)
		{
			LoadShaderConfiguration();
			if (Window::GetKeyboard()->KeyDown(KEYBOARD_T))
			{
				glDisable(GL_DEPTH_TEST);
				glPointSize(11.0f);
				glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "isTextured"), 0.0f);
				glDrawArrays(GL_POINTS, 0, scenario->m_ClothSim->m_NumPhyxels);
				glEnable(GL_DEPTH_TEST);
			}
			else
			{
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				glUniform1f(glGetUniformLocation(m_RenderShader->GetProgram(), "isTextured"), 0.0f);
				glDrawElements(GL_TRIANGLES, m_NumIndicies, GL_UNSIGNED_INT, 0);
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			}
		}
		glBindVertexArray(0);

		if (cullface)
			glEnable(GL_CULL_FACE);

		glUseProgram(pid);
	}
}

void ClothScenarioViewer::OnUpdateObject(float dt)
{
	GameObject::OnUpdateObject(dt);




	//CUTTING
	if (Window::GetKeyboard()->KeyDown(KEYBOARD_CONTROL))
	{
		bool mouseDown = Window::GetMouse()->ButtonDown(MOUSE_LEFT);
		bool mouseHeld = Window::GetMouse()->ButtonHeld(MOUSE_LEFT);
		
		m_CutEnd = Window::GetWindow().GetMousePosition();
		if (mouseDown && !mouseHeld)
		{
			m_CutState = CUTSTATE_INPROGRESS;
			m_CutStart = m_CutEnd;
		}
		else if (!mouseDown && m_CutState == CUTSTATE_INPROGRESS)
		{
			m_CutState = CUTSTATE_FINISHED;

			Ray start, end;
			int idx_start, idx_end;
			GetMouseRayCasts(start, end);
			
			Vector3 crossDir = (end.wsPos - start.wsPos); // crossDir.Normalise();
			Vector3 normal = Vector3::Cross(start.wsDir, crossDir); normal.Normalise();

			Plane plane;
			plane.SetNormal(normal);
			plane.SetDistance(-Vector3::Dot(start.wsPos, normal));

			//Trigger On Cut Method!!!
			bool has_cut = false;
			for (int i = 0; i < m_NumScenariosInView; ++i)
			{
				int scenarioIdx = (m_CurrentScenario + i) % m_Scenarios.size();

				ClothScenario* scenario = GetScenario(scenarioIdx);
				if (scenario != NULL)
				{
					if (scenario->GetCloth()->GetClosestPointsForCutting(start, end, &idx_start, &idx_end))
					{
						scenario->OnClothCut(idx_start, idx_end, start, end, plane);
						has_cut = true;
					}
				}
			}

			if (has_cut) m_ArrayObjInvalidated = true;


		}
		else if (!mouseDown)
		{
			m_CutState = CUTSTATE_IDLE;
		}
	}
	else
	{
		m_CutState = CUTSTATE_IDLE;
	}

	//SCENARIO UPDATES
	for (int i = 0; i < m_NumScenariosInView; ++i)
	{
		int scenarioIdx = (m_CurrentScenario + i) % m_Scenarios.size();

		ClothScenario* scenario = GetScenario(scenarioIdx);
		if (scenario != NULL)
		{
			if (!m_Paused)
				scenario->OnUpdateScene(1.0f / 60.0f);
		}
	}

	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_X))
		m_DrawStrain = !m_DrawStrain;

	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_F))
		m_DrawElements = !m_DrawElements;

	BufferGLArrays();
}



void ClothScenarioViewer::BufferGLArrays()
{
	int totalVerts		= 0;
	int totalIndices	= 0;

	for (int i = 0; i < m_NumScenariosInView; ++i)
	{
		int scenarioIdx = (m_CurrentScenario + i) % m_Scenarios.size();
		ClothScenario* scenario = GetScenario(scenarioIdx);
		if (scenario != NULL)
		{
			ClothBase* base = scenario->m_ClothSim;
			totalVerts += base->m_NumPhyxels;
			totalIndices += base->m_NumTris * 3;
		}
	}

	if (totalVerts != 0)
	{
		glBindVertexArray(m_glArrayObj);
		if (m_ArrayObjInvalidated || totalVerts != m_NumVertices)
		{
			m_NumVertices = totalVerts;

			//Buffer vertex data
			glBindBuffer(GL_ARRAY_BUFFER, m_glBufferVerts);
			glBufferData(GL_ARRAY_BUFFER, m_NumVertices*sizeof(Vector3), NULL, GL_STATIC_DRAW);
			glVertexAttribPointer(VERTEX_BUFFER, 3, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(VERTEX_BUFFER);

			//Buffer vertex colour data	
			glBindBuffer(GL_ARRAY_BUFFER, m_glBufferColours);
			glBufferData(GL_ARRAY_BUFFER, m_NumVertices*sizeof(float), NULL, GL_STATIC_DRAW);
			glVertexAttribPointer(COLOUR_BUFFER, 1, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(COLOUR_BUFFER);

			//Buffer texture coord data
			glBindBuffer(GL_ARRAY_BUFFER, m_glBufferTexCoords);
			glBufferData(GL_ARRAY_BUFFER, m_NumVertices*sizeof(Vector2), NULL, GL_STATIC_DRAW);
			glVertexAttribPointer(TEXTURE_BUFFER, 2, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(TEXTURE_BUFFER);

			//Buffer normal data
			glBindBuffer(GL_ARRAY_BUFFER, m_glBufferNormals);
			glBufferData(GL_ARRAY_BUFFER, m_NumVertices*sizeof(Vector3), NULL, GL_STATIC_DRAW);
			glVertexAttribPointer(NORMAL_BUFFER, 3, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(NORMAL_BUFFER);
		}

		if (m_ArrayObjInvalidated || totalIndices != m_NumIndicies)
		{
			m_NumIndicies = totalIndices;

			m_IndicesArray.resize(m_NumIndicies);
		
			totalIndices = 0;
			totalVerts = 0;
			for (int i = 0; i < m_NumScenariosInView; ++i)
			{
				int scenarioIdx = (m_CurrentScenario + i) % m_Scenarios.size();
				ClothScenario* scenario = GetScenario(scenarioIdx);
				if (scenario != NULL)
				{
					ClothBase* base = scenario->m_ClothSim;
					uint numIndices = base->m_NumTris * 3;
					uint numVerts = base->m_NumPhyxels;

					for (uint i = 0; i < numIndices; ++i)
						m_IndicesArray[totalIndices+i] = base->m_TriIndicies[i] + totalVerts;

					totalIndices += numIndices;
					totalVerts += numVerts;
				}
			}

			//buffer index data
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_glBufferIndices);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_NumIndicies*sizeof(GLuint), &m_IndicesArray[0], GL_STATIC_DRAW);
		}


		totalVerts = 0;
		for (int i = 0; i < m_NumScenariosInView; ++i)
		{
			int scenarioIdx = (m_CurrentScenario + i) % m_Scenarios.size();
			ClothScenario* scenario = GetScenario(scenarioIdx);
			if (scenario != NULL)
			{
				ClothBase* base = scenario->m_ClothSim;
				int numVerts = base->m_NumPhyxels;

				//Buffer vertex data	
				glBindBuffer(GL_ARRAY_BUFFER, m_glBufferVerts);
				glBufferSubData(GL_ARRAY_BUFFER, totalVerts * sizeof(Vector3), numVerts * sizeof(Vector3), &base->m_ClothStatesPos[0]);

				//Buffer vertex colour data	
				glBindBuffer(GL_ARRAY_BUFFER, m_glBufferColours);
				glBufferSubData(GL_ARRAY_BUFFER, totalVerts * sizeof(float), numVerts * sizeof(float), &base->m_RenderValue[0]);

				//Buffer texture coord data
				glBindBuffer(GL_ARRAY_BUFFER, m_glBufferTexCoords);
				glBufferSubData(GL_ARRAY_BUFFER, totalVerts * sizeof(Vector2), numVerts * sizeof(Vector2), &base->m_PhyxelTexCoords[0]);

				//Buffer normal data
				glBindBuffer(GL_ARRAY_BUFFER, m_glBufferNormals);
				glBufferSubData(GL_ARRAY_BUFFER, totalVerts * sizeof(Vector3), numVerts * sizeof(Vector3), &base->m_PhyxelNormals[0]);

				totalVerts += numVerts;
			}
		}

		glBindVertexArray(0);
	}

	m_ArrayObjInvalidated = false;
}

void ClothScenarioViewer::BuildGradientMap()
{
	const uint num_vals = 360;

	glGenTextures(1, &m_glClothOutGradientMap);
	glBindTexture(GL_TEXTURE_1D, m_glClothOutGradientMap);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	Vector3* gradient = new Vector3[num_vals];
	float step = 1.0f / num_vals * 0.25f;
	for (uint i = 0; i < num_vals; ++i)
	{
		hsv2rgb(gradient[i], Vector3(i * 0.25f + 270.0f, 1.0f, 1.0f));
	}

	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB32F, num_vals, 0, GL_RGB, GL_FLOAT, gradient);

	delete[] gradient;

	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_1D, 0);
}

void ClothScenarioViewer::GetMouseRayCasts(Ray& start, Ray& end)
{
	Vector2 screen_size = Window::GetWindow().GetScreenSize();
	Matrix4 invProjView = Matrix4::Inverse(m_Scene->GetProjViewMatrix());

	auto buildRay = [&](Vector2& mPos, Ray* out) {
		Vector2 clip_pos = Vector2(mPos.x / screen_size.x, mPos.y / screen_size.y);
		clip_pos.y = 1.0f - clip_pos.y; //Opengl has flipped Y
		clip_pos = (clip_pos * 2.0f) - 1.0f;

		Vector3 ws_pos_far = invProjView * Vector3(clip_pos.x, clip_pos.y, -1.0f);
		Vector3 ws_pos_near = invProjView * Vector3(clip_pos.x, clip_pos.y, 0.0f);

		out->wsPos = ws_pos_near;
		out->wsDir = (ws_pos_near - ws_pos_far);
		out->wsDir.Normalise();
	};

	buildRay(m_CutStart, &start);
	buildRay(m_CutEnd, &end);
}