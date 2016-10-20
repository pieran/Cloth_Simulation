#pragma once

#include "ClothScenario.h"
#include <nclgl\Shader.h>
#include <nclgl\Mesh.h>
#include <ncltech\GameObject.h>

enum CutState { CUTSTATE_IDLE, CUTSTATE_INPROGRESS, CUTSTATE_FINISHED };


class ClothScenarioViewer : public GameObject
{
public:
	ClothScenarioViewer(int numSimultaneousViews);
	virtual ~ClothScenarioViewer();


	void AddScenario(const std::string& name, ClothScenario* scenario);

	void GotoNextScenario();
	void GotoPreviousScenario();

	void GotoScenario(const std::string& name);
	void GotoScenario(int idx);

	const std::string GetScenarioName() const;
	const std::string& GetScenarioName(int idx) const;

	int GetScenarioIndex() const;
	int GetScenarioIndex(const std::string& name) const;

	ClothScenario* GetScenario();
	ClothScenario* GetScenario(const std::string& name);
	ClothScenario* GetScenario(int idx);

	const ClothScenario* GetScenario() const;
	const ClothScenario* GetScenario(const std::string& name) const;
	const ClothScenario* GetScenario(int idx) const;

	bool DeleteScenario(const std::string& name);
	bool DeleteScenario(int idx);

	inline int GetNumScenarios() const { return (int)m_Scenarios.size(); }
	inline int GetNumViewedScenarios() const { return (int)m_NumScenariosInView; }

	inline void SetPaused(bool paused) { m_Paused = paused; }
	inline void TogglePaused() { m_Paused = !m_Paused; }
	inline bool IsPaused() const { return m_Paused; }


	void SingleTimeStep()
	{
		m_Paused = true;
		for (int i = 0; i < m_NumScenariosInView; ++i)
		{
			int scenarioIdx = (m_CurrentScenario + i) % m_Scenarios.size();

			ClothScenario* scenario = GetScenario(scenarioIdx);
			if (scenario != NULL)
			{
				scenario->OnUpdateScene(scenario->GetCloth()->Timestep());
			}
		}
	}
protected:
	void OnRenderObject() override;				//Handles OpenGL calls to Render the object
	void OnUpdateObject(float dt) override;		//Override to handle things like AI etc on update loop

	void LoadShaderConfiguration();

	void BufferGLArrays();
	void BuildGradientMap();

	void GetMouseRayCasts(Ray& start, Ray& end);
protected:
	bool m_Paused;
	bool m_DrawElements;
	bool m_DrawStrain;

	int m_NumScenariosInView;

	int m_CurrentScenario;
	std::vector<std::string> m_ScenarioLookups;
	std::vector<ClothScenario*> m_Scenarios;

	Shader* m_RenderShader;
	GLuint  m_glClothOutGradientMap;

	GLuint m_glClothTexture;

	bool m_ArrayObjInvalidated;
	PArray<int> m_IndicesArray;

	uint   m_NumVertices;
	uint   m_NumIndicies;
	GLuint m_glArrayObj;
	GLuint m_glBufferVerts;
	GLuint m_glBufferColours;
	GLuint m_glBufferTexCoords;
	GLuint m_glBufferNormals;
	GLuint m_glBufferIndices;

	CutState m_CutState;
	Vector2  m_CutStart;
	Vector2  m_CutEnd;
};

