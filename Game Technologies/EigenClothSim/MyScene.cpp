
#include "MyScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\DistanceConstraint.h>

#include "PBD.h"
#include "XPBD.h"
#include "FEDefault.h"
#include "FEImprovedRotations.h"

#include "ScenarioStretchVertical.h"
#include "ScenarioParachute.h"
#include "ScenarioCutting.h"
#include "ScenarioCollisionYusef.h"

MyScene::MyScene(Window& window) : Scene(window), m_ClothRenderer(NULL)
{

	m_InvLightDirection = Vector3(0.5f, 1.0f, 1.0f);
	m_InvLightDirection.Normalise();

	if (init == true)
		init = InitialiseGL();

	UpdateWorldMatrices(m_RootGameObject, Matrix4());

	m_Encoder = new VideoEncoder();

	ifstream myfile;
	myfile.open("videos\\index.txt");
	if (!myfile.is_open())
	{
		m_EncoderVideoNum = 0;
	}
	else
	{
		myfile >> m_EncoderVideoNum;
		myfile.close();
	}

	m_HudVisibility = 0;
	m_ClothUpdateMs = 0.0f;

	//Profiling
	m_GraphObject = new GraphObject();
	this->AddGameObject(m_GraphObject);
	m_GraphObject->AddSeries("Rotations", Vector4(0.7f, 0.5f, 1.0f, 1.0f));
	m_GraphObject->AddSeries("K_Matrix", Vector4(1.0f, 0.5f, 0.8f, 1.0f));
	m_GraphObject->AddSeries("Solver", Vector4(0.5f, 0.8f, 1.0f, 1.0f));
	m_GraphObject->AddSeries("Collisions", Vector4(0.8f, 1.0f, 0.5f, 1.0f));
	m_GraphObject->Parameters().title_text = "Cloth Computation Time";
	m_GraphObject->Parameters().axis_y_span_label_interval = 100;
	m_GraphObject->Parameters().axis_y_span_sub_interval = 25;
	m_GraphObject->Parameters().axis_y_span_max = 200;
	m_GraphObject->Parameters().axis_y_span_max_isdynamic = true;
	m_GraphObject->Parameters().axis_y_label_subfix = "ms";
	

	m_GraphObjectSolver = new GraphObject();
	this->AddGameObject(m_GraphObjectSolver);
	m_GraphObjectSolver->AddSeries("Init", Vector4(0.7f, 0.5f, 1.0f, 1.0f));
	m_GraphObjectSolver->AddSeries("A Mtx Mult", Vector4(1.0f, 0.5f, 0.8f, 1.0f));
	m_GraphObjectSolver->AddSeries("Dot Products", Vector4(0.5f, 0.8f, 1.0f, 1.0f));
	m_GraphObjectSolver->Parameters().title_text = "Conjugate Gradient";
	
	
	OnResize();
}

MyScene::~MyScene()
{
	delete m_Encoder;
}

void MyScene::OnResize()
{
	//Profiling
	m_GraphObject->Parameters().graph_width = 300;
	m_GraphObject->Parameters().graph_height = 150;
	m_GraphObject->SetPosition(width - 10, 30);


	m_GraphObjectSolver->Parameters().graph_width = 200;
	m_GraphObjectSolver->Parameters().graph_height = 100;
	m_GraphObjectSolver->SetPosition(width - 10, 280);
}


void MyScene::InitialiseCloth()
{
	m_ClothViewer = new ClothScenarioViewer(1);
	//m_ClothViewer->AddScenario("Impr Rotations (smth) - Timestep: 1/120", new ScenarioParachute(new FEImprovedRotations(false), 11, 11, Vector2(2.0f, 2.0f), Vector3(-2.0f, 1.0f, 0.0f)));
	//m_ClothViewer->AddScenario("Impr Rotations        - Timestep: 1/120", new ScenarioParachute(new FEImprovedRotations(false), 11, 11, Vector2(2.0f, 2.0f), Vector3(0.0f, 1.0f, 0.0f)));
	//m_ClothViewer->AddScenario("Default               - Timestep: 1/120", new ScenarioParachute(new FEDefault(), 11, 11, Vector2(2.0f, 2.0f), Vector3(2.0f, 1.0f, 0.0f)));

	//m_ClothViewer->GetScenario(0)->GetCloth()->SetTimestep(1.0f / 120.0f);
	//m_ClothViewer->GetScenario(1)->GetCloth()->SetTimestep(1.0f / 120.0f);
	//m_ClothViewer->GetScenario(2)->GetCloth()->SetTimestep(1.0f / 120.0f);



	//m_ClothViewer->AddScenario("PBD 51x51 0.125dt", new ScenarioCollisionYusef(new PBD(), 51, 51, Vector2(1.0f, 1.0f), Vector3(-1.5f, 3.0f, -1.5f)));
	//m_ClothViewer->AddScenario("PBD 51x51 0.25dt", new ScenarioCollisionYusef(new PBD(), 101, 101, Vector2(1.0f, 1.0f), Vector3(0.0f, 3.0f, -1.5f)));
	//m_ClothViewer->AddScenario("PBD 51x51 0.55dt", new ScenarioCollisionYusef(new PBD(), 151, 151, Vector2(1.0f, 1.0f), Vector3(1.5f, 3.0f, -1.5f)));

	//m_ClothViewer->AddScenario("FEM 11x11 0.5dt", new ScenarioParachute(new FEDefault(), 11, 11, Vector2(2.0f, 2.0f), Vector3(-1.0f, 1.0f, -1.5f)));
	//m_ClothViewer->AddScenario("FEM 51x51 0.25dt", new ScenarioCollisionYusef(new FEDefault(), 51, 51, Vector2(1.0f, 1.0f), Vector3(1.0f, 1.0f, -1.5f)));
	m_ClothViewer->AddScenario("FEM 51x51 0.25dt", new ScenarioCollisionYusef(new XPBD(), 64, 64, Vector2(1.0f, 1.0f), Vector3(0.0f, 3.0f, -1.5f)));
	//m_ClothViewer->AddScenario("FEM 101x101 0.25dt", new ScenarioCollisionYusef(new FEImprovedRotations(true), 101, 101, Vector2(1.0f, 1.0f), Vector3(0.0f, 3.0f, -1.5f)));
	//m_ClothViewer->AddScenario("FEM 151x151 0.25dt", new ScenarioCollisionYusef(new FEImprovedRotations(true), 201, 201, Vector2(1.0f, 1.0f), Vector3(1.5f, 3.0f, -1.5f)));

	//m_ClothViewer->AddScenario("FEM 51x51 0.5dt", new ScenarioCutting(new FEDefault(), 141, 141, Vector2(2.0f, 2.0f), Vector3(-1.0f, 1.0f, -1.5f)));
	//m_ClothViewer->AddScenario("FEM 51x51 0.125dt", new ScenarioParachute(new FEDefault(), 101, 101, Vector2(1.0f, 1.0f), Vector3(1.0f, 1.0f, -1.5f)));
	
	//m_ClothViewer->AddScenario("FEM 51x51 0.25dt", new ScenarioParachute(new FEDefault(), 101, 101, Vector2(1.0f, 1.0f), Vector3(1.0f, 1.0f, -1.5f)));



	m_ClothViewer->GetScenario(0)->GetCloth()->Timestep() = (1.0f / 60.0f);
	//m_ClothViewer->GetScenario(1)->GetCloth()->Timestep() = (1.0f / 1800.0f);
	//m_ClothViewer->GetScenario(2)->GetCloth()->Timestep() = (1.0f / 1800.0f);
	//m_ClothViewer->GetScenario(1)->GetCloth()->Timestep() = (1.0f / 900.0f);
	//m_ClothViewer->GetScenario(2)->GetCloth()->Timestep() = (1.0f / 120.0f);
	//m_ClothViewer->GetScenario(3)->GetCloth()->Timestep() = (1.0f / 120.0f);

	//m_ClothViewer->GetScenario(0)->GetCloth()->SetIntegrator(CLOTH_INTEGRATION_RK4);
	//m_ClothViewer->GetScenario(1)->GetCloth()->AdvanceVelocityHalfTimestep();
	//m_ClothViewer->GetScenario(2)->GetCloth()->SetIntegrator(CLOTH_INTEGRATION_RK4);
	//m_ClothViewer->GetScenario(0)->GetCloth()->SetUseRK4Integrator(true);
	/*m_ClothViewer->GetScenario(0)->GetCloth()->SetUseRK4Integrator(false);
	m_ClothViewer->GetScenario(1)->GetCloth()->SetUseRK4Integrator(true);
	m_ClothViewer->GetScenario(2)->GetCloth()->SetUseRK4Integrator(false);*/

	this->AddGameObject(m_ClothViewer);



//	//Matrix4 transform = Matrix4::Translation(Vector3(0, 1.0f, -0.0f));// *Matrix4::Rotation(90.0f, Vector3(1, 0, 0));
//	Matrix4 transform = Matrix4::Translation(Vector3(0, 1.0f, -5.0f));
//	//transform.ToIdentity();
//
//	Mat44 tmp;
//	for (int i = 0; i < 4; ++i)
//	{
//		for (int j = 0; j < 4; ++j)
//		{
//			tmp(i, j) = transform[j * 4 + i];
//		}
//	}
//
//#if (MODE==3)
//	ClothDesignEntity<6> cloth_design;
//#else
//	ClothDesignEntity<3> cloth_design;
//#endif
//	cloth_design.GenerateDefaultClothVertices(n_divisions, 1.0f, 1.0f, &tmp); //1m x 1m (subdivided 33 times along x/y
//	cloth_design.GenerateDefaultClothTriangles(n_divisions);
//
//	//Set Static Vertices
//	for (Vertex& v : cloth_design.Vertices())
//	{
//		v.flags = 0;
//	}
//	cloth_design.Vertices()[0].flags = VERTEXFLAGS_IS_STATIC;
//	cloth_design.Vertices()[n_divisions - 1].flags = VERTEXFLAGS_IS_STATIC;
//
//	cloth_design.Vertices()[(n_divisions - 1) * n_divisions + 0].flags = VERTEXFLAGS_IS_STATIC;
//	cloth_design.Vertices()[(n_divisions - 1) * n_divisions + n_divisions - 1].flags = VERTEXFLAGS_IS_STATIC;
//
//	for (uint i = 0; i < n_divisions; ++i)
//	{
//	cloth_design.Vertices()[i].flags = VERTEXFLAGS_IS_STATIC;
//	cloth_design.Vertices()[(n_divisions - 1) * n_divisions + i].flags = VERTEXFLAGS_IS_STATIC;
//	}
//
//
//	CLOTH_SIM_TYPE* cloth = new CLOTH_SIM_TYPE();
//	{
//		Ellipsoid e;
//		e.Transform = Matrix4::Translation(Vector3(0.0f, 0.5f, 0.0f));
//		e.radius = 0.5f;
//		e.invTransform = Matrix4::Inverse(e.Transform);
//		cloth->m_Ellipsoids.push_back(e);
//
//		GameObject* e_obj = BuildSphereObject("", Vector3(0.0f, 0.5f, 0.0f), 0.45f, 0.0f);
//		e_obj->SetColour(Vector4(1.0f, 0.8f, 1.0f, 0.5f));
//		this->AddGameObject(e_obj);
//	}
//	cloth->simulation_OnClothDesignChanged(&cloth_design);
//
//	if (!m_ClothRenderer)
//	{
//		m_ClothRenderer = new ClothRenderObject(cloth);
//		this->AddGameObject(m_ClothRenderer);
//	}
}

bool MyScene::InitialiseGL()
{
	InitialiseCloth();

	m_Camera->SetPosition(Vector3(-3.0f, 10.0f, 10.0f));
	m_Camera->SetPitch(-20.f);


	//Create Ground
	GameObject* ground = BuildCuboidObject("Ground", Vector3(0.0f, -1.001f, 0.0f), Vector3(20.0f, 1.0f, 20.0f), 0.0f, Vector4(0.2f, 0.5f, 1.0f, 1.0f));

	this->AddGameObject(ground);


	LoadCameraData();

	return true;
}

void MyScene::RenderScene()
{
	Scene::RenderScene();

	m_Encoder->EncodeFrame();
}

bool m_OrthoMode = false;
float m_OrthoScale = 2.5f;
void MyScene::UpdateScene(float dt)
{
	Scene::UpdateScene(dt);

	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_O))
	{
		m_OrthoMode = !m_OrthoMode;

		float aspect = (float)width / (float)height;
		if (m_OrthoMode)
			projMatrix = Matrix4::Orthographic(0.01f, 1000.0f, m_OrthoScale * aspect, -m_OrthoScale * aspect, m_OrthoScale, -m_OrthoScale);
		else
			projMatrix = Matrix4::Perspective(0.01f, 1000.0f, aspect, 45.0f);
	}

	if (m_OrthoMode)
	{
		int wheel = Window::GetMouse()->GetWheelMovement();
		if (wheel != 0)
		{
			m_OrthoScale += wheel * 0.1f;
			float aspect = (float)width / (float)height;
			projMatrix = Matrix4::Orthographic(0.01f, 1000.0f, m_OrthoScale * aspect, -m_OrthoScale * aspect, m_OrthoScale, -m_OrthoScale);
		}
	}

	ClothScenario* scenario = m_ClothViewer->GetScenario();
	if (scenario != NULL)
	{
		ClothBase* cloth = scenario->GetCloth();

		float time_step = cloth->Timestep();
		int n_steps = (floor)(1.0f / 60.0f / time_step);

		float total_time = (m_ClothViewer->IsPaused()) ? 0.0f : cloth->Profiler().GetTimingMS(PROFILERID_CLOTH_TOTALTIME);
		if (m_HudVisibility < 3)
		{
			NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "No. Elements    : %d", cloth->m_NumTris);
			NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Time Step       : %fms  (~%dfps)", time_step, n_steps);
			NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Cloth Update    : %5.2fms", total_time);
			NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Scenario %d-%d of %d (Ctrl A/D to change)", m_ClothViewer->GetScenarioIndex() + 1, m_ClothViewer->GetScenarioIndex() + m_ClothViewer->GetNumViewedScenarios(), m_ClothViewer->GetNumScenarios(), m_ClothViewer->GetScenarioName().c_str());
		
			for (int i = 0; i < m_ClothViewer->GetNumViewedScenarios(); ++i)
			{
				int idx = (m_ClothViewer->GetScenarioIndex() + i) % m_ClothViewer->GetNumScenarios();
				NCLDebug::AddStatusEntry(Vector4(1.0f, 0.8f, 1.0f, 1.0f), "Scenario#%d: %s", idx + 1, m_ClothViewer->GetScenarioName(idx).c_str());
			}
		
		}

		if (!m_ClothViewer->IsPaused())
		{
			m_GraphObject->UpdateYScale(1.0f / 60.0f);
			m_ClothUpdateMs = total_time;

			float rotations_time = cloth->Profiler().GetTimingMS(PROFILERID_CLOTH_BUILDROTATIONANDNORMALS);
			float matricies_time = cloth->Profiler().GetTimingMS(PROFILERID_CLOTH_COMPUTEFORCES);
			float cg_time = cloth->Profiler().GetTimingMS(PROFILERID_CLOTH_SOLVER);
			float col_time = cloth->Profiler().GetTimingMS(PROFILERID_CLOTH_EXTERNALCOLLISIONS);

			float val = rotations_time;
			m_GraphObject->AddDataEntry(0, val);
			val += matricies_time; m_GraphObject->AddDataEntry(1, val);
			val += cg_time;  m_GraphObject->AddDataEntry(2, val);
			val += col_time; m_GraphObject->AddDataEntry(3, val);



			m_GraphObjectSolver->UpdateYScale(1.0f / 60.0f);
			float s_total_time = cg_time;
			float s_init_time = cloth->Profiler().GetTimingMS(PROFILERID_SOLVER_INITIALIZATION);
			float s_upper_time = cloth->Profiler().GetTimingMS(PROFILERID_SOLVER_UPPERAMATRIX);
			float s_lower_time = cloth->Profiler().GetTimingMS(PROFILERID_SOLVER_LOWERAMATRIX);

			float s_p_init = (s_init_time / s_total_time) * 100.0f;
			float s_p_upper = (s_upper_time / s_total_time) * 100.0f;
			float s_p_lower = (s_lower_time / s_total_time) * 100.0f;


			val = s_p_init;
			m_GraphObjectSolver->AddDataEntry(0, val);
			val += s_p_upper; m_GraphObjectSolver->AddDataEntry(1, val);
			val += s_p_lower;  m_GraphObjectSolver->AddDataEntry(2, val);


			//NCLDebug::Log(Vector3(1.0f, 1.0f, 1.0f), "Solver Time: %5.2fms\t Iterations: %5.2f", s_total_time, m_ClothRenderer->m_Cloth->m_Solver.GetAverageIterations());
			//NCLDebug::Log(Vector3(1.0f, 1.0f, 1.0f), "KMatrix Time: %5.2fms", matricies_time);
		}
	}
	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_G))
	{
		m_ClothViewer->TogglePaused();
		if (m_ClothViewer->IsPaused())
			m_Encoder->EndEncoding();
	}

	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_R))
	{
		m_ClothViewer->SingleTimeStep();
		m_Encoder->EndEncoding();
	}

	if (Window::GetKeyboard()->KeyDown(KEYBOARD_CONTROL))
	{
		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_A))
		{
			m_ClothViewer->GotoPreviousScenario();
		}
		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_D))
		{
			m_ClothViewer->GotoNextScenario();
		}

		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_P))
		{

			m_ClothViewer->TogglePaused();

			if (!m_ClothViewer->IsPaused())
			{
				stringstream ss;
				ss << "Videos\\screenCapture_";
				ss << m_EncoderVideoNum;
				ss << ".mp4";
				m_Encoder->BeginEncoding(width, height, ss.str());
				m_EncoderVideoNum++;

				ofstream myfile;
				myfile.open("videos\\index.txt");
				if (myfile.is_open())
				{
					myfile << m_EncoderVideoNum;
					myfile.close();
				}
			}
			else
				m_Encoder->EndEncoding();
		}

		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_C))
		{
			LoadCameraData();
			NCLDebug::Log(Vector3(0.0f, 0.0f, 0.0f), "Loaded Camera Location");
		}
		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_S))
		{
			SaveCameraData();
			NCLDebug::Log(Vector3(1.0f, 0.0f, 0.0f), "Saved Camera Location");
		}

		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_H))
		{
			m_HudVisibility = (m_HudVisibility + 1) % MAX_HUD_VISIBILITY_TYPES;
			m_GraphObject->SetVisibility(m_HudVisibility < 2);
			m_GraphObjectSolver->SetVisibility(m_HudVisibility < 1);
		}
	}
}

void MyScene::SaveCameraData()
{
	ofstream myfile;
	myfile.open("screen_info.txt");
	if (myfile.is_open())
	{
		Vector3 pos = m_Camera->GetPosition();
		float pitch = m_Camera->GetPitch();
		float yaw   = m_Camera->GetYaw();

		int width, height;
		Window::GetWindow().GetWindowSize(&width, &height);

		myfile << width << "\n";
		myfile << height << "\n";
		myfile << pos.x << "\n";
		myfile << pos.y << "\n";
		myfile << pos.z << "\n";
		myfile << pitch << "\n";
		myfile << yaw;

		myfile.close();
	}
}

void MyScene::LoadCameraData()
{
	ifstream myfile;
	myfile.open("screen_info.txt");
	if (myfile.is_open())
	{
		Vector3 pos;
		float pitch, yaw;

		int width, height;
		myfile >> width;
		myfile >> height;

		Window::GetWindow().SetWindowSize(width, height);

		myfile >> pos.x;
		myfile >> pos.y;
		myfile >> pos.z;
		myfile >> pitch;
		myfile >> yaw;

		m_Camera->SetPosition(pos);
		m_Camera->SetPitch(pitch);
		m_Camera->SetYaw(yaw);
		myfile.close();
	}
}

GameObject* MyScene::BuildSphereObject(const std::string& name, const Vector3& pos, float radius, float invmass, const Vector4& color)
{
	//Create Ground
	SimpleMeshObject* sphere = new SimpleMeshObject(name);
	sphere->SetMesh(CommonMeshes::Sphere(), false);
	sphere->SetLocalTransform(Matrix4::Scale(Vector3(radius, radius, radius)));
	sphere->SetColour(color);
	sphere->SetBoundingRadius(radius);

	sphere->Physics()->SetPosition(pos);
	sphere->Physics()->SetInverseMass(invmass);

	//Build Inverse Inertia Matrix of a sphere
	Matrix3 inertia;
	{
		float i = 2.5f * invmass * radius * radius;

		inertia._11 = i;
		inertia._22 = i;
		inertia._33 = i;
	}
	sphere->Physics()->SetInverseInertia(inertia);

	return sphere;
}

GameObject* MyScene::BuildCuboidObject(const std::string& name, const Vector3& pos, const Vector3& halfdims, float invmass, const Vector4& color)
{
	SimpleMeshObject* cuboid = new SimpleMeshObject(name);
	cuboid->SetMesh(CommonMeshes::Cube(), false);
	cuboid->SetLocalTransform(Matrix4::Scale(halfdims));
	cuboid->SetColour(color);
	cuboid->SetBoundingRadius(80.0f * 80.f);

	GLuint tex = SOIL_load_OGL_texture(TEXTUREDIR"white.bmp", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);
	glBindTexture(GL_TEXTURE_2D, tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST); //No linear interpolation to get crisp checkerboard no matter the scalling
	glBindTexture(GL_TEXTURE_2D, 0);
	cuboid->SetTexture(tex, true);

	cuboid->Physics()->SetPosition(pos);

	cuboid->Physics()->SetInverseMass(invmass);
	//Build Inverse Inertia Matrix of a cuboid
	Matrix3 inertia;
	{
		Vector3 dimsSq = (halfdims + halfdims);
		dimsSq = dimsSq * dimsSq;

		inertia._11 = 12.f * invmass * 1.f / (dimsSq.y + dimsSq.z);
		inertia._22 = 12.f * invmass * 1.f / (dimsSq.x + dimsSq.z);
		inertia._33 = 12.f * invmass * 1.f / (dimsSq.x + dimsSq.y);
	}
	cuboid->Physics()->SetInverseInertia(inertia);

	return cuboid;
}
