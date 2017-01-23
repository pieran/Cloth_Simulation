
#include "MyScene.h"

#include <nclgl/OBJMesh.h>

#include <ncltech\SimpleMeshObject.h>
#include <ncltech\PhysicsEngine.h>
#include <ncltech\CommonMeshes.h>
#include <ncltech\DistanceConstraint.h>


MyScene::MyScene(Window& window) : Scene(window), m_Sim(NULL), m_SimTimestep(1.0f / 60.0f), m_SimPaused(true)
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
	const Matrix4 transform = Matrix4();// Matrix4::Scale(Vector3(0.75f, 0.75f, 1.0f));
	const uint visual_subdivisions = 4;
	const uint subdivisions = visual_subdivisions * 2 + 1;

	uint num_verts = subdivisions * subdivisions;
	FEVertDescriptor* vertices = new FEVertDescriptor[num_verts];

	uint num_quads = (subdivisions - 1) * (subdivisions - 1) / 4;
	uint num_tris = num_quads * 2;
	FETriangle* tris = new FETriangle[num_tris];

	float width = 1.f;
	float height = 1.f;

	const float divisor_scalar = 1.f / (subdivisions - 1);

	unsigned int count = 0;
	for (unsigned int y = 0; y < subdivisions; ++y)
	{
		for (unsigned int x = 0; x < subdivisions; ++x)
		{
			Vector3 pos = transform * Vector3(
				((float(x) * divisor_scalar) - 0.5f) *  width,
				((float(y) * divisor_scalar) - 0.5f)  * -height,
				0.0f);

			vertices[count].pos = pos;
			vertices[count].ipos = Vector3(
				((float(x) * divisor_scalar) - 0.5f) *  width,
				((float(y) * divisor_scalar) - 0.5f)  * -height,
				0.0f);
			vertices[count].tCoord = Vector2(
				(float(x) * divisor_scalar),
				(float(y) * divisor_scalar)
				);
			vertices[count].isStatic = (y == 0 || y == subdivisions - 1);// && (x == 0 || x == subdivisions - 1);
			count++;
		}
	}



	uint tcount = 0;
	uint qcount = 0;
//	num_tris /= 2;
	for (uint x = 2; x < subdivisions; x += 2)
	{
		for (uint y = 2; y < subdivisions; y += 2)
		{
			uint a = subdivisions*(y - 2) + x - 2;
			uint b = subdivisions*(y - 2) + x;
			uint c = subdivisions*(y)+x;
			uint d = subdivisions*(y)+x - 2;

			uint ab = subdivisions*(y - 2) + x - 1;
			uint dc = subdivisions*(y)+x - 1;
			uint ad = subdivisions*(y - 1) + x - 2;
			uint bc = subdivisions*(y - 1) + x;

			uint mid = subdivisions*(y - 1) + x - 1;



			//Alternate triangle orientation to allow FEM to better represent the forces
			if (false)//((y ^ x) & 0x2) == 0)
			{
				//tris[tcount++] = FETriangle(a, d, c, ad, dc, mid);
				//tris[tcount++] = FETriangle(a, c, b, mid, bc, ab);

				tris[tcount++] = FETriangle(a, d, b, ad, mid, ab);
				tris[tcount++] = FETriangle(d, c, b, dc, bc, mid);
			}
			else
			{
				tris[tcount++] = FETriangle( d, b, a, mid, ab, ad);
				tris[tcount++] = FETriangle( d, c, b, dc, bc, mid);
			}
		}
	}




	if (m_Sim)
		delete m_Sim;

	m_Sim = new FEM6Noded();
	m_Sim->simulation_OnClothDesignChanged(num_tris, tris, num_verts, vertices);

	delete[] tris;
	delete[] vertices;
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
	m_Sim->Render_DrawingToVisualiser();
	Scene::RenderScene();

	m_Encoder->EncodeFrame();
}

bool m_OrthoMode = false;
float m_OrthoScale = 2.5f;
void MyScene::UpdateScene(float dt)
{
	Scene::UpdateScene(dt);

	if (!m_SimPaused)
		m_Sim->Simulation_StepSimulation(m_SimTimestep);

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

	{
		float time_step = m_SimTimestep;
		int n_steps = (floor)(1.0f / 60.0f / time_step);

		float total_time = (m_SimPaused) ? 0.0f : m_Sim->m_ProfilingTotalTime.GetTimedMilliSeconds();
		if (m_HudVisibility < 3)
		{
			NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "No. Elements    : %d", m_Sim->m_NumTriangles);
			NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Time Step       : %fms  (~%dfps)", time_step, n_steps);
			NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Cloth Update    : %5.2fms", total_time);

		}

		if (!m_SimPaused)
		{
			m_GraphObject->UpdateYScale(1.0f / 60.0f);
			m_ClothUpdateMs = total_time;

			float rotations_time = m_Sim->m_ProfilingRotationNormals.GetTimedMilliSeconds();
			float matricies_time = m_Sim->m_ProfilingBuildMatricies.GetTimedMilliSeconds();
			float cg_time = m_Sim->m_ProfilingConjugateGradient.GetTimedMilliSeconds();
			float col_time = m_Sim->m_ProfilingExternalCollisions.GetTimedMilliSeconds();

			float val = rotations_time;
			m_GraphObject->AddDataEntry(0, val);
			val += matricies_time; m_GraphObject->AddDataEntry(1, val);
			val += cg_time;  m_GraphObject->AddDataEntry(2, val);
			val += col_time; m_GraphObject->AddDataEntry(3, val);



			m_GraphObjectSolver->UpdateYScale(1.0f / 60.0f);
			float s_total_time = cg_time;
			float s_init_time = m_Sim->m_Solver.m_ProfilingInitialization.GetTimedMilliSeconds();
			float s_upper_time = m_Sim->m_Solver.m_ProfilingUpper.GetTimedMilliSeconds();
			float s_lower_time = m_Sim->m_Solver.m_ProfilingLower.GetTimedMilliSeconds();

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
		m_SimPaused = !m_SimPaused;
		if (m_SimPaused)
			m_Encoder->EndEncoding();
	}

	if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_R))
	{
		m_Sim->Simulation_StepSimulation(m_Sim->m_TimeStep);
		m_Encoder->EndEncoding();
	}

	if (Window::GetKeyboard()->KeyDown(KEYBOARD_CONTROL))
	{
		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_R))
		{
			this->InitialiseCloth();
		}


		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_P))
		{

			m_SimPaused = !m_SimPaused;

			if (!m_SimPaused)
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
		float yaw = m_Camera->GetYaw();

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
