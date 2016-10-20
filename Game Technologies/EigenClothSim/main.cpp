#include <nclgl\Window.h>
#include "MyScene.h"
#include <ncltech\PhysicsEngine.h>
#include <ncltech\NCLDebug.h>
#include <iomanip>

MyScene* scene = NULL;

int Quit(bool pause = false, const string &reason = "") {
	if (scene)
	{
		delete scene;
		scene = NULL;
	}

	Window::Destroy();

	if (pause) {
		std::cout << reason << std::endl;
		system("PAUSE");
	}

	return 0;
}

int main()
{
	//-------------------
	//--- MAIN ENGINE ---
	//-------------------

	//Initialise the Window
	if (!Window::Initialise("Cloth Simultaion", 1366, 768, false))
	{
		return Quit(true, "Window failed to initialise!");
	}

	//Initialise the PhysicsEngine
	PhysicsEngine::Instance();

	//Initialise the Scene
	scene = new MyScene(Window::GetWindow());
	if (!scene->HasInitialised())
	{
		return Quit(true, "Renderer failed to initialise!");
	}

	GameTimer engine_timer;

	const float title_update_seconds = 0.33f;
	float time_accum = 0.0f, time_accum_cloth = 0.0f;
	int frames = 0;
	//Create main game-loop
	while (Window::GetWindow().UpdateWindow() && !Window::GetKeyboard()->KeyDown(KEYBOARD_ESCAPE)){
		float dt = Window::GetWindow().GetTimer()->GetTimedMS() * 0.001f;	//How many milliseconds since last update?
	
		//NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Cloth Simulation: %s (Press P to toggle)", PhysicsEngine::Instance()->IsPaused() ? "Paused" : "Enabled");
		//NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "----------------------------------------");


		if (Window::GetKeyboard()->KeyTriggered(KEYBOARD_P))
			PhysicsEngine::Instance()->SetPaused(!PhysicsEngine::Instance()->IsPaused());


		engine_timer.GetTimedMS();

		PhysicsEngine::Instance()->Update(dt);
		float physics_ms = engine_timer.GetTimedMS();

		scene->UpdateScene(dt);
		float update_ms = engine_timer.GetTimedMS();

		//NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Scene Update    : %5.2fms", update_ms);

		scene->RenderScene();

		time_accum += dt;
		time_accum_cloth += scene->m_ClothUpdateMs;
		frames++;

		if (time_accum > title_update_seconds)
		{
			time_accum /= (float)frames;
			time_accum_cloth /= (float)frames;

			float fps = 1.0f / time_accum;
			stringstream ss;
			ss.setf(ios::fixed);
			ss << "Cloth Simultaion " << std::setprecision(2) << fps << "fps    - Cloth: ~" << time_accum_cloth << "ms";
			Window::GetWindow().SetWindowTitle(ss.str().c_str());

			time_accum = 0.0f;
			time_accum_cloth = 0.0f;
			frames = 0;
		}

	}

	//Cleanup
	return Quit();
}