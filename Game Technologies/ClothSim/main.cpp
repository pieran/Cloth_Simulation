#include <nclgl\Window.h>
#include "ClothScene.h"
#include <ncltech\PhysicsEngine.h>
#include <ncltech\NCLDebug.h>

#pragma comment(lib, "nclgl.lib")
#pragma comment(lib, "ncltech.lib")


Scene* scene = NULL;

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
	if (!Window::Initialise("Yusef Cloth Simulation", 1280, 800, false))
	{
		return Quit(true, "Window failed to initialise!");
	}

	//Initialise the Scene
	scene = new ClothScene(Window::GetWindow());
	if (!scene->HasInitialised())
	{
		return Quit(true, "Renderer failed to initialise!");
	}

	GameTimer engine_timer;

	//Create main game-loop
	while (Window::GetWindow().UpdateWindow() && !Window::GetKeyboard()->KeyDown(KEYBOARD_ESCAPE)){
		float dt = Window::GetWindow().GetTimer()->GetTimedMS() * 0.001f;	//How many milliseconds since last update?

		engine_timer.GetTimedMS();
		float physics_ms = engine_timer.GetTimedMS();

		scene->UpdateScene(dt);
		float update_ms = engine_timer.GetTimedMS();

		NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Cloth Simulation");
		NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "--------------------------------");
		NCLDebug::AddStatusEntry(Vector4(1.0f, 1.0f, 1.0f, 1.0f), "Cloth Update  : %5.2fms", update_ms);

		scene->RenderScene();
	}

	//Cleanup
	return Quit();
}