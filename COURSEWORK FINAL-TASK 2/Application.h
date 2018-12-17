#pragma once



// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// project includes
#include "Camera.h"
#include "Mesh.h"


class Application
{

public:
	Application();
	~Application();

	// static variables for callback function
	// window
	static const GLuint WIDTH = 800;
	static const GLuint HEIGHT = 600;
	static int SCREEN_WIDTH, SCREEN_HEIGHT;

	// Camera
	static Camera camera;
	static double lastX;
	static double lastY;

	static bool firstMouse; 
	static bool keys[1024];


	// get and set methods
	GLFWwindow* getWindow() { return m_window; }

	// Function prototypes
	
	//static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mode);
	//void scrollCallback(GLFWwindow *window, double xOffset, double yOffset);
	//void mouseCallback(GLFWwindow *window, double xPos, double yPos);
	void doMovement(GLfloat deltaTime);
	int initRender();

	// other functions
	void clear();
	void draw(const Mesh &mesh);
	void display();
	void terminate(){ glfwTerminate(); }


	void Application::showFPS() {

		static double previousSeconds = 0.0;

		static int frameCount;

		double elapsedSeconds;

		double currentSeconds = glfwGetTime();



		elapsedSeconds = currentSeconds - previousSeconds;



		// limit FPS refresh rate to 4 times per second

		if (elapsedSeconds > 0.25) {

			previousSeconds = currentSeconds;

			double FPS = (double)frameCount / elapsedSeconds;

			double msPerFrame = 1000.0 / FPS;



			std::ostringstream outs;

			outs.precision(3);

			outs << std::fixed

				<< "Physics Simulation  "

				<< "FPS: " << FPS << "  "

				<< "Frame time: " << msPerFrame << "ms" << std::endl;

			glfwSetWindowTitle(m_window, outs.str().c_str());



			frameCount = 0;

		}

		frameCount++;

	}

private:

	// view and projection matrices
	glm::mat4 m_view = glm::mat4(1.0f);
	glm::mat4 m_projection = glm::mat4(1.0f);

	// window
	GLFWwindow* m_window = NULL;

};

