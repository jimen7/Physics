#pragma once
// Math constants
#define _USE_MATH_DEFINES
#include <cmath>  
#include <random>

// Std. Includes
#include <string>
#include <time.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include "glm/ext.hpp"

// Other Libs
#include "SOIL2/SOIL2.h"

// project includes
#include "Application.h"
#include "Shader.h"
#include "Mesh.h"
using namespace glm;


// time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;


float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}
/*
struct particle {                                          //////////////////////////////////////////////////////////////////////////QUESTION 6 STUCT
	vec3 velocity;
	GLfloat firstFrame;
	GLfloat currentFrame;
	GLfloat lastFrame;
	Mesh mesh; //mesh
	vec3 a; //acceleration
	vec3 Ftotal; //Force applied to the particle
	const float m = 1.0f;
	vec3 Fg; //Force applied by gravity
	vec3 rFin; //Displacement Final
	vec3 rIn; //Displacement Initial

	particle::particle() {
		rIn = glm::vec3(RandomFloat(0.0f, 5.0f), RandomFloat(0.0f, 5.0f), RandomFloat(0.0f, 5.0f));
		velocity = glm::vec3(RandomFloat(0.0f, 5.0f), RandomFloat(0.0f, 5.0f), RandomFloat(0.0f, 5.0f));
		firstFrame = (GLfloat)glfwGetTime();
		mesh = Mesh::Mesh(Mesh::QUAD);
		a = glm::vec3(0.0f, -9.8f, 0.0f);
		currentFrame = 0.0f;
		lastFrame = 0.0f;
		rIn = glm::vec3(0.0f, 5.0f, 0.0f);

		//scale it down (x.1), translate it up by 2.5 and rotate it by 90 degrees around the x axis
		mesh.translate(glm::vec3(0.0f, 2.5f, 0.0f));
		mesh.scale(glm::vec3(.1f, .1f, .1f));
		mesh.rotate((GLfloat)M_PI_2, glm::vec3(1.0f, 0.0f, 0.0f));
		// allocate shader
		mesh.setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));
	}
};
*/


// main function
int main()
{
	// create application
	Application app = Application::Application();
	app.initRender();
	Application::camera.setCameraPosition(glm::vec3(0.0f, 5.0f, 20.0f));
			
	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	Shader lambert = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	plane.setShader(lambert);


	// create particle
	Mesh particle1 = Mesh::Mesh(Mesh::QUAD);
	//scale it down (x.1), translate it up by 2.5 and rotate it by 90 degrees around the x axis
	particle1.translate(glm::vec3(0.0f, 2.5f, 0.0f));
	particle1.scale(glm::vec3(.1f, .1f, .1f));
	particle1.rotate((GLfloat) M_PI_2, glm::vec3(1.0f, 0.0f, 0.0f));
	particle1.setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));
	
	/*std::vector<particle> allPart;
const int MAXPARTICLES = 10;

for (unsigned int i = 0; i < MAXPARTICLES; ++i) {
	allPart.push_back(particle());
}*/


	// create demo objects (a cube and a sphere)
	Mesh sphere = Mesh::Mesh("resources/models/sphere.obj");
	sphere.translate(glm::vec3(-1.0f, 1.0f, 0.0f));
	sphere.setShader(lambert);
	Mesh cube = Mesh::Mesh("resources/models/cube.obj");
	cube.translate(glm::vec3(1.0f, .5f, 0.0f));
	cube.setShader(lambert);

	// time
	GLfloat firstFrame = (GLfloat) glfwGetTime();


	//Declaring particle 
	glm::vec3 a; //acceleration
	//glm::vec3 inpos = glm::vec3(0.0f, 5.0f, 0.0f); //initial position
	glm::vec3 u = glm::vec3(5.0f, 20.0f, 3.0f); //initisal velocity
	glm::vec3 Ftotal; //Force applied to the particle
	//glm::vec3 uFin; //Final velocity
	static float m = 1.0f;
	glm::vec3 Fg;
	glm::vec3 rFin; //Displacement Final
	glm::vec3 rIn = glm::vec3(0.0f, 5.0f, 0.0f); //Displacement Initial

	//TASK 3 VARIABLES
	glm::vec3 cubecorner = glm::vec3(-5.0f, 0.0f, -5.0f);
	glm::vec3 d = glm::vec3(10.0f);

	//Task 4 variables

	float density = 1.225f;
	float coefficient = 1.05f;
	vec3 e;
	vec3 area = vec3(0.1f, 0.1f, 0.0f);
	vec3 Fa;
	vec3 absoluteu;
	glm::vec3 g = glm::vec3(0.0f, -9.8f, 0.0f);
	
	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{
		// Set frame time
		GLfloat currentFrame = (GLfloat)  glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);


		/*
		**	SIMULATION
		*/

		//Task 1

		u = u + deltaTime * g;
		/*Fg = m * g;
		absoluteu = abs(u);
		e = -u / absoluteu;
		Fa = 0.5*density*absoluteu*absoluteu*coefficient*area*e;
		Ftotal = Fg - Fa;
		a = Ftotal / m;*/

		//rIn = rFin;
		//rFin = rIn + deltaTime * u;



		//Task 2
	//	if (particle1.getPos().x < (cubecorner-d).x || particle1.getPos().y < (cubecorner - d).y || particle1.getPos().z < (cubecorner + d).z) {

		//	firstFrame = (GLfloat)glfwGetTime();
	//	}

		for (unsigned int i = 0; i < 3; ++i) {
			if (particle1.getPos()[i] < cubecorner[i]) {
				//rIn[i] = particle1.getPos()[i];
				u[i] = -u[i] * 0.9f;
				particle1.setPos(i, cubecorner[i]);

			}
			else if (particle1.getPos()[i] > cubecorner[i] + d[i]) {
				//rIn[i] = particle1.getPos()[i];
				u[i] = -u[i] * 0.9f;
				particle1.setPos(i, cubecorner[i] + d[i]);
			}
		}
		particle1.translate(deltaTime * u);
		


		/*
		**	RENDER 
		*/		
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		// draw particles
		app.draw(particle1);	

		// draw demo objects
		app.draw(cube);
		app.draw(sphere);

		app.display();
	}

	app.terminate();

	return EXIT_SUCCESS;
}

