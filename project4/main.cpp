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
#include "Body.h"
#include "Particle.h"
#include "Force.h"
using namespace glm;

// time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

Gravity g = Gravity(glm::vec3(0.0f, -9.8f, 0.0f));

//Drag d = Drag();


//Method tyhat creates random variables between 2 values
float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

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


	//// create particle
	//Particle particle1 = Particle();
	////scale it down (x.1), translate it up by 2.5 and rotate it by 90 degrees around the x axis
	//particle1.translate(glm::vec3(0.0f, 2.5f, 0.0f));
	////particle1.scale(glm::vec3(.1f, .1f, .1f));
	////particle1.rotate((GLfloat) M_PI_2, glm::vec3(1.0f, 0.0f, 0.0f));
	//particle1.getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));

	//Create particle via class:

	const int particlenum = 2;
	std::vector<Particle> allPart;
	Shader blue = Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag");

	for (unsigned int i = 0; i < particlenum; ++i) {
		allPart.push_back(Particle::Particle());
		allPart[i].getMesh().setShader(blue);
	}

	for (unsigned int i = 0; i < particlenum; ++i) {
		//allPart[i].setVel(vec3(RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 8.0f)));
		allPart[i].setVel(vec3(0.0f));
		//allPart[i].setPos(vec3(RandomFloat(0.0f, 5.0f), RandomFloat(5.0f, 9.0f), RandomFloat(0.0f, 8.0f)));
	}
	allPart[0].setPos(vec3(0.0f, 10.0f, 0.0f));
	allPart[1].setPos(vec3(0.0f, 8.0f, 0.0f));

	// create demo objects (a cube and a sphere)
	Mesh sphere = Mesh::Mesh("resources/models/sphere.obj");
	sphere.translate(glm::vec3(-1.0f, 1.0f, 0.0f));
	sphere.setShader(lambert);
	Mesh cube = Mesh::Mesh("resources/models/cube.obj");
	cube.translate(glm::vec3(1.0f, .5f, 0.0f));
	cube.setShader(lambert);

	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();


	//Timestep variables
	float t = 0.0f;
	const float dt = 0.01f;
	float currentTime = (GLfloat)glfwGetTime();
	float accumulator = 0.0f;

	//Cube Variables
	glm::vec3 cubecorner = glm::vec3(-5.0f, 0.0f, -5.0f);
	glm::vec3 d = glm::vec3(10.0f);





	//Hooke force

	float ks = 10.0f; //spring stiffness 
	float rest = 4.0f;//spring rest length
	float kd = 10.0f;//damping coefficient

	Hooke* h = new Hooke(&(allPart[0]), &(allPart[1]), ks, kd, rest);
	
	//h->setks(ks);
	//h->setrest(rest);
	//h->setkd(kd);
	

	//Adding the forces applied to the particle
	//for (unsigned int i = 0; i < particlenum; i++) {
		//allPart[1].addForce(&g);
		allPart[1].addForce(h);
	//}


	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{
		//// Set frame time
		//GLfloat currentFrame = (GLfloat)  glfwGetTime() - firstFrame;
		//// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		//currentFrame *= 1.5f;
		//deltaTime = currentFrame - lastFrame;
		//lastFrame = currentFrame;
		float newTime = (GLfloat)glfwGetTime();
		float frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(dt);

		while (accumulator >= dt) {
			/*
			**	SIMULATION
			*/
			for (unsigned int i = 0; i < particlenum; i++) {

				allPart[i].setAcc(allPart[i].applyForces(allPart[i].getPos(), allPart[i].getVel(), t, dt));

				/*for (unsigned int j = 0; j < 3; j++) {
					if (allPart[i].getPos()[j] < cubecorner[j]) {
						allPart[i].setVel(j, allPart[i].getVel()[j] * -0.5f);
						allPart[i].setPos(j, cubecorner[j]);
					}
					else if (allPart[i].getPos()[j] > cubecorner[j] + d[j]) {
						allPart[i].setVel(j, allPart[i].getVel()[j] * -0.5f);
						allPart[i].setPos(j, cubecorner[j] + d[j]);
					}

				}*/

				allPart[i].setVel(allPart[i].getVel() + dt * allPart[i].getAcc());
				allPart[i].translate(dt*allPart[i].getVel());
			}



			accumulator -= dt;
			t += dt;
		}

		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		// draw particles
		//app.draw(particle1.getMesh());
		for (unsigned int i = 0; i < particlenum; i++) {
			app.draw(allPart[i].getMesh());
		}

		// draw demo objects
		app.draw(cube);
		app.draw(sphere);

		app.display();
	}

	app.terminate();

	return EXIT_SUCCESS;
}

