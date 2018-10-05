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
using namespace glm;

float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

// time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

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
	particle1.rotate((GLfloat)M_PI_2, glm::vec3(1.0f, 0.0f, 0.0f));
	particle1.setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));


	//Create particle via class:

	const int particlenum = 10;
	std::vector<Particle> allPart;
	for (unsigned int i = 0; i < particlenum; ++i) {
		allPart.push_back(Particle::Particle());
		allPart[i].getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));
	}

	//Particle *particle2 = new Particle();
	//particle2->getVel();

	// create demo objects (a cube and a sphere)
	Mesh sphere = Mesh::Mesh("resources/models/sphere.obj");
	sphere.translate(glm::vec3(-1.0f, 1.0f, 0.0f));
	sphere.setShader(lambert);
	Mesh cube = Mesh::Mesh("resources/models/cube.obj");
	cube.translate(glm::vec3(1.0f, .5f, 0.0f));
	cube.setShader(lambert);

	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	///////////////////////////////////////////////////////Forces

	glm::vec3 Ftotal; //Force applied to the particle
	glm::vec3 Fg;



	//TASK 1 Variables
	glm::vec3 cubecorner = glm::vec3(-5.0f, 0.0f, -5.0f);
	glm::vec3 d = glm::vec3(10.0f);


	float density = 1.225f;
	float coefficient = 1.05f;
	vec3 e;
	vec3 area = vec3(0.1f, 0.1f, 0.0f);
	vec3 Fa;
	vec3 absoluteu;
	glm::vec3 g = glm::vec3(0.0f, -9.8f, 0.0f);


	glm::vec3 ParticleVelocity[particlenum]; //initial; velocity for multiple particles
	for (unsigned int i = 0; i < particlenum; ++i) {
		allPart[i].setVel(vec3(RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 8.0f)));
		allPart[i].setPos(vec3(RandomFloat(0.0f, 5.0f), RandomFloat(0.0f, 8.0f), RandomFloat(0.0f, 8.0f)));
	}




	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{
		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
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
		for (unsigned int i = 0; i < particlenum; i++) {

			Fg = allPart[i].getMass() * g;
			absoluteu = abs(allPart[i].getVel());
			e = -allPart[i].getVel() / absoluteu;
			Fa = 0.5*density*absoluteu*absoluteu*coefficient*area*e;
			Ftotal = Fg + Fa ;
			allPart[i].setAcc(Ftotal/allPart[i].getMass());


			allPart[i].setVel(allPart[i].getVel()+deltaTime*allPart[i].getAcc());
			//


			for (unsigned int j = 0; j < 3; j++) {
				if (allPart[i].getPos()[j] < cubecorner[j]) {
					//rIn[i] = particle1.getPos()[i];
					allPart[i].setVel(j, allPart[i].getVel()[j]*-0.9f);
					/*N = m * g;
					Ff = N * frictioncoef;*/  //Was trying to add friction here			
					allPart[i].setPos(j, cubecorner[j]);

				}
				else if (allPart[i].getPos()[j] > cubecorner[j] + d[j]) {
					//rIn[i] = particle1.getPos()[i];
					allPart[i].setVel(j, allPart[i].getVel()[j]*-0.9f);
					/*N = m * g;
					Ff = N * frictioncoef;*/  //Was trying to add friction here
					allPart[i].setPos(j, cubecorner[j]+d[j]);
				}
			}

			allPart[i].translate(deltaTime*allPart[i].getVel());


		}



		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		// draw particles
		//app.draw(particle1);	////////////////////////////////////////////////////put partile here
		for (unsigned int i = 0; i < particlenum; ++i) {
			app.draw(allPart[i].getMesh());
		}
		// draw demo objects
		//app.draw(cube);
		//app.draw(sphere);

		app.display();
	}

	app.terminate();

	return EXIT_SUCCESS;
}

