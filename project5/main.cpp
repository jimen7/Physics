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

float ks = 10.0f; //spring stiffness 
float rest = 0.5f;//spring rest length
float kd = 5.0f;//damping coefficient


const int particlenum = 10;//particle num per vertex
const int vertnum = 10;
//Drag d = Drag();


//Method tyhat creates random variables between 2 values
float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}


void seflagPos(std::vector<std::vector<Particle>> &pvert) {
	vec3 inpos = vec3(-5.0f, 10.0f, 0.0f);
	float gap = 1.0f;
	for (unsigned int x = 0; x <10 ; ++x) {
		for (unsigned int y = 0; y < 10; ++y) {
			pvert[x][y].setPos(inpos + vec3(x*gap,0.0f,y*gap));
		}		
	}
}

void addHookeForces(std::vector<std::vector<Particle>> &pvert) {
	for (unsigned int j = 0; j < vertnum; j++) {
		for (int i = 0; i < particlenum-1; i++)
		{
			if (i != 0 && j != 0 && i != 9 && j != 9) {

				Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
				Hooke* left = new Hooke(&(pvert[j][i]), &(pvert[j][i - 1]), ks, kd, rest);
				Hooke* bottom = new Hooke(&(pvert[j][i]), &(pvert[j + 1][i]), ks, kd, rest);
				Hooke* right = new Hooke(&(pvert[j][i]), &(pvert[j][i + 1]), ks, kd, rest);

				pvert[j][i].addForce(&g);

				pvert[j][i].addForce(top);
				pvert[j][i].addForce(left);
				pvert[j][i].addForce(bottom);
				pvert[j][i].addForce(right);

				vec3 topnormal = (cross((pvert[i-1][j-1].getPos()- pvert[i][j].getPos()),(pvert[i-1][j].getPos()- pvert[i][j].getPos())))/glm::length(cross((pvert[i - 1][j - 1].getPos() - pvert[i][j].getPos()), (pvert[i - 1][j].getPos() - pvert[i][j].getPos())));
				float toparea = 0.5f * glm::length(cross((pvert[i - 1][j - 1].getPos() - pvert[i][j].getPos()), (pvert[i - 1][j].getPos() - pvert[i][j].getPos())));
			}
			


		}
	}
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


	Shader blue = Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag");
	//Create particles via class:
	//std::vector<Particle> allPart;
	//for (unsigned int i = 0; i < particlenum; ++i) {
	//	allPart.push_back(Particle::Particle());
	//	allPart[i].getMesh().setShader(blue);
	//}

	//Create vertice of vertices
	std::vector<std::vector<Particle>> vecvec(vertnum);
	for (unsigned int i = 0; i < particlenum; ++i) {
		for (unsigned int j = 0; j < vertnum; ++j) {
			vecvec[j].push_back(Particle::Particle());
			vecvec[j][i].getMesh().setShader(blue);
		}
	}



	//Set particles initial velocity
	for (unsigned int i = 0; i < particlenum; ++i) {
		for (unsigned int j = 0; j < vertnum; ++j) {
			//vecvec[j][i].setVel(vec3(RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 8.0f)));
			vecvec[j][i].setVel(vec3(0.0f));
			//vecvec[j][i].setPos(vec3(RandomFloat(0.0f, 5.0f), RandomFloat(5.0f, 9.0f), RandomFloat(0.0f, 8.0f)));
		}
	}

	//Set particles' initial position
	//allPart[0].setPos(vec3(-5.0f, 10.0f, 0.0f));
	//for (unsigned int i = 1; i < particlenum; ++i) {
	//	allPart[i].setPos(vec3(allPart[i-1].getPos())+vec3(2.0f,0.0f,0.0f));
	//}


	//Set particles initial position
	seflagPos(vecvec);             //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////FUNCTION HERE


	addHookeForces(vecvec);

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
	const int HookeForces = (particlenum*vertnum)-1;
	std::vector<Hooke*> Hookes;

	//Adding forces to the Hook vector


	//for (unsigned int j = 0; j < vertnum; j++) {
	//	for (int i = 0; i < HookeForces; i++)
	//	{
	//		if (i!=0) {

	//			Hookes.push_back(new Hooke(&(vecvec[j][i]), &(vecvec[j - 1][i]), ks, kd, rest));
	//			Hookes.push_back(new Hooke(&(vecvec[j][i]), &(vecvec[j][i - 1]), ks, kd, rest));
	//			Hookes.push_back(new Hooke(&(vecvec[j][i]), &(vecvec[j + 1][i]), ks, kd, rest));
	//			Hookes.push_back(new Hooke(&(vecvec[j][i]), &(vecvec[j][i + 1]), ks, kd, rest));
	//		}
	//		else {
	//			Hookes.push_back(new Hooke(&(vecvec[j][i]), &(vecvec[j - 1][i]), ks, kd, rest));
	//		}
	//		
	//	}
	//}

	//Adding the forces applied to the particle
	//for (unsigned int j = 0; j < vertnum; j++) {
	//	vecvec[j][1].addForce(Hookes[0]);
	//	for (int i = 1; i < particlenum-1; i++) {
	//		vecvec[j][i].addForce(&g);
	//		vecvec[j][i].addForce(Hookes[i]);
	//		vecvec[j][i+1].addForce(Hookes[i]);
	//	}
	//	//vecvec[j][8].addForce(Hookes[8]);
	//}
	
	for (unsigned int j = 0; j < vertnum; j++) {
		for (int i = 1; i < particlenum - 1; i++) {




		}
	}





	/*Leaving this for reference of how the spring forces are applied

	Hooke* h1 = new Hooke(&(allPart[0]), &(allPart[1]), ks, kd, rest);
	Hooke* h2 = new Hooke(&(allPart[1]), &(allPart[2]), ks, kd, rest);
	Hooke* h3 = new Hooke(&(allPart[2]), &(allPart[3]), ks, kd, rest);
	Hooke* h4 = new Hooke(&(allPart[3]), &(allPart[4]), ks, kd, rest);


	allPart[1].addForce(h1);
	allPart[1].addForce(h2);
	allPart[2].addForce(h2);
	allPart[2].addForce(h3);
	allPart[3].addForce(h3);
	allPart[3].addForce(h4);
	allPart[4].addForce(h4);
	}

	for (unsigned int i = 0; i < particlenum; i++) {
		allPart[1].addForce(&g);
		allPart[2].addForce(&g);
		allPart[3].addForce(&g);
		allPart[4].addForce(&g);
	}*/

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
			for (unsigned int j = 0; j < vertnum; j++) {
				for (unsigned int i = 1; i < particlenum; i++) {

					vecvec[j][i].setAcc(vecvec[j][i].applyForces(vecvec[j][i].getPos(), vecvec[j][i].getVel(), t, dt));

					for (unsigned int k = 0; k < 3; k++) {
						if (vecvec[j][i].getPos()[k] < cubecorner[k]) {
							vecvec[j][i].setVel(k, vecvec[j][i].getVel()[k] * -1.0f);
							vecvec[j][i].setPos(k, cubecorner[k]);
						}
						else if (vecvec[j][i].getPos()[k] > cubecorner[k] + d[k]) {
							vecvec[j][i].setVel(k, vecvec[j][i].getVel()[k] * -1.0f);
							vecvec[j][i].setPos(k, cubecorner[k] + d[k]);
						}

					}
				}
			}

			for (unsigned int j = 0; j < vertnum; j++) {
				for (unsigned int i = 1; i < particlenum; i++) {

					vecvec[j][i].setVel(vecvec[j][i].getVel() + dt * vecvec[j][i].getAcc());
					vecvec[j][i].translate(dt*vecvec[j][i].getVel());
				}
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
		for (unsigned int j = 0; j < vertnum; j++) {
			for (unsigned int i = 0; i < particlenum; i++) {
				app.draw(vecvec[j][i].getMesh());
			}
		}
		// draw demo objects
		//app.draw(cube);
		//app.draw(sphere);

		app.display();
	}

	app.terminate();

	return EXIT_SUCCESS;
}

