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

//Task 4 Variables


//float coneheight = 20.0f; //Height of cone
//vec3 coneorigin = vec3(0.0f,-21.0f,0.0f); //Tip of the cone
//vec3 coneaxis = vec3(0.0f, 1.0f, 0.0f); //Central axis in the middle of the copne so that we can find thye projections from it
//float coneradius = 5.0f; //radius at the top of the cylinder 

float coneheight = 20.0f; //Height of cone
vec3 coneorigin = vec3(0.0f, -15.0f, 0.0f); //Tip of the cone(origin)
vec3 coneaxis = vec3(0.0f, 1.0f, 0.0f); //Central axis in the middle of the copne so that we can find thye projections from it
float coneradius = 4.0f; //radius at the top of the cylinder 


glm::vec3 CalculateWindForce(glm::vec3 pos) {

	//Variables:
	//coneorigin
	//coneaxis
	//coneheight
	//coneradius

	float distanceforYaxis = dot(pos - coneorigin, coneaxis); //Calculates distance from particle to height of the cone. We are doing this for the "if" in the y axis(see below)

	//Set the force to 0 if it's within the height of the axis
	std::cout << "distanceforYaxis" << distanceforYaxis << std::endl;
	//std::cout << "coneheight" << coneheight << std::endl;

	if (distanceforYaxis < 0 || distanceforYaxis > coneheight) {
		return vec3(0.0f);
	}


	//Now we will calculate the radius of the cone, at the height of the current possition of the particle by applying a math formula
	float currentr = (distanceforYaxis / coneheight)*coneradius;
	std::cout << "currentr" << currentr << std::endl;
	//Calculate distance from current positionb to the cemntral axis for the cone, in orser touse it for the if statement for the x axis
	float distanceforXaxis = length((pos - coneorigin) - distanceforYaxis * coneaxis);
	//float distanceforXaxis = ;




	//Set the force to 0 if the distance from the central axis is greater than the radius of the cone at the current y position of the particle

	if (distanceforXaxis > currentr) {
		return vec3(0.0f);
	}



	//std::cout << "distanceforXaxis" << distanceforXaxis << std::endl;

	//Now we will calculate the wind force the particle receives within the cube, based on its current position  (WIND FORCE CALCULATION)

	//Calculate the force so that it generates is maximal at the bottom of the cone and null at the top.
	float ymagnitude = 1.0f - (distanceforYaxis / coneheight);
	//float ymagnitude = 1.0f - (distanceforYaxis / coneheight);

	//Calculate the force so that it's maximal along the central vertical axis of the cone at any given height and decreases radially until it is null at the edge of the cone.
	float xmagnitude = 1.0f - (distanceforXaxis / coneradius);

	//Calculate total decrease of the force due to distance from central axis and height, by multiplying them together
	float magnitude = ymagnitude * xmagnitude;

	//Now we will calculate the direction that the force will be applied to the particle, by geting the vector from the particle position tpo the cube origin
	vec3 forceangle = pos - coneorigin;

	//We calculate the force of the wnd, and normalise it depending on the forceangle(above)
	vec3 Fcurrent = normalize(forceangle)*magnitude;

	//Return the force

	std::cout << "Fcurrent" << to_string(Fcurrent) << std::endl;

	return Fcurrent;

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


	// create particle
	Mesh particle1 = Mesh::Mesh(Mesh::QUAD);
	//scale it down (x.1), translate it up by 2.5 and rotate it by 90 degrees around the x axis
	particle1.translate(glm::vec3(0.0f, 2.5f, 0.0f));
	particle1.scale(glm::vec3(.1f, .1f, .1f));
	particle1.rotate((GLfloat)M_PI_2, glm::vec3(1.0f, 0.0f, 0.0f));
	particle1.setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));


	//Create particle via class:

	const int particlenum = 50;
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
	//	Mesh cylinder = Mesh::Mesh("resources/models/cone2.obj");
	//	cylinder.translate(glm::vec3(1.0f, .5f, 0.0f));
	//	cylinder.setShader(lambert);

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
	float absoluteu;
	glm::vec3 g = glm::vec3(0.0f, -9.8f, 0.0f);


	//glm::vec3 ParticleVelocity[particlenum]; //initial; velocity for multiple particles

	////Task 1&2
	for (unsigned int i = 0; i < particlenum; ++i) {
		allPart[i].setVel(vec3(RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 8.0f)));
		allPart[i].setPos(vec3(RandomFloat(0.0f, 5.0f), RandomFloat(0.0f, 8.0f), RandomFloat(0.0f, 8.0f)));
	}

	//Task 2
	//for (unsigned int i = 0; i < particlenum; ++i) {
		//allPart[i].setVel(vec3(0.0f));
		//allPart[i].setPos(vec3(0.0f, 2.5f, 0.0f));
		//std::cout << "mass: " << allPart[i].getMass() << std::endl;
	//}

	/*allPart[0].setPos(vec3(0.0f, 2.0f, 0.0f));
	allPart[1].setPos(vec3(2.0f, 2.0f, 0.0f));
	allPart[2].setPos(vec3(-2.0f, 2.0f, 0.0f));*/

	float t = 0.0f;
	const float dt = 0.01f;
	float currentTime = (GLfloat)glfwGetTime();
	float accumulator = 0.0f;


	////task 4 lefover variables


	float l = 5.0f;


	//Fwind is P*S, where P is pressure and S is the area of the cone at a particular point, also we know the P1S1=P2S2
	float Pmax = 50.0f; //Pressure at the bottom of the cone
	float conearea;
	//float currentr; //current r for paryticles
	float Pcurrent;


	float cheight; //current height of particle


	vec3 Fwind; //The force that will be applied to the particle if it's within the cone



	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{
		float newTime = (GLfloat)glfwGetTime();
		float frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;



		////// Set frame time
		//GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		////// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		//currentFrame *= 1.5f;
		//deltaTime = currentFrame - lastFrame;
		//lastFrame = currentFrame;



		/*
				**	INTERACTION
				*/
				// Manage interaction
		app.doMovement(dt);


		while (accumulator >= dt)
		{
			/*
			**	SIMULATION
			*/
			for (unsigned int i = 0; i < particlenum; i++) {

				// compute forces
				// gravity
				Fg = allPart[i].getMass() * g;
				Fa = vec3(0.0f);

				Fwind = CalculateWindForce(allPart[i].getPos()) * 200;
				//std::cout << "Fwind" << to_string(Fwind) << std::endl;


				//Fwind = vec3(0.0f); ///////////////////SET FWIND TO 0 FOR TESTINGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

				//Calculate Fwind.We know that Fwind = P*S, where S is pressure at a point and S is the area.

				//currentr = glm::length(vec3(allPart[i].getPos().x ,0.0f, allPart[i].getPos().z));
				//conearea = M_PI * currentr *currentr;
				////We kbnow that P1A1=P2A2, so 
				//Pcurrent = Fwindymax / conearea;
				////SO, Fwind = P*Area, so				////////////////////////MY VERSION OF FWIND
				////std::cout << r << std::endl;
				//cheight = h - allPart[i].getPos().y;
				//Fwind = vec3(0.0f, currentr*Pcurrent, 0.0f);


				//Wind force
				//if (allPart[i].getAcc().y==0){

				//}


				// aerodynamic drag
				if (glm::length(allPart[i].getVel()) == 0.0f) {	//If velocity is 0, set the aerodynamic drag to 0
					Fa = vec3(0.0f);
				}
				else {
					//std::cout << glm::to_string(allPart[i].getVel()) << std::endl;
					absoluteu = length(allPart[i].getVel()); /////////////////////////THIS LINE IS BREAKING MY CODE
					//std::cout << "Abs: " << absoluteu << std::endl;


					e = -allPart[i].getVel() / absoluteu;
					Fa = 0.5*density*absoluteu*absoluteu*coefficient*area*e;
					//Fa = vec3(0.0f);
				}




				Ftotal = Fg + Fa + Fwind;

				//std::cout << "Ftotal" << to_string(Ftotal) << std::endl;
				//std::cout << "ftotal: " << to_string(Ftotal) << std::endl;


				//TESTING FOR TASK 4
				// acceleration
				allPart[i].setAcc(Ftotal / allPart[i].getMass());
				//std::cout << "ftotal: " << to_string(Fa) << std::endl;


				for (unsigned int j = 0; j < 3; j++) {
					if (allPart[i].getPos()[j] < cubecorner[j]) {
						allPart[i].setVel(j, allPart[i].getVel()[j] * -1.0f);
						allPart[i].setPos(j, cubecorner[j]);

					}
					else if (allPart[i].getPos()[j] > cubecorner[j] + d[j]) {
						allPart[i].setVel(j, allPart[i].getVel()[j] * -1.0f);
						allPart[i].setPos(j, cubecorner[j] + d[j]);
					}

				}

				// integrate

				allPart[i].translate(dt*allPart[i].getVel());
				allPart[i].setVel(allPart[i].getVel() + dt * allPart[i].getAcc());

				

				/*                                            TASK 2

				if (i == 1) {  //Forward Euler
					// acceleration
					allPart[i].setAcc(Ftotal / allPart[i].getMass());
					//std::cout << "ftotal: " << to_string(Fa) << std::endl;


					for (unsigned int j = 0; j < 3; j++) {
						if (allPart[i].getPos()[j] < cubecorner[j]) {
							allPart[i].setVel(j, allPart[i].getVel()[j] * -1.0f);
							allPart[i].setPos(j, cubecorner[j]);

						}
						else if (allPart[i].getPos()[j] > cubecorner[j] + d[j]) {
							allPart[i].setVel(j, allPart[i].getVel()[j] * -1.0f);
							allPart[i].setPos(j, cubecorner[j] + d[j]);
						}

					}

					// integrate


					allPart[i].setVel(allPart[i].getVel() + dt * allPart[i].getAcc());

					allPart[i].translate(dt*allPart[i].getVel());
				}

				if (i == 2) {
					allPart[i].setAcc(Ftotal / allPart[i].getMass());
					//std::cout << "ftotal: " << to_string(Fa) << std::endl;


					for (unsigned int j = 0; j < 3; j++) {
						if (allPart[i].getPos()[j] < cubecorner[j]) {
							allPart[i].setVel(j, allPart[i].getVel()[j] * -1.0f);
							allPart[i].setPos(j, cubecorner[j]);

						}
						else if (allPart[i].getPos()[j] > cubecorner[j] + d[j]) {
							allPart[i].setVel(j, allPart[i].getVel()[j] * -1.0f);
							allPart[i].setPos(j, cubecorner[j] + d[j]);
						}

					}

					// integrate
					allPart[i].translate(dt*allPart[i].getVel());

					allPart[i].setVel(allPart[i].getVel() + dt * allPart[i].getAcc());


				}
				*/

				//if (allPart[i].getVel() == vec3(0.0f)) {
				//	allPart[i].translate(vec3(0.0f,2.5f,0.0f));
				//}
				//else {
				//	allPart[i].translate(dt*allPart[i].getVel());
				//}
				//std::cout << glm::to_string(allPart[i].getVel());


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
		//app.draw(particle1);	////////////////////////////////////////////////////put partile here
		for (unsigned int i = 0; i < particlenum; ++i) {
			app.draw(allPart[i].getMesh());
		}
		// draw demo objects
		//app.draw(cube);
		//app.draw(sphere);
		//app.draw(cylinder);

		app.display();
	}

	app.terminate();

	return EXIT_SUCCESS;
}
