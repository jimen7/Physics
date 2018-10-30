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
#include "RigidBody.h"
using namespace glm;

// time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

Gravity g = Gravity(glm::vec3(0.0f, -9.8f, 0.0f));

float ks = 50.0f; //spring stiffness 
float rest = 0.1f;//spring rest length
float kd = 30.0f;//damping coefficient


const int particlenum = 10;//particle num per vertex
const int vertnum = 10;
//Drag d = Drag();

int test = 0;

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
	for (unsigned int x = 0; x < 10; ++x) {
		for (unsigned int y = 0; y < 10; ++y) {
			pvert[x][y].setPos(inpos + vec3(x*gap, 0.0f, -1.0f*y*gap));
		}
	}
}

void seflagPoscloth(std::vector<std::vector<Particle>> &pvert) {
	vec3 inpos = vec3(-5.0f, 10.0f, 0.0f);
	float gap = 1.0f;
	for (unsigned int x = 0; x < 10; ++x) {
		for (unsigned int y = 0; y < 10; ++y) {
			pvert[x][y].setPos(inpos + vec3(x*gap, -1.0f*y*gap, 0.0f));
			//std::cout << pvert[x][y].getPos()[0] <<";"<< pvert[x][y].getPos()[1] << ";" << pvert[x][y].getPos()[2] << std::endl;
		}
	}
}


void addHookeForces(std::vector<std::vector<Particle>> &pvert) {
	for (unsigned int j = 0; j < vertnum; j++) {
		for (int i = 0; i < particlenum - 1; i++)
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

				vec3 topnormal = (cross((pvert[i - 1][j - 1].getPos() - pvert[i][j].getPos()), (pvert[i - 1][j].getPos() - pvert[i][j].getPos()))) / glm::length(cross((pvert[i - 1][j - 1].getPos() - pvert[i][j].getPos()), (pvert[i - 1][j].getPos() - pvert[i][j].getPos())));
				float toparea = 0.5f * glm::length(cross((pvert[i - 1][j - 1].getPos() - pvert[i][j].getPos()), (pvert[i - 1][j].getPos() - pvert[i][j].getPos())));
			}



		}
	}
}

void addClothForces(std::vector<std::vector<Particle>> &pvert) {

	for (unsigned int j = 0; j < vertnum; j++) {

		for (int i = 0; i < particlenum; i++)
		{

			if (j != 0) {

				if (i == 0) {   ////LEFT SIDE
					if (j == 9) { //Bottom left particle
						Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
						Hooke* right = new Hooke(&(pvert[j][i]), &(pvert[j][i + 1]), ks, kd, rest);

						pvert[j][i].addForce(&g);

						pvert[j][i].addForce(top);
						pvert[j][i].addForce(right);
					}
					else {
						Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
						Hooke* right = new Hooke(&(pvert[j][i]), &(pvert[j][i + 1]), ks, kd, rest);
						Hooke* bottom = new Hooke(&(pvert[j][i]), &(pvert[j + 1][i]), ks, kd, rest);

						Wind * wind1 = new Wind(&(pvert[j][i]), &(pvert[j][i + 1]), &(pvert[j + 1][i])); //particle, particle right, particle down

						pvert[j][i].addForce(&g);

						pvert[j][i].addForce(top);
						pvert[j][i].addForce(bottom);
						pvert[j][i].addForce(right);


						pvert[j][i].addForce(wind1);
						pvert[j][i + 1].addForce(wind1);
						pvert[j + 1][i].addForce(wind1);
					}
				}



				if (j == 9) { //////BOTTOM SIDE
					if (i == 0) { //Bottom left particle
						/*Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
						Hooke* right = new Hooke(&(pvert[j][i]), &(pvert[j][i + 1]), ks, kd, rest);

						pvert[j][i].addForce(&g);
						pvert[j][i].addForce(top);
						pvert[j][i].addForce(right);*/
					}
					else if (i == 9) { //bottom right particle
						Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
						Hooke* left = new Hooke(&(pvert[j][i]), &(pvert[j][i - 1]), ks, kd, rest);

						pvert[j][i].addForce(&g);
						pvert[j][i].addForce(top);
						pvert[j][i].addForce(left);
					}
					else {
						Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
						Hooke* left = new Hooke(&(pvert[j][i]), &(pvert[j][i - 1]), ks, kd, rest);
						Hooke* right = new Hooke(&(pvert[j][i]), &(pvert[j][i + 1]), ks, kd, rest);

						pvert[j][i].addForce(&g);

						pvert[j][i].addForce(top);
						pvert[j][i].addForce(left);
						pvert[j][i].addForce(right);
					}
				}

				if (i == 9) { //RIGHT SIDE
					if (j == 9) { //Bottom right particle
						/*Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
						Hooke* left = new Hooke(&(pvert[j][i]), &(pvert[j][i - 1]), ks, kd, rest);*/
					}
					else {
						Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
						Hooke* left = new Hooke(&(pvert[j][i]), &(pvert[j][i - 1]), ks, kd, rest);
						Hooke* bottom = new Hooke(&(pvert[j][i]), &(pvert[j + 1][i]), ks, kd, rest);

						pvert[j][i].addForce(&g);

						pvert[j][i].addForce(top);
						pvert[j][i].addForce(left);
						pvert[j][i].addForce(bottom);
					}
				}


				if (i != 0 && j != 0 && i != 9 && j != 9) {  //INSIDE PARTICLES
					Hooke* top = new Hooke(&(pvert[j][i]), &(pvert[j - 1][i]), ks, kd, rest);
					Hooke* left = new Hooke(&(pvert[j][i]), &(pvert[j][i - 1]), ks, kd, rest);
					Hooke* bottom = new Hooke(&(pvert[j][i]), &(pvert[j + 1][i]), ks, kd, rest);
					Hooke* right = new Hooke(&(pvert[j][i]), &(pvert[j][i + 1]), ks, kd, rest);

					pvert[j][i].addForce(&g);

					pvert[j][i].addForce(top);
					pvert[j][i].addForce(left);
					pvert[j][i].addForce(bottom);
					pvert[j][i].addForce(right);

					Wind * wind1 = new Wind(&(pvert[j][i]), &(pvert[j][i + 1]), &(pvert[j + 1][i])); //particle, particle right, particle down

					Wind * wind2 = new Wind(&(pvert[j][i + 1]), &(pvert[j + 1][i + 1]), &(pvert[j + 1][i])); //particle right, particle right and down, particle down

					pvert[j][i].addForce(wind1);
					pvert[j][i + 1].addForce(wind1);
					pvert[j + 1][i].addForce(wind1);

					pvert[j][i + 1].addForce(wind2);
					pvert[j + 1][i + 1].addForce(wind2);
					pvert[j + 1][i].addForce(wind2);

				}

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
	std::cout << "Size of plane: " << plane.getVertices().size() << std::endl;
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


	/*
	////////////////////////////////////////Create vertice of vertices
	std::vector<std::vector<Particle>> vecvec(vertnum);
	for (unsigned int j = 0; j < particlenum; ++j) {
		for (unsigned int i = 0; i < vertnum; ++i) {
			vecvec[j].push_back(Particle::Particle());
			vecvec[j][i].getMesh().setShader(blue);
		}
	}

	*/




	//Set particles initial velocity
	//for (unsigned int j = 0; j < particlenum; ++j) {
	//	for (unsigned int i = 0; i < vertnum; ++i) {
	//		//vecvec[j][i].setVel(vec3(RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 10.0f), RandomFloat(-10.0f, 8.0f)));
	//		vecvec[j][i].setVel(vec3(0.0f));
	//		//vecvec[j][i].setPos(vec3(RandomFloat(0.0f, 5.0f), RandomFloat(5.0f, 9.0f), RandomFloat(0.0f, 8.0f)));
	//	}
	//}


	//Set particles initial position
	//seflagPos(vecvec); 
	//seflagPoscloth(vecvec);

	//addHookeForces(vecvec);
	//addClothForces(vecvec);

	// create demo objects (a cube and a sphere)
	Mesh sphere = Mesh::Mesh("resources/models/sphere.obj");
	sphere.translate(glm::vec3(-1.0f, 1.0f, 0.0f));
	sphere.setShader(lambert);
	Mesh cube = Mesh::Mesh("resources/models/cube.obj");
	cube.translate(glm::vec3(1.0f, .5f, 0.0f));
	cube.setShader(lambert);

	//Set up a cubic rigid body
	RigidBody rb = RigidBody();
	Mesh m = Mesh::Mesh(Mesh::CUBE);
	rb.setMesh(m);
	Shader rbShader = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	rb.getMesh().setShader(rbShader);


	// rigid body motion values
	rb.translate(vec3(0.0f, 5.0f, 0.0f));
	rb.setVel(vec3(2.0f, 2.0f, 0.0f));
	rb.setAngVel(vec3(0.0f, 2.0f, 5.0f));

	//add forces to Rigid body
	rb.addForce(&g);




	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();


	//Timestep variables
	float t = 0.0f;
	const float dt = 0.001f;
	float currentTime = (GLfloat)glfwGetTime();
	float accumulator = 0.0f;

	//Cube Variables
	glm::vec3 cubecorner = glm::vec3(-5.0f, 0.0f, -5.0f);
	glm::vec3 d = glm::vec3(10.0f);

	std::vector<glm::vec3> dir;
	dir.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	dir.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
	dir.push_back(glm::vec3(0.0f, 0.0f, 1.0f));

	//Hooke force
	const int HookeForces = (particlenum*vertnum) - 1;
	std::vector<Hooke*> Hookes;

	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{
		float newTime = (GLfloat)glfwGetTime();
		float frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		
		while (accumulator >= dt) {
			/*
			**	SIMULATION
			*/

			//RIGID BODIES
			
			/*
		**	INTERACTION
		*/
		// Manage interaction
			app.doMovement(dt);



			/*if (test == 0) {*/
				rb.setAcc(rb.applyForces(rb.getPos(), rb.getVel(), t, dt));
		/*	}
			else if (test == 1) {
				rb.setAcc(rb.applyForces(rb.getPos() - 1.8f, rb.getVel(), t, dt));
			}
			else if (test == 2) {
				rb.setAcc(rb.applyForces(rb.getPos() + 1.8f, rb.getVel(), t, dt));
			}*/
			
			
			//rb.setAngAccl(vec3(0.0f,1.0f,0.0f));
			//integration (rotation)
			rb.setAngVel(rb.getAngVel() + dt * rb.getAngAcc());
			//crete skew symmetric matrix for w
			glm::mat3 angVelSkew = glm::matrixCross3(rb.getAngVel());
			//create 3x3 rotation matrix from rb rotation matrix
			glm::mat3 R = glm::mat3(rb.getRotate());
			//update rotation matrix
			R += dt * angVelSkew * R;
			R = glm::orthonormalize(R);
			rb.setRotate(glm::mat4(R));



			//WORKS WITH CENTER OF MASS
			//for (unsigned int k = 0; k < 3; k++) {
			//	if (rb.getPos()[k] - 1.0f < cubecorner[k]) {
			//		rb.setVel(k, rb.getVel()[k] * -1.0f);
			//		rb.setPos(k - 1.0f, cubecorner[k]);
			//	}
			//	else if (rb.getPos()[k] + 1.0f > cubecorner[k] + d[k]) {
			//		rb.setVel(k, rb.getVel()[k] * -1.0f);
			//		rb.setPos(k + 1.0f, cubecorner[k] + d[k]);
			//	}
			//}
			test = 0;
			//Will try to make it work withy vertices
			for (Vertex vert : rb.getMesh().getVertices())
			{
				for (unsigned int k = 0; k < 3; k++) {
					vec4 worldcoord = rb.getMesh().getModel()*vec4(vert.getCoord(), 1.0f);
					vec3 pointofcontact;
					float sign;
					if (worldcoord[k] < cubecorner[k] ) {
						test = 1;
						pointofcontact = vec3(worldcoord);
						rb.setVel(k, rb.getVel()[k] * -0.5f);
						sign = rb.getVel()[k] / length(rb.getVel()[k]);
						//rb.setPos(k, cubecorner[k] + 1.8f);
						//rb.translate((0.1f * dir[k] * sign));

						rb.translate(((worldcoord[k]-cubecorner[k]) * dir[k] )* -1.0f);

						rb.setAngVel(rb.getAngVel()*(-0.5f));
						//std::cout << "Collision1: ("  << pointofcontact.x << "," << pointofcontact.y << "," << pointofcontact.z << ")" << std::endl;
						//std::cout << "worldpos: (" << worldcoord[k] << ")" << std::endl;
					}
					else if (worldcoord[k] > cubecorner[k] + d[k] ) {
						test = 2;
						pointofcontact = vec3(worldcoord);
						rb.setVel(k, rb.getVel()[k] * -0.5f);
						sign = rb.getVel()[k] / length(rb.getVel()[k]);
						//rb.setPos(k, cubecorner[k] + d[k] -1.8f);
						//rb.translate(((cubecorner[k] + d[k]) * dir[k] )* sign);

						rb.translate(((worldcoord[k] - (cubecorner[k]+d[k]) ) * dir[k])* -1.0f);
						//rb.translate((0.1f * dir[k] * sign));
						rb.setAngVel(rb.getAngVel()*(-0.5f));
						//std::cout << "Collision2: (" << pointofcontact.x << "," << pointofcontact.y << "," << pointofcontact.z << ")" << std::endl;
					}
				}
			}
			
			rb.setVel(rb.getVel() + dt * rb.getAcc());
			rb.translate(dt*rb.getVel());



			for (unsigned int j = 0; j < vertnum; j++) {
				for (unsigned int i = 0; i < particlenum; i++) {



					//vecvec[j][i].setAcc(vecvec[j][i].applyForces(vecvec[j][i].getPos(), vecvec[j][i].getVel(), t, dt));

				/*	for (unsigned int k = 0; k < 3; k++) {
						if (vecvec[j][i].getPos()[k] < cubecorner[k]) {
							vecvec[j][i].setVel(k, vecvec[j][i].getVel()[k] * -1.0f);
							vecvec[j][i].setPos(k, cubecorner[k]);
						}
						else if (vecvec[j][i].getPos()[k] > cubecorner[k] + d[k]) {
							vecvec[j][i].setVel(k, vecvec[j][i].getVel()[k] * -1.0f);
							vecvec[j][i].setPos(k, cubecorner[k] + d[k]);
						}

					}*/




				}
			}

			for (unsigned int j = 0; j < vertnum; j++) {
				for (unsigned int i = 0; i < particlenum; i++) {



					//vecvec[j][i].setVel(vecvec[j][i].getVel() + dt * vecvec[j][i].getAcc());
					//vecvec[j][i].translate(dt*vecvec[j][i].getVel());
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
		//for (unsigned int j = 0; j < vertnum; j++) {
		//	for (unsigned int i = 0; i < particlenum; i++) {
		//		app.draw(vecvec[j][i].getMesh());
		//	}
		//}
		// draw demo objects
		//app.draw(cube);
		//app.draw(sphere);
		app.draw(rb.getMesh());
		app.display();
	}

	app.terminate();

	return EXIT_SUCCESS;
}

