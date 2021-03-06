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
#include "Sphere.h"
using namespace glm;

// time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

Gravity g = Gravity(glm::vec3(0.0f, -9.8f, 0.0f));

float ks = 50.0f; //spring stiffness 
float rest = 0.1f;//spring rest length
float kd = 30.0f;//damping coefficient

///////TASK 2

const int spherenum = 7200;

const float sphereradius = 1.0f;


const int particlenum = 10;//particle num per vertex
const int vertnum = 10;
//Drag d = Drag();

int test = 0;

const int gridnum = 100;
const float planeScale = 500;

//Method tyhat creates random variables between 2 values
float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

//vec3 RandomPosition(float a, float b, Mesh plane) {
//	std::vector<int> thing;
//	bool test = true;
//	float xvalue;
//	float yvalue;
//	if (test) {
//		for (int i = 0; i < plane.getScale[0][0]-1; i++) {
//			thing.push_back(i);
//			if (i < plane.getScale[0][0] - 2) {
//				test = false;
//				xvalue = RandomFloat();
//			}
//		}
//	}
//	
//	float x = thing[rand()];
//	float z = thing[rand()];
//	vec3 finalpos = vec3(x,0.0f,z);
//	return finalpos;
//}


void clearGrid(std::vector<std::vector<std::vector<Sphere*>>> &optSpheres) {
	for (int i = 0; i < gridnum; i++) {
		for (int j = 0; j < gridnum; j++) {
			optSpheres[i][j].clear();
		}
	}
}

void updateGrid(std::vector<std::vector<std::vector<Sphere*>>> &optSpheres, Sphere *current, float numberofcolumns, int sidelength) {


	float xpos = current->getMesh().getPos()[0];
	float zpos = current->getMesh().getPos()[2];

	int column = floor(xpos / sidelength) + numberofcolumns / 2;
	int row = floor(zpos / sidelength) + numberofcolumns / 2;
	optSpheres[column][row].push_back(current);
	//std::cout << "size of cell " << column << "  " << row << " :" << optSpheres[column][row].size() << std::endl;

	//Right and front


	int column1 = floor((xpos + current->getRadius())/ sidelength) + numberofcolumns / 2;
	int row1 = floor((zpos + current->getRadius() )/ sidelength) + numberofcolumns / 2;
	if (column1 > 0 && row1 > 0 && column1 < numberofcolumns && row1 < numberofcolumns) {
		if (column1 != column || row1 != row) {
			optSpheres[column1][row1].push_back(current);
		}
	}
	

	//BVack and left

	int column2 = floor((xpos - current->getRadius()) / sidelength) + numberofcolumns / 2;
	int row2 = floor((zpos - current->getRadius()) / sidelength) + numberofcolumns / 2;
	if (column2 > 0 && row2 > 0 && column2 < numberofcolumns && row2 < numberofcolumns) {
		if ((column2 != column || row2 != row) && (column2 != column1 || row2 != row1)) {
			optSpheres[column2][row2].push_back(current);
		}
	}

	

	//Right and back
	int column3 = floor((xpos + current->getRadius()) / sidelength) + numberofcolumns / 2;
	int row3 = floor((zpos - current->getRadius()) / sidelength) + numberofcolumns / 2;
	if (column3 > 0 && row3 > 0 && column3 < numberofcolumns && row3 < numberofcolumns) {
		if ((column3 != column || row3 != row) && (column3 != column1 || row3 != row1 ) && (column3 != column2 || row3 != row2) ) {
			optSpheres[column3][row3].push_back(current);
		}
	}

	

	//Left and back
	//if (column2 != NULL && row2 != NULL) { optSpheres[column1][row1].clear(); }
	int column4 = floor((xpos - current->getRadius()) / sidelength) + numberofcolumns / 2;
	int row4 = floor((zpos + current->getRadius()) / sidelength) + numberofcolumns / 2;
	if (column4 > 0 && row4 > 0 && column4 < numberofcolumns && row4 < numberofcolumns) {
		if ((column4 != column || row4 != row) && (column4 != column1 || row4 != row1) && (column4 != column2 || row4 != row2) && (column4 != column3 || row4 != row3) ) {
			optSpheres[column4][row4].push_back(current);
		}
	}
	


}



std::vector<Vertex> giveColVertices(float  y, RigidBody &r) {

	std::vector<Vertex> colVert;


	//CHecking if any vertice is added to the plane, and if it is add it to the vector of vrtices
	for (Vertex v : r.getMesh().getVertices()) {
		//World position of given vertex:
		vec3 worlpos = mat3(r.getMesh().getModel()) * v.getCoord() + r.getPos();
		//Check if the world position of the vertices are coplliding with the plane
		if (worlpos[1] < y) {
			colVert.push_back(worlpos);
		}
	}

	return colVert;
}


// main function
int main()
{


	double framerateaverage;
	// create application
	Application app = Application::Application();
	app.initRender();
	Application::camera.setCameraPosition(glm::vec3(0.0f, 5.0f, 20.0f));

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	//std::cout << "Size of plane: " << plane.getVertices().size() << std::endl;
	// scale it up x5
	plane.scale(glm::vec3(planeScale, 0.0f, planeScale));
	//plane.scale(glm::vec3(1000.0f, 0.0f, 1000.0f));
	Shader lambert = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	plane.setShader(lambert);

	Shader blue = Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag");


	std::vector<std::vector<std::vector<Sphere*>>> optSpheres;

	for (int i = 0; i < gridnum; i++) {
		std::vector<std::vector<Sphere*>> columns;
		for (int j = 0; j < gridnum; j++) {
			std::vector<Sphere*> rows;
			columns.push_back(rows);
		}
		optSpheres.push_back(columns);
	}




	//Set up a cubic rigid body
	RigidBody rb = RigidBody();
	Mesh m = Mesh::Mesh(Mesh::CUBE);
	rb.setMesh(m);
	Shader rbShader = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	rb.getMesh().setShader(rbShader);
	rb.scale(vec3(1.0f, 3.0f, 1.0f));
	rb.setMass(1.0f);

	//Set up spheres
	std::vector<Sphere*> Spheres;

	//Sphere Spheres[spherenum];
	Mesh sphere = Mesh::Mesh("resources/models/sphere2.obj");
	Shader sphereshader = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	for (unsigned int i = 0; i < spherenum; i++) {
		Sphere *sp = new Sphere();

		sp->setMesh(sphere);
		sp->getMesh().setShader(sphereshader);

		//sp.scale(glm::vec3(sphereradius, sphereradius, sphereradius));
		sp->setRadius(1.0f);
		sp->setMass(1.0f);

		//updateGrid(optSpheres, sp, gridnum, (planeScale*2)/gridnum));
		Spheres.push_back(sp);
		//Spheres[i] = sp;



	}

	std::vector<std::pair<Sphere*,Sphere*>> collisionCheck;
	bool hasImpulseBeenAdded;

	//float r1;
	//float r2;
	for (unsigned int i = 0; i < spherenum; i++) {
		Spheres[i]->translate(glm::vec3(RandomFloat(-planeScale+1.0f, planeScale - 1.0f), Spheres[i]->getRadius(), RandomFloat(-planeScale + 1.0f, planeScale - 1.0f)));		//////////////////////////////////this i need to fix
		Spheres[i]->setVel(vec3(RandomFloat(-20.0f, 20.0f), 0.0f, RandomFloat(-20.0f, 20.0f)));
		//Spheres[i]->setAngVel(vec3(0.0f, 0.0f, 0.0f));
	}



	//Vertex




	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();


	//Timestep variables
	float t = 0.0f;
	const float dt = 0.01f;
	float currentTime = (GLfloat)glfwGetTime();
	float accumulator = 0.0f;

	//Cube Variables
	glm::vec3 cubecorner = glm::vec3(-29.0f, 0.0f, -29.0f);
	glm::vec3 d = glm::vec3(58.0f, 0.0f, 58.0f);

	std::vector<glm::vec3> dir;
	dir.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	dir.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
	dir.push_back(glm::vec3(0.0f, 0.0f, 1.0f));

	//Initial impulse
	bool impulseAppplied = false;
	vec3 impulseF = vec3(2.0f, 1.0f, 0.0f);
	//vec3 impulseF2 = vec3(3.0f, 0.0f, 0.0f);
	float impAppTime = 2.0f;



	//Ipulse pooint variables
	float e = 0.6f;

	bool stop = false;

	//// test
	//Sphere s1 = Sphere();
	//Sphere *s2 = &(Sphere());
	//s1.setMesh(sphere);
	//s1.getMesh().setShader(sphereshader);
	//s2->setMesh(sphere);
	//s2->getMesh().setShader(sphereshader);

	////sp.scale(glm::vec3(sphereradius, sphereradius, sphereradius));
	//s1.setRadius(1.0f);
	//s1.setMass(1.0f);
	//s2->setRadius(1.0f);
	//s2->setMass(1.0f);

	//std::vector<Sphere*> balls;


	////balls.push_back(&s1);			//MIRACLE
	//balls.push_back(s2);
	//balls[0]->translate(glm::vec3(0.0f, 10.0f, 0.0f));
	////s2->translate(glm::vec3(0.0f, 10.0f, 0.0f));
	//std::cout << optSpheres.size() << std::endl;



	//std::pair<Sphere*, Sphere*> BallPair;  //Check collision for pairs

	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{

		app.showFPS();
		float newTime = (GLfloat)glfwGetTime();
		float frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		//optSpheres.clear();
		while (accumulator >= dt) {
			/*
			**	SIMULATION
			*/

			//RIGID BODIES

			/*
		**	INTERACTION
		*/
		// Manage interaction

			app.doMovement(dt*10);

			clearGrid(optSpheres);

			for (unsigned int i = 0; i < spherenum; i++) {

				//Set Acceleration for Spheres
				//Spheres[i]->setAcc(Spheres[i]->applyForces(Spheres[i]->getPos(), Spheres[i]->getVel(), t, dt));
				Spheres[i]->setVel(Spheres[i]->getVel() + dt * Spheres[i]->getAcc());
				Spheres[i]->translate(dt*Spheres[i]->getVel());


				//Collisions with cushion UNOPTIMISED
				//for (unsigned int k = 0; k < 3; k++) {
				//	if (k != 1) {

				//		if (Spheres[i]->getPos()[k] < plane.getPos()[k] - plane.getScale()[k][k] + Spheres[i]->getRadius()) {
				//			//Spheres[i]->setVel(glm::vec3(0.0f));
				//			Spheres[i]->setVel(k, Spheres[i]->getVel()[k] * -1.0f);
				//			//Spheres[i]->translate(Spheres[i]->getPos() - plane.getScale()[k][k] + Spheres[i]->getRadius());
				//			Spheres[i]->setPos(k, plane.getPos()[k] - plane.getScale()[k][k] + Spheres[i]->getRadius());
				//		}

				//		else if (Spheres[i]->getPos()[k] > plane.getPos()[k] + plane.getScale()[k][k] - Spheres[i]->getRadius()) {
				//			//Spheres[i]->setVel(glm::vec3(0.0f));
				//			Spheres[i]->setVel(k, Spheres[i]->getVel()[k] * -1.0f);
				//			//Spheres[i]->translate(Spheres[i]->getPos() + plane.getScale()[k][k] - Spheres[i]->getRadius());
				//			Spheres[i]->setPos(k, plane.getPos()[k] + plane.getScale()[k][k] - Spheres[i]->getRadius());

				//		}

				//	}

				//}


				//for (unsigned int j = i + 1; j < spherenum; j++) {
				//	if (glm::distance(Spheres[i]->getPos(), Spheres[j]->getPos()) < Spheres[i]->getRadius() + Spheres[j]->getRadius()) {
				//		vec3 n = glm::normalize(Spheres[j]->getPos() - Spheres[i]->getPos());

				//		float displacement = Spheres[j]->getRadius() + Spheres[i]->getRadius() - glm::distance(Spheres[j]->getPos(), Spheres[i]->getPos());

				//		float ball1portion = glm::length(Spheres[i]->getVel()) / (glm::length(Spheres[i]->getVel()) + glm::length(Spheres[j]->getVel()));
				//		float ball2portion = glm::length(Spheres[j]->getVel()) / (glm::length(Spheres[i]->getVel()) + glm::length(Spheres[j]->getVel()));

				//		Spheres[i]->translate(displacement * -n * ball1portion);
				//		Spheres[j]->translate(displacement * n * ball2portion);

				//		vec3 vr = Spheres[j]->getVel() - Spheres[i]->getVel();
				//		//vec3 vr = Spheres[i]->getVel() - Spheres[j]->getVel();
				//		//vec3 n = glm::normalize(Spheres[i]->getVel());

				//		float jn = (-(1.0f + e)*dot(vr, n)) / (1.0f / Spheres[i]->getMass() + 1.0f / Spheres[j]->getMass());

				//		Spheres[i]->setVel(Spheres[i]->getVel() - (jn*n / Spheres[i]->getMass()));
				//		Spheres[j]->setVel(Spheres[j]->getVel() + (jn*n / Spheres[j]->getMass()));

				//	}
				//}

				updateGrid(optSpheres, Spheres[i], gridnum, (planeScale * 2) / gridnum);
			}


			
			for (int column = 0; column < gridnum; column++) {
				for (int row = 0; row < gridnum; row++) {
						for (int i = 0; i < optSpheres[column][row].size(); i++) {

							//Collision with Cushion
							if (column <= 0 || row <= 0 || column >= gridnum - 1 || row >= gridnum - 1) {
								for (unsigned int k = 0; k < 3; k++) {
									if (k != 1) {
										if (optSpheres[column][row][i]->getPos()[k] < plane.getPos()[k] - plane.getScale()[k][k] + optSpheres[column][row][i]->getRadius()) {
											//Spheres[i]->setVel(glm::vec3(0.0f));
											optSpheres[column][row][i]->setVel(k, optSpheres[column][row][i]->getVel()[k] * -1.0f);
											//Spheres[i]->translate(Spheres[i]->getPos() - plane.getScale()[k][k] + Spheres[i]->getRadius());
											optSpheres[column][row][i]->setPos(k, plane.getPos()[k] - plane.getScale()[k][k] + optSpheres[column][row][i]->getRadius());
										}

										else if (optSpheres[column][row][i]->getPos()[k] > plane.getPos()[k] + plane.getScale()[k][k] - optSpheres[column][row][i]->getRadius()) {
											//Spheres[i]->setVel(glm::vec3(0.0f));
											optSpheres[column][row][i]->setVel(k, optSpheres[column][row][i]->getVel()[k] * -1.0f);
											//Spheres[i]->translate(Spheres[i]->getPos() + plane.getScale()[k][k] - Spheres[i]->getRadius());
											optSpheres[column][row][i]->setPos(k, plane.getPos()[k] + plane.getScale()[k][k] - optSpheres[column][row][i]->getRadius());
										}
									}
								}
							}

							//Sphere Collision
							if (optSpheres[column][row].size() > 1) {
								for (unsigned int j = i + 1; j < optSpheres[column][row].size(); j++) {

									if (glm::distance(optSpheres[column][row][i]->getPos(), optSpheres[column][row][j]->getPos()) < optSpheres[column][row][i]->getRadius() + optSpheres[column][row][j]->getRadius()){

										if ((glm::length(optSpheres[column][row][i]->getVel()) > 0 || glm::length(optSpheres[column][row][j]->getVel()) > 0)) {



											/*hasImpulseBeenAdded = true;
											BallPair = std::pair<Sphere*, Sphere*>(optSpheres[column][row][i], optSpheres[column][row][j]);
											for (int t = 0; t < collisionCheck.size(); t++) {
												if (collisionCheck[t] == BallPair) {
													hasImpulseBeenAdded = false;
												}
											}*/

											//if (hasImpulseBeenAdded) {
												vec3 n = glm::normalize(optSpheres[column][row][j]->getPos() - optSpheres[column][row][i]->getPos());

												float displacement = optSpheres[column][row][j]->getRadius() + optSpheres[column][row][i]->getRadius() - glm::distance(optSpheres[column][row][j]->getPos(), optSpheres[column][row][i]->getPos());

												float ball1portion = glm::length(optSpheres[column][row][i]->getVel()) / (glm::length(optSpheres[column][row][i]->getVel()) + glm::length(optSpheres[column][row][j]->getVel()));
												float ball2portion = glm::length(optSpheres[column][row][j]->getVel()) / (glm::length(optSpheres[column][row][i]->getVel()) + glm::length(optSpheres[column][row][j]->getVel()));

												optSpheres[column][row][i]->translate(displacement * -n * ball1portion);
												optSpheres[column][row][j]->translate(displacement * n * ball2portion);

												vec3 vr = optSpheres[column][row][j]->getVel() - optSpheres[column][row][i]->getVel();
												//vec3 vr = optSpheres[column][row][i]->getVel() - optSpheres[column][row][j]->getVel();
												//vec3 n = glm::normalize(optSpheres[column][row][i]->getVel());


												float jn = (-(1.0f + e)*dot(vr, n)) / (1.0f / optSpheres[column][row][i]->getMass() + 1.0f / optSpheres[column][row][j]->getMass());
												optSpheres[column][row][i]->setVel(optSpheres[column][row][i]->getVel() - (jn*n / optSpheres[column][row][i]->getMass()));
												optSpheres[column][row][j]->setVel(optSpheres[column][row][j]->getVel() + (jn*n / optSpheres[column][row][j]->getMass()));






												//collisionCheck.push_back(BallPair);



											//}
									}

									}

								}
							}

					}

				}
			}

			//collisionCheck.clear();
			

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

		// test
		//app.draw(balls[0]->getMesh());
		//app.draw(s2->getMesh());



		// draw particles
		//app.draw(particle1.getMesh());
		//for (unsigned int j = 0; j < vertnum; j++) {
		//	for (unsigned int i = 0; i < particlenum; i++) {
		//		app.draw(vecvec[j][i].getMesh());
		//	}
		//}

		for (unsigned int i = 0; i < spherenum; i++) {
			app.draw(Spheres[i]->getMesh());
		}

		// draw demo objects
		//app.draw(cube);
		//app.draw(sphere);
		app.display();
	}

	app.terminate();

	return EXIT_SUCCESS;
}

