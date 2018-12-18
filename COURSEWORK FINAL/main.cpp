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






const float sphereradius = 1.0f;


const int particlenum = 10;//particle num per vertex
const int vertnum = 10;
//Drag d = Drag();

int test = 0;

const int gridnum = 100;
const float planeScale = 30;


//DEMO VALUES
bool triangle = true;
//const int spherenum = 6; //For DEMO
//const int spherenum = 15; //For DEMO(opssite spin)
//const int spherenum = 16; //For TRIANGLE, we can have 16, 22, 29, 37,46, 56,67
const int spherenum = 16;
//float mue = 0.7f;		//For Ricochet Testing
float mue = 0.12f;	//For Straight Testing


//Method tyhat creates random variables between 2 values
float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

std::vector<Vertex> giveBallColVertices(float  y, Sphere *s) {

	std::vector<Vertex> colVert;


	//CHecking if any vertice is added to the plane, and if it is add it to the vector of vrtices
	for (Vertex v : s->getMesh().getVertices()) {
		//World position of given vertex:
		vec3 worlpos = mat3(s->getMesh().getModel()) * v.getCoord() + s->getPos();
		//Check if the world position of the vertices are coplliding with the plane
		if (worlpos[1] < y) {
			colVert.push_back(worlpos);
		}
	}

	return colVert;
}

void setTriangleBallPos(std::vector<Sphere*> &balls, vec3 &givenPos) {

	int oddxiteration;
	int oddziteration;
	int evenziteration;
	int evenxiteration;

	int ballline = 2;
	for (int i = 0; i < spherenum; i++) {
		if (i == 0) {
			balls[i]->setPos(givenPos);
		}
		else if (i == 1) {
			balls[i]->translate(balls[0]->getPos() + vec3(0, 0, 25));
		}
		else {
			if (ballline % 2 == 0) {

				evenxiteration = ballline -1;

				if (ballline == 2) {
					evenziteration = 0;
				}
				else {
					evenziteration = ballline;
				}



				for (int j=0, h = i; j < ballline; j++,h++) {
					balls[h]->translate(balls[0]->getPos() + vec3(evenxiteration*balls[0]->getRadius(),0, (-2-evenziteration)*balls[0]->getRadius()));
					evenxiteration = evenxiteration-2;
					i = h;
				}
				ballline++;
			}
			else {

				oddxiteration = ballline - 1;

				if (ballline == 3) {
					oddziteration = 0;
				}
				else {
					oddziteration = ballline-1;
				}


				for (int j = 0, h = i; j < ballline; j++, h++) {

					balls[h]->translate(balls[0]->getPos() + vec3(oddxiteration*balls[0]->getRadius(), 0, (-4 - oddziteration)*balls[0]->getRadius()));
					oddxiteration = oddxiteration - 2;
					i = h;
				}


				ballline++;
			}

		}
	}
}

// main function
int main()
{


	double framerateaverage;
	// create application
	Application app = Application::Application();
	app.initRender();
	Application::camera.setCameraPosition(glm::vec3(0.0f, 5.0f, 50.0f));

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up *planescale variable
	plane.scale(glm::vec3(planeScale, 0.0f, planeScale));
	Shader lambert = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	plane.setShader(lambert);

	//Set up spheres
	std::vector<Sphere*> Spheres;

	//Sphere Spheres[spherenum];
	Mesh sphere = Mesh::Mesh("resources/models/sphere.obj");
	Shader sphereshader = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	Shader whiteball = Shader("resources/shaders/physics.vert", "resources/shaders/solid_grey.frag");
	for (unsigned int i = 0; i < spherenum; i++) {
		Sphere *sp = new Sphere();

		sp->setMesh(sphere);
		if (triangle) {
			if (i == 1) {
				sp->getMesh().setShader(whiteball);		//SET WHITE BALL SHADER
			}
			else{
				sp->getMesh().setShader(sphereshader);		//SET OTHER SHADERS
			}
		}
		else{
			sp->getMesh().setShader(sphereshader);
		}
		
		sp->setRadius(1.0f);
		sp->setMass(1.0f);
		Spheres.push_back(sp);
	}




	//SET POSITIONS DEPENDING ON TRIANGLE BOOLEAN
	if (triangle) {
		setTriangleBallPos(Spheres, vec3(0.0f, 1.0f, -14.0f));
	}
	else {
		//SET RANDOM POSITIONS
		for (unsigned int i = 0; i < spherenum; i++) {
			Spheres[i]->translate(glm::vec3(RandomFloat(-planeScale + 1.0f, planeScale - 1.0f), Spheres[i]->getRadius(), RandomFloat(-planeScale + 1.0f, planeScale - 1.0f)));
		}
	}

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

	//Initial impulse Values
	bool impulseAppplied = false;
	vec3 impulseF;
	float impAppTime = 4.0f;		//Timer before impulse is applied
	int impulsenumcheck = 0;


	//Ipulse point variables
	float e = 0.6f;


	for (int i = 0; i < spherenum; i++) {
		//Spheres[i]->addForce(&g);		//Gravity doesn't work, but it's not needed for this part
	}


	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{

		app.showFPS();
		float newTime = (GLfloat)glfwGetTime();
		float frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;
		bool AccVel = false;


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

			for (unsigned int i = 0; i < spherenum; i++) {

				//Set Acceleration for Spheres
				Spheres[i]->setAcc(Spheres[i]->applyForces(Spheres[i]->getPos(), Spheres[i]->getVel(), t, dt));
				Spheres[i]->setVel(Spheres[i]->getVel() + dt * Spheres[i]->getAcc());
				Spheres[i]->translate(dt*Spheres[i]->getVel());


				// rotate body
			//rb.setAngAccl(vec3(0.0f,1.0f,0.0f));
			//integration (rotation)
				Spheres[i]->setAngVel(Spheres[i]->getAngVel() + dt * Spheres[i]->getAngAcc());
				//crete skew symmetric matrix for w
				glm::mat3 angVelSkew = glm::matrixCross3(Spheres[i]->getAngVel());
				//create 3x3 rotation matrix from rb rotation matrix
				glm::mat3 R = glm::mat3(Spheres[i]->getRotate());
				//update rotation matrix
				R += dt * angVelSkew * R;
				R = glm::orthonormalize(R);
				Spheres[i]->setRotate(glm::mat4(R));

				
				if (currentTime > impAppTime && !impulseAppplied) {
					if (triangle) {
						if (i == 1) {		//We only apply the impulse to the white ball
							if (mue>0.6){ impulseF = vec3(16.0, 0.0f, -8.0f); }		//Depending on mue, we want to test different kind of impulses for the demo, thus they are changed
							else{ impulseF = vec3(0, 0.0f, -8.0f); }
							
							Vertex inapppoint = (Spheres[i]->getPos() + vec3(0.0f, 0.01f, 0.0f)); //initial application point
							Spheres[i]->setVel(Spheres[i]->getVel() + (impulseF) / Spheres[i]->getMass());
							Spheres[i]->setAngVel(Spheres[i]->getAngVel() + Spheres[i]->getItinverse()*glm::cross(inapppoint.getCoord() - Spheres[i]->getPos(), impulseF));
							impulsenumcheck = impulsenumcheck + 1;
							if (impulsenumcheck >= spherenum) {
								impulseAppplied = true;
							}
						}
					}
					else {
						if (i < spherenum / 3) {	//For random positions, we apply the impulses to a third of the ball number
							impulseF = vec3(RandomFloat(10.0f, 15.0f), 0.0f, RandomFloat(10.0f, 15.0f));
							Vertex inapppoint = (Spheres[i]->getPos() + vec3(0.0f, 0.008f, 0.0f)); //initial application point
							Spheres[i]->setVel(Spheres[i]->getVel() + (impulseF) / Spheres[i]->getMass());
							Spheres[i]->setAngVel(Spheres[i]->getAngVel() + Spheres[i]->getItinverse()*glm::cross(inapppoint.getCoord() - Spheres[i]->getPos(), impulseF));
							impulsenumcheck = impulsenumcheck + 1;
							if (impulsenumcheck >= spherenum) {
								impulseAppplied = true;
							}
						}
					}
				}


				//PLane Collsion with When Falling(Will need to chanage to bal;l collision)
				std::vector<Vertex> colVert = giveBallColVertices(plane.getPos()[1], Spheres[i]);

				vec3 sumpoint = vec3(0.0f);

				float planeBallDistance = Spheres[i]->getPos().y - plane.getPos().y;

				if (planeBallDistance<=1.0f) {
					sumpoint = Spheres[i]->getPos()-vec3(0.0f,Spheres[i]->getRadius(),0.0f);		//This will be the point of impulse
				}
				if (sumpoint != vec3(0.0f)) {

					vec3 pointofimpulse = /*Vertex(sumpoint / colVert.size())*/ sumpoint;		//We don't calculate lowest verices, as balls don't bounce, so the lowest point will be position - radius




				//FRICTION


				//Calculate Friction
					vec3 r = pointofimpulse - Spheres[i]->getPos(); // r in j formula
					vec3 vr = Spheres[i]->getVel() + glm::cross(Spheres[i]->getAngVel(), r);
					vec3 n = glm::normalize(vec3(0.0f, 1.0f, 0.0f));
					float e = 0.8f;
					vec3 vt = vr - dot(vr, n)*n;
					//n = normalize(Spheres[i]->getVel())*-1.0f;
					vec3 jFriction;
					vec3 tan = normalize(vt);

					if (vt != vec3(0.0f)) {

						//jFriction = -mue * glm::length(jn) * normalize(vt);	//Attempt 1, jn was removed as there was no bouncing
						jFriction = (-mue * tan) / ((1.0f + Spheres[i]->getMass()) + dot(glm::cross(Spheres[i]->getItinverse()*cross(r, tan), r), tan)); //WORKING FRICTION
						//jFriction = -tan * (mue*dot(vr, tan)) / ((1.0f + Spheres[i]->getMass()) + dot(glm::cross(Spheres[i]->getItinverse()*cross(r, tan), r), tan));
						//jFriction = -1.0f * mue * glm::length(jn)*(vt / glm::length(vt));  //2nd implementation, MUE MUST BE SMALLER THAN 0.2


					}
					else {
						jFriction = vec3(0.0f);
					}

					float minus = -1.0f;

					//Set VEL and ANG.VEl
					if (glm::length(Spheres[i]->getVel()) < 10.0f) {
						Spheres[i]->setVel(Spheres[i]->getVel()*0.99f + (jFriction / Spheres[i]->getMass()));
						AccVel = true;
						Spheres[i]->setAngVel((Spheres[i]->getAngVel()*0.95 + Spheres[i]->getItinverse()*glm::cross(r, jFriction)));
						if (mue > 0.6) {
							if (glm::length(Spheres[i]->getAngVel()) < 0.16) {	//Balls stopped when mue was high, but the were moving very liitle bit  constantly after they stopped, thus this was implemented
								Spheres[i]->setVel(vec3(0.0f));
							}
						}
						
					}
					else {
						Spheres[i]->setVel(Spheres[i]->getVel() + (jFriction / Spheres[i]->getMass()));
					}

					if (currentTime < impAppTime) {
						Spheres[i]->setAngVel((Spheres[i]->getAngVel() + Spheres[i]->getItinverse()*glm::cross(r, jFriction)));

						//create skew symmetric matrix for w
						angVelSkew = glm::matrixCross3(Spheres[i]->getAngVel());
						//create 3x3 rotation matrix from rb rotation matrix
						R = glm::mat3(Spheres[i]->getRotate());
						//update rotation matrix
						R += dt * angVelSkew * R;
						R = glm::orthonormalize(R);
						Spheres[i]->setRotate(glm::mat4(R));
					}
					else {

						if (!AccVel) {
							Spheres[i]->setAngVel((Spheres[i]->getAngVel() + Spheres[i]->getItinverse()*glm::cross(r, jFriction)));
						}

							//create skew symmetric matrix for w
							angVelSkew = glm::matrixCross3(Spheres[i]->getAngVel());
							//create 3x3 rotation matrix from rb rotation matrix
							R = glm::mat3(Spheres[i]->getRotate());
							//update rotation matrix
							R += dt * angVelSkew * R;
							R = glm::orthonormalize(R);
							Spheres[i]->setRotate(glm::mat4(R));
					}

				}

					//Collisions with cushion UNOPTIMISED
					for (unsigned int k = 0; k < 3; k++) {
						if (k != 1) {

							if (Spheres[i]->getPos()[k] < plane.getPos()[k] - plane.getScale()[k][k] + Spheres[i]->getRadius()) {
								//Spheres[i]->setVel(glm::vec3(0.0f));
								Spheres[i]->setVel(k, Spheres[i]->getVel()[k] * -0.9f);
								Spheres[i]->setAngVel(Spheres[i]->getAngVel()*0.8f);	//Balls will lose a bit of Angular Velocity and speed when bouncing of the wall to make it look more realistic
								Spheres[i]->setPos(k, plane.getPos()[k] - plane.getScale()[k][k] + Spheres[i]->getRadius());
							}

							else if (Spheres[i]->getPos()[k] > plane.getPos()[k] + plane.getScale()[k][k] - Spheres[i]->getRadius()) {
								//Spheres[i]->setVel(glm::vec3(0.0f));
								Spheres[i]->setVel(k, Spheres[i]->getVel()[k] * -0.9f);
								Spheres[i]->setAngVel(Spheres[i]->getAngVel()*0.8f);		//Balls will lose a bit of Angular Velocity and speed when bouncing of the wall to make it look more realistic
								Spheres[i]->setPos(k, plane.getPos()[k] + plane.getScale()[k][k] - Spheres[i]->getRadius());

							}


							//create skew symmetric matrix for w
							angVelSkew = glm::matrixCross3(Spheres[i]->getAngVel());
							//create 3x3 rotation matrix from rb rotation matrix
							R = glm::mat3(Spheres[i]->getRotate());
							//update rotation matrix
							R += dt * angVelSkew * R;
							R = glm::orthonormalize(R);
							Spheres[i]->setRotate(glm::mat4(R));
						}

					}


					for (unsigned int j = i + 1; j < spherenum; j++) {
						if (glm::distance(Spheres[i]->getPos(), Spheres[j]->getPos()) < Spheres[i]->getRadius() + Spheres[j]->getRadius()) {
							vec3 n = glm::normalize(Spheres[j]->getPos() - Spheres[i]->getPos());

							float displacement = Spheres[j]->getRadius() + Spheres[i]->getRadius() - glm::distance(Spheres[j]->getPos(), Spheres[i]->getPos());

							float ball1portion = glm::length(Spheres[i]->getVel()) / (glm::length(Spheres[i]->getVel()) + glm::length(Spheres[j]->getVel()));
							float ball2portion = glm::length(Spheres[j]->getVel()) / (glm::length(Spheres[i]->getVel()) + glm::length(Spheres[j]->getVel()));

							Spheres[i]->translate(displacement * -n * ball1portion);
							Spheres[j]->translate(displacement * n * ball2portion);

							vec3 vr = Spheres[j]->getVel() - Spheres[i]->getVel();
							//vec3 vr = Spheres[i]->getVel() - Spheres[j]->getVel();

							float jn = (-(1.0f + e)*dot(vr, n)) / (1.0f / Spheres[i]->getMass() + 1.0f / Spheres[j]->getMass());

							if (Spheres[i]->getVel() == vec3(0.0f)) {
								Spheres[i]->setAngVel(Spheres[j]->getAngVel());
							}
							else if (Spheres[j]->getVel() == vec3(0.0f)) {
								Spheres[j]->setAngVel(Spheres[i]->getAngVel());
							}

							//create skew symmetric matrix for w
							angVelSkew = glm::matrixCross3(Spheres[i]->getAngVel());
							//create 3x3 rotation matrix from rb rotation matrix
							R = glm::mat3(Spheres[i]->getRotate());
							//update rotation matrix
							R += dt * angVelSkew * R;
							R = glm::orthonormalize(R);
							Spheres[i]->setRotate(glm::mat4(R));

							Spheres[i]->setVel(Spheres[i]->getVel() - (jn*n / Spheres[i]->getMass()));
							Spheres[j]->setVel(Spheres[j]->getVel() + (jn*n / Spheres[j]->getMass()));

						}
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

