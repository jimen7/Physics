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


const int spherenum = 10;

const float sphereradius = 1.0f;


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



std::vector<Vertex> giveColVertices(float  y, RigidBody &r) {

	std::vector<Vertex> colVert;

	//CHecking if any vertice is added to the plane, and if it is add it to the vector of vrtices
	for (Vertex v : r.getMesh().getVertices()) {
		//World position of given vertex:
		vec3 worlpos = mat3(r.getMesh().getModel()) * v.getCoord() + r.getPos();
		//Check if the world position of the vertices are coplliding with the plane
		if (worlpos[1]<y) {
			colVert.push_back(worlpos);
		}
	}

	return colVert;
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
	//std::cout << "Size of plane: " << plane.getVertices().size() << std::endl;
	// scale it up x5
	plane.scale(glm::vec3(30.0f, 0.0f, 30.0f));
	Shader lambert = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	plane.setShader(lambert);

	Shader blue = Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag");





	//Set up a cubic rigid body
	RigidBody rb = RigidBody();
	Mesh m = Mesh::Mesh(Mesh::CUBE);
	rb.setMesh(m);
	Shader rbShader = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	rb.getMesh().setShader(rbShader);
	rb.scale(vec3(1.0f,3.0f,1.0f));
	rb.setMass(1.0f);

	//Set up spheres
	std::vector<Sphere> Spheres;
	Mesh sphere = Mesh::Mesh("resources/models/sphere.obj");
	for (unsigned int i = 0; i < spherenum; i++){
		Sphere sp = Sphere();

		sp.setMesh(sphere);
		sp.getMesh().setShader(rbShader);

		//sp.scale(glm::vec3(sphereradius, sphereradius, sphereradius));
		sp.setRadius(1.0f);
		sp.setMass(1.0f);

		Spheres.push_back(sp);


		
	}

	for (unsigned int i = 0; i < spherenum; i++) {
		Spheres[i].translate(glm::vec3(RandomFloat(1.0f, 29.0f), Spheres[i].getRadius(), RandomFloat(1.0f, 29.0f)));
		Spheres[i].setVel(vec3(20.0f, 0.0f, 20.0f));
		Spheres[i].setAngVel(vec3(0.0f, 0.0f, 0.0f));
		//spheres.addforce?
	}


	





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

	std::vector<glm::vec3> dir;
	dir.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
	dir.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
	dir.push_back(glm::vec3(0.0f, 0.0f, 1.0f));

	//Initial impulse
	bool impulseAppplied = false;
	vec3 impulseF = vec3(2.0f,1.0f,0.0f);
	//vec3 impulseF2 = vec3(3.0f, 0.0f, 0.0f);
	float impAppTime = 2.0f;

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


			//Translate Spheres
			for (unsigned int i = 0; i < spherenum; i++) {
				Spheres[i].setAcc(Spheres[i].applyForces(Spheres[i].getPos(), Spheres[i].getVel(), t, dt));
				Spheres[i].setVel(Spheres[i].getVel() + dt * Spheres[i].getAcc());
				Spheres[i].translate(dt*Spheres[i].getVel());
			}
			
			
			for (unsigned int i = 0; i < spherenum; i++) {
				if (Spheres[i].getPos().x > 30.0f || Spheres[i].getPos().z > 30.0f || Spheres[i].getPos().x < -30.0f || Spheres[i].getPos().z < -30.0f) {
					Spheres[i].setVel(-Spheres[i].getVel());
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
		for (unsigned int i = 0; i < spherenum; i++) {
			app.draw(Spheres[i].getMesh());
		}

		// draw demo objects
		//app.draw(cube);
		//app.draw(sphere);
		app.display();
	}

	app.terminate();

	return EXIT_SUCCESS;
}

