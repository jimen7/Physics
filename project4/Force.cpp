#include <iostream>
#include <cmath>
#include "Force.h"
#include "Body.h"
#include "glm/ext.hpp"

float density = 1.225f;
float coefficient = 1.05f;
glm::vec3 area = glm::vec3(0.1f, 0.1f, 0.0f);

Force::Force()
{
}


Force::~Force()
{
}

glm::vec3 Force::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {				
	return glm::vec3(0.0f);
}

/*
** Gravity
*/
glm::vec3 Gravity::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	//complete. should return the acceleration from gravity		
	glm::vec3 ForceG = mass * (glm::vec3(0.0f, -9.8f, 0.0f));
	glm::vec3 acceleration = ForceG / mass; //
	return acceleration;
}

/*
** Drag
*/
glm::vec3 Drag::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	//complete. should return the acceleration from aerodynamic Drag
	float absoluteu = glm::length(vel); 
	glm::vec3 e = -vel / absoluteu;
	glm::vec3 Fa = 0.5*density*absoluteu*absoluteu*coefficient*area*e;
	glm::vec3 acceleration = Fa / mass; // 
	return acceleration;
}

/*
** Hooke's Law
*/
glm::vec3 Hooke::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	//complete. should return the acceleration from the Spring force
	
	return glm::vec3(0.0f);
}