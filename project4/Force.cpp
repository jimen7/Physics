#include <iostream>
#include <cmath>
#include "Force.h"
#include "Body.h"
#include "glm/ext.hpp"

float density = 1.225f;
float aircoefficient = 1.05f;
glm::vec3 area = glm::vec3(0.1f, 0.1f, 0.0f);
//float ks = 0.25f; //Spring constant, how loose the sppring is
//float kd = 1.0f; //Damping Factor
//float springrest = 1.5f;


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
	glm::vec3 Fa = 0.5*density*absoluteu*absoluteu*aircoefficient*area*e;
	glm::vec3 acceleration = Fa / mass; // 
	return acceleration;
}

/*
** Hooke's Law
*/
glm::vec3 Hooke::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	//complete. should return the acceleration from the Spring force

	float Fh = (-1.0f)* (this->getks()) * (this->getrest() - /*pos[1]*/glm::length(pos)); //x in this case is rest lenth -current length
	float Fd = (-1.0f)*(this->getkd()) * ( - glm::length(vel));
	glm::vec3 acceleration = (glm::vec3(0.0f, Fh+Fd, 0.0f))/mass;
	return acceleration;
}