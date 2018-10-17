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
	//complete. should return the force from gravity		
	glm::vec3 ForceG = mass * (glm::vec3(0.0f, -9.8f, 0.0f));
	//glm::vec3 acceleration = ForceG / mass; //
	return ForceG;
}

/*
** Drag
*/
glm::vec3 Drag::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	//complete. should return the force from aerodynamic Drag
	float absoluteu = glm::length(vel);
	glm::vec3 e = -vel / absoluteu;
	glm::vec3 Fa = 0.5*density*absoluteu*absoluteu*aircoefficient*area*e;
	//glm::vec3 acceleration = Fa / mass; // 
	return Fa;
}

/*
** Hooke's Law
*/
glm::vec3 Hooke::apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel) {
	//complete. should return the force from the Spring force

	float length = glm::length(m_b2->getPos() - m_b1->getPos());
	glm::vec3 e;
	if (length == 0) {
		e = glm::vec3(0.0f);
	}
	else {
		e = (m_b2->getPos() - m_b1->getPos()) / length;
	}
	


	float Fh = (-1.0f) * getks() * (getrest() - length); //x in this case is rest lenth -current length


	float Fd = (-1.0f)*getkd() * (dot(m_b1->getVel(), e) - dot(m_b2->getVel(), e));
	//glm::vec3 acceleration = (glm::vec3(0.0f, (Fh+Fd), 0.0f))/mass;

	if (m_b2->getPos() == pos) {
		return -(Fh+Fd)*e;
	}
	else {
		return (Fh+Fd)*e;
	}
	
	//return glm::vec3(0.0f, (Fh + Fd), 0.0f);
}