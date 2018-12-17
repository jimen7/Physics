#include "Sphere.h"



Sphere::Sphere()
{
}


Sphere::~Sphere()
{
}


void Sphere::calculateInertia() {
	glm::mat3 tensor = glm::mat3(0.0f);
	float r = getRadius();
	//float w = glm::distance(getMesh().getVertices()[1].getCoord(), getMesh().getVertices()[0].getCoord()) * getScale()[1][1];
	//float h = glm::distance(getMesh().getVertices()[1].getCoord(), getMesh().getVertices()[3].getCoord()) * getScale()[2][2];
	//float d = glm::distance(getMesh().getVertices()[1].getCoord(), getMesh().getVertices()[5].getCoord()) * getScale()[3][3];
	float w = getMesh().getScale()[1][1] * 2.0f;
	float h = getMesh().getScale()[2][2] * 2.0f;
	float d = getMesh().getScale()[3][3] * 2.0f;
	float m = getMass();
	tensor[0][0] = (2.0f / 3.0f)*m*r*r;
	tensor[1][1] = (2.0f / 3.0f)*m*r*r;
	tensor[2][2] = (2.0f / 3.0f)*m*r*r;
	m_invInertia = glm::inverse(tensor);
	//return tensor;
}