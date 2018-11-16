#include "RigidBody.h"




RigidBody::RigidBody()
{
}


RigidBody::~RigidBody()
{
}

void RigidBody::calculateInertia() {
	glm::mat3 tensor = glm::mat3(0.0f);
	//float w = glm::distance(getMesh().getVertices()[1].getCoord(), getMesh().getVertices()[0].getCoord()) * getScale()[1][1];
	//float h = glm::distance(getMesh().getVertices()[1].getCoord(), getMesh().getVertices()[3].getCoord()) * getScale()[2][2];
	//float d = glm::distance(getMesh().getVertices()[1].getCoord(), getMesh().getVertices()[5].getCoord()) * getScale()[3][3];
	float w = getMesh().getScale()[1][1] * 2.0f;
	float h = getMesh().getScale()[2][2] * 2.0f;
	float d = getMesh().getScale()[3][3] * 2.0f;
	float m = getMass();
	tensor[0][0] = (1.0f / 12.0f)*m*(h * h + d * d);
	tensor[1][1] = (1.0f / 12.0f)*m*(w * w + d * d);
	tensor[2][2] = (1.0f / 12.0f)*m*(w * w + h * h);
	m_invInertia = glm::inverse(tensor);
	//return tensor;
}

void RigidBody::scale(const glm::vec3 &vect) {
	//__super::scale(vect);
	Body::scale(vect);
	this->calculateInertia();
}
