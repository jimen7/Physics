#include "RigidBody.h"




RigidBody::RigidBody()
{
}


RigidBody::~RigidBody()
{
}

void RigidBody::calculateInertia() {
	glm::mat3 tensor = glm::mat3();
	float w = glm::distance(getMesh().getVertices()[1], getMesh().getVertices()[1]);
	float h;
	float d;
	float m;
	tensor[0][0] = (1.0f / 12.0f)*m*(h * h + d * d);
	tensor[1][1] = (1.0f / 12.0f)*m*(w * w + d * d);
	tensor[2][2] = (1.0f / 12.0f)*m*(w * w + h * h);
	//return tensor;
}
