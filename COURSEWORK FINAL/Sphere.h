#pragma once
#include "RigidBody.h"

class Sphere:
	public RigidBody
{
public:
	Sphere();
	~Sphere();

	void setRadius(const float &r) { m_radius = r*glm::normalize(getScale()[1][1]);  }

	float getRadius() { return m_radius; }


private:
	float m_radius;
};

