#pragma once
#include "RigidBody.h"

class Sphere :
	public RigidBody
{
public:
	Sphere();
	~Sphere();

	void setRadius(const float &r) { getMesh().scale(glm::vec3(r, r, r)); m_radius = r; }

	float getRadius() { return m_radius; }


private:
	float m_radius;
};

