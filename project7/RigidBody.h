#pragma once
#include "Body.h"

class RigidBody :
	public Body
{
public:
	RigidBody();
	~RigidBody();

	//set and get methods
	void setAngVel(const glm::vec3 &omega) { m_angVel = omega; }
	void setAngAccl(const glm::vec3 &alpha) { m_angAcc = alpha; }
	void setAngVel(const glm::mat3 &invIneria) { m_invInertia = invIneria; }

	void scale(const glm::vec3 &vect) override;

	void calculateInertia();


	glm::vec3 getAngVel() { return m_angVel; }
	glm::vec3 getAngAcc() { return m_angAcc; }
	glm::mat3 getInvInertia() { return m_invInertia; }
	//void getScale();
	//glm::mat3 



private:
	float m_density;
	glm::mat3 m_invInertia; //Inverse Inertia
	glm::vec3 m_angVel; //Angular velocity
	glm::vec3 m_angAcc; //Angular Acceleration
};

