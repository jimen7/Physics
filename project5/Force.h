#pragma once
#include <glm/glm.hpp>
#include <iostream>

class Body; //Forward declaration to avoid dependencies


class Force
{
public:
	Force();
	~Force();

	virtual glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);
};

/*
**Gravity Class
*/

class Gravity : public Force {

public:
	//constructors
	Gravity() {}
	Gravity(const glm::vec3 gravity) { m_gravity = gravity; }

	//get and set the methods
	glm::vec3 getGravity() const { return m_gravity; }
	void setGravity(glm::vec3 gravity) { m_gravity = gravity; }

	//physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:
	glm::vec3 m_gravity = glm::vec3(0.0f, 9.8f, 0.0f);
};

/*
**DRAG CLASS
*/
class Drag : public Force {
public:
	Drag() {}

	// physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:

};

class Hooke : public Force {
public:
	Hooke() {}
	Hooke(Body* b1, Body* b2, float ks, float kd, float rest) {
		m_ks = ks; m_kd = kd; m_rest = rest; m_b1 = b1; m_b2 = b2;
	}
	//get and set methods
	//you can write these yourself as necessary
	float getks() { return m_ks; }
	float getkd() { return m_kd; }
	float getrest() { return m_rest; }


	void setks(float & spconst) { m_ks = spconst; } //Spring constant
	void setkd(float & dampfac) { m_kd = dampfac; } //Damping factor
	void setrest(float & restpos) { m_rest = restpos; } //Rest legth



	//physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);


private:
	float m_ks; //spring stiffness
	float m_kd; //damping coefficient
	float m_rest; //spring rest length

	Body* m_b1; //pointer to the body connected to one extremity of the spring
	Body* m_b2; //pointer to the body connected to the other extremity
};


class Wind : public Force {
public:
	Wind() {}
	Wind(Body* b1, Body* b2, Body* b3) {
		m_b1 = b1; m_b2 = b2; m_b3 = b3;
	}

	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);


private:
	Body* m_b1; //pointer to the body connected to one extremity of the spring
	Body* m_b2; //pointer to the body connected to the other extremity
	Body* m_b3; //pointer to the body connected to the other extremity
};

