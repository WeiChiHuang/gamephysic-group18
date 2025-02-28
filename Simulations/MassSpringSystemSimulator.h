#ifndef MASSSPRINGSYSTEMSIMULATOR_h
#define MASSSPRINGSYSTEMSIMULATOR_h
#include "Simulator.h"

// Do Not Change
#define EULER 0
#define LEAPFROG 1
#define MIDPOINT 2
// Do Not Change


class MassSpringSystemSimulator:public Simulator{
public:
	// Construtors
	MassSpringSystemSimulator();
	
	// UI Functions
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
	void externalForcesCalculations(float timeElapsed);
	void simulateTimestep(float timeStep);
	void onClick(int x, int y);
	void onMouse(int x, int y);

	// Specific Functions
	void setMass(float mass);
	void setStiffness(float stiffness);
	void setDampingFactor(float damping);
	int addMassPoint(Vec3 position, Vec3 Velocity, bool isFixed);
	void addSpring(int masspoint1, int masspoint2, float initialLength);
	int getNumberOfMassPoints();
	int getNumberOfSprings();
	Vec3 getPositionOfMassPoint(int index);
	Vec3 getVelocityOfMassPoint(int index);
	void applyExternalForce(Vec3 force);

	// Do Not Change
	void setIntegrator(int integrator) {
		m_iIntegrator = integrator;
	}

private:
	// Data Attributes
	float m_fMass;
	float m_fStiffness;
	float m_fDamping;
	int m_iIntegrator;

	// Extra Attributes
	struct Masspoint {
		Vec3 position;
		Vec3 Velocity;
		bool isFixed;
	};
	std::vector<Masspoint> masspointVector;
	struct Spring {
		int masspoint1;
		int masspoint2;
		float initialLength;
	};
	std::vector<Spring> springVector;
	float massPointSize;
	Vec3 lineColor;
	float defalut_initlength;
	float m_fGravity;
	float m_fTimestep;
	int m_numberOfSpring;

	// UI Attributes
	Vec3 m_externalForce;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;

	// Extra Functions
	void explicitEuler(float timeStep);
	std::tuple <Vec3, Vec3>  acceleration(Masspoint mp1, Masspoint mp2);
	Vec3 direction(Vec3 p1, Vec3 p2);
	void midpoint(float timeStep);
	void leapfrog(float timeStep);
	void wallCollision(Masspoint& mp);
	
};
#endif