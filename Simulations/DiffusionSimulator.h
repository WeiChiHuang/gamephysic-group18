#ifndef DIFFUSIONSIMULATOR_h
#define DIFFUSIONSIMULATOR_h

#include "Simulator.h"
#include "vectorbase.h"

class Grid {
public:
	// Construtors
	Grid();
	Grid(size_t const m, size_t const n, size_t const d);
	// Functions
	size_t get_m() const ;
	size_t get_n() const;
	size_t get_d() const;
private:
	// Attributes
	size_t m_ = 0;
	size_t n_ = 0;
	size_t d_ = 0;
	// to be implemented
};



class DiffusionSimulator:public Simulator{
public:
	// Construtors
	DiffusionSimulator();

	// Functions
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
	void simulateTimestep(float timeStep);
	void externalForcesCalculations(float timeElapsed) {};
	void onClick(int x, int y);
	void onMouse(int x, int y);
	// Specific Functions
	void drawObjects();

	// Feel free to change the signature of these functions, add arguments, etc.
	void diffuseTemperatureExplicit(float timeStep);
	void diffuseTemperatureImplicit(float timeStep);

	void diffuseTemperatureImplicit2D(float timeStep);

private:
	// Attributes
	Vec3  m_vfMovableObjectPos;
	Vec3  m_vfMovableObjectFinalPos;
	Vec3  m_vfRotate;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
	Grid T;
	vector<vector<vector<double>>> temp_vector_;
	size_t x_ = 16;
	size_t y_ = 16;
	size_t z_ = 16;
	double alpha_ = 0.3;
	double boundary_temperature_ = 0.0;
	double max_temperature_ = 100.0;
};

#endif