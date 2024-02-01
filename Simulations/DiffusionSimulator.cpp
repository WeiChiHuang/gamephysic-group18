#include "DiffusionSimulator.h"
#include "pcgsolver.h"
using namespace std;


DiffusionSimulator::DiffusionSimulator()
{
	m_iTestCase = 0;
	m_vfMovableObjectPos = Vec3();
	m_vfMovableObjectFinalPos = Vec3();
	m_vfRotate = Vec3();
	T = Grid(x_, y_, z_);
	
	// rest to be implemented
}

const char * DiffusionSimulator::getTestCasesStr(){
	return "Explicit_solver, Implicit_solver";
}

void DiffusionSimulator::reset(){
		m_mouse.x = m_mouse.y = 0;
		m_trackmouse.x = m_trackmouse.y = 0;
		m_oldtrackmouse.x = m_oldtrackmouse.y = 0;
		temp_vector_.clear();
		T = Grid(x_, y_,z_);
		temp_vector_.resize(T.get_m(), std::vector<std::vector<double>>(T.get_n(), std::vector<double>(T.get_d(), 0.0)));
		for (int i = T.get_m()-1; i  > T.get_m()/2; i--) {
			for (int j = T.get_n() - 1; j > T.get_n()/2; j--) {
				if (T.get_d() > 2)
				{
					for (int d = 0; d < T.get_d()/2 +1 ; d++) {
						temp_vector_[i][j][d] = 80.0;
					}
				}else if (T.get_d() == 2)
				{
					temp_vector_[i][j][0] = 80.0;
					temp_vector_[i][j][1] = 80.0;
				}
				else {
					temp_vector_[i][j][0] = 80.0;
				}
			}
		}

}

void DiffusionSimulator::initUI(DrawingUtilitiesClass * DUC)
{
	this->DUC = DUC;
	TwAddVarRW(DUC->g_pTweakBar, "nx", TW_TYPE_INT32, &x_, "min=1 max=100 step=1");
	TwAddVarRW(DUC->g_pTweakBar, "ny", TW_TYPE_INT32, &y_, "min=1 max=100 step=1");
	TwAddVarRW(DUC->g_pTweakBar, "nz", TW_TYPE_INT32, &z_, "min=1 max=100 step=1");
	TwAddVarRW(DUC->g_pTweakBar, "alpha", TW_TYPE_DOUBLE, &alpha_, "min=0 max=1 step=0.1");
	// to be implemented
}

void DiffusionSimulator::notifyCaseChanged(int testCase)
{
	m_iTestCase = testCase;
	reset();
	m_vfMovableObjectPos = Vec3(0, 0, 0);
	m_vfRotate = Vec3(0, 0, 0);
	m_vfMovableObjectPos = m_vfMovableObjectPos - Vec3((T.get_m()-1) * 0.05, (T.get_n()-1) * 0.05, 0);
	//
	// to be implemented
	//
	switch (m_iTestCase)
	{
	case 0:
		cout << "Explicit solver!\n";
		break;
	case 1:
		cout << "Implicit solver!\n";
		break;
	default:
		cout << "Empty Test!\n";
		break;
	}
}

void DiffusionSimulator::diffuseTemperatureExplicit(float timeStep) {
	vector<vector<vector<double>>> newTemperature;
	double initialValue = 0.0; 
	newTemperature.resize(T.get_m(), std::vector<std::vector<double>>(T.get_n(), std::vector<double>(T.get_d(), initialValue)));
		for (int i = 1; i < T.get_m() - 1; ++i) {
			for (int j = 1; j < T.get_n() - 1; ++j) {
				if(T.get_d() > 2)
				{
					for (int d = 1; d < T.get_d() - 1; ++d) {
						newTemperature[i][j][d] = temp_vector_[i][j][d] + alpha_ * timeStep * (
							(temp_vector_[i + 1][j][d] + temp_vector_[i - 1][j][d] - 2 * temp_vector_[i][j][d]) +
							(temp_vector_[i][j + 1][d] + temp_vector_[i][j - 1][d] - 2 * temp_vector_[i][j][d]) +
							(temp_vector_[i][j][d + 1] + temp_vector_[i][j][d - 1] - 2 * temp_vector_[i][j][d])
							);
					}
				}
				else if(T.get_d() == 2)
				{
					newTemperature[i][j][0] = temp_vector_[i][j][0] + alpha_ * timeStep * (
						(temp_vector_[i + 1][j][0] + temp_vector_[i - 1][j][0] - 2 * temp_vector_[i][j][0]) +
						(temp_vector_[i][j + 1][0] + temp_vector_[i][j - 1][0] - 2 * temp_vector_[i][j][0])
						);
					newTemperature[i][j][1] = temp_vector_[i][j][1] + alpha_ * timeStep * (
						(temp_vector_[i + 1][j][1] + temp_vector_[i - 1][j][1] - 2 * temp_vector_[i][j][1]) +
						(temp_vector_[i][j + 1][1] + temp_vector_[i][j - 1][1] - 2 * temp_vector_[i][j][1])
						);
				}
				else
				{
					newTemperature[i][j][0] = temp_vector_[i][j][0] + alpha_ * timeStep * (
						(temp_vector_[i + 1][j][0] + temp_vector_[i - 1][j][0] - 2 * temp_vector_[i][j][0]) +
						(temp_vector_[i][j + 1][0] + temp_vector_[i][j - 1][0] - 2 * temp_vector_[i][j][0])
						);
				}
			}
		}

	 ////bound
		////front back x y
		//for (int i = 0; i < T.get_m(); ++i) {
		//	for (int j = 0; j < T.get_n(); ++j) {
		//		if (j == 0 || j == T.get_n() - 1 || i == 0 || i == T.get_m() - 1) {
		//			newTemperature[i][j][0] = boundary_temperature_;
		//			newTemperature[i][j][T.get_d() - 1] = boundary_temperature_;
		//		}
		//	}
		//}
		////left right y z
		//for (int j = 0; j < T.get_n(); ++j) {
		//	for (int k = 0; k < T.get_d()	; ++k) {
		//		if (j == 0 || j == T.get_n() - 1 || k == 0 || k == T.get_d() - 1) {
		//			newTemperature[0][j][k] = boundary_temperature_;
		//			newTemperature[T.get_m() - 1][j][k] = boundary_temperature_;
		//		}
		//	}
		//}
		//// down top  x z
		//for (int i = 0; i < T.get_m(); ++i) {
		//	for (int k = 0; k < T.get_d(); ++k) {
		//		if (k == 0 || k == T.get_d() - 1 || i == 0 || i == T.get_m() - 1) {
		//		newTemperature[i][0][k] = boundary_temperature_; 
		//		newTemperature[i][T.get_n() - 1][k] = boundary_temperature_;
		//			}
		//	}
		//}
// to be implemented
		//update temperature
		for (int i = 0; i < T.get_m(); ++i) {
			for (int j = 0; j < T.get_n(); ++j) {
				for (int d = 0; d < T.get_d(); ++d) {
					temp_vector_[i][j][d] = newTemperature[i][j][d];
				}
			}
		}
}


void DiffusionSimulator::diffuseTemperatureImplicit(float timeStep) {
	// solve A T = b

	// Flatten the 2D grid into a 1D vector

	// This is just an example to show how to work with the PCG solver,
	/*const int nx = 5;
	const int ny = 5;
	const int nz = 5;*/
	const int N = T.get_m()*T.get_n()*T.get_d();

	SparseMatrix<Real> A(N);
	A.zero();
	std::vector<Real> b(N);
	for (int i = 0; i < b.size(); i++) {
		b[i] = 0;
	}
	// This is the part where you have to assemble the system matrix A and the right-hand side b!
	for (int i = 0; i < T.get_m(); ++i) {
		for (int j = 0; j < T.get_n(); ++j) {
			for (int k = 0; k < T.get_d(); ++k) {
				int index = i * T.get_n() * T.get_d() + j * T.get_d() + k;
				if (i == 0 || j == 0 || k == 0 || (i == T.get_m() - 1) || (j == T.get_n() - 1) || (k == T.get_d() - 1)) {
					A.set_element(index, index, 1);
				}
				else
				{
					double center = 1 + 6 * alpha_ * timeStep;
					double left = -alpha_ * timeStep;
					double right = -alpha_ * timeStep;
					double top = -alpha_ * timeStep;
					double down = -alpha_ * timeStep;
					double front = -alpha_ * timeStep;
					double back = -alpha_ * timeStep;
					A.set_element(index, index, center);
					A.set_element(index - 1, index, left);
					A.set_element(index + 1, index, right);
					A.set_element(index - T.get_d(), index, top);
					A.set_element(index + T.get_d(), index, down);
					A.set_element(index - T.get_n() * T.get_d(), index, front);
					A.set_element(index + T.get_n() * T.get_d(), index, back);
				}
			}
		}
	}

	for (int i = 0; i < T.get_m(); ++i) {
		for (int j = 0; j < T.get_n(); ++j) {
			for (int k = 0; k < T.get_d(); ++k) {
				b[i * T.get_n() * T.get_d() + j * T.get_d()+k] = temp_vector_[i][j][k];
			}
		}
	}

	// perform solve
	Real pcg_target_residual = 1e-05;
	Real pcg_max_iterations = 1000;
	Real ret_pcg_residual = 1e10;
	int  ret_pcg_iterations = -1;

	SparsePCGSolver<Real> solver;
	solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97, 0.25);

	std::vector<Real> x(N,0.0);
	for (int j = 0; j < N; ++j) { x[j] = 0.; }

	// preconditioners: 0 off, 1 diagonal, 2 incomplete cholesky
	solver.solve(A, b, x, ret_pcg_residual, ret_pcg_iterations, 0);

	// Final step is to extract the grid temperatures from the solution vector x
	// to be implemented
	int idx = 0;
	for (int i = 0; i < T.get_m(); ++i) {
		for (int j = 0; j < T.get_n(); ++j) {
			for (int k = 0; k < T.get_d(); ++k) {
				temp_vector_[i][j][k] = x[idx++];
			}
		}
	}

	//front back x y
	for (int i = 0; i < T.get_m(); ++i) {
		for (int j = 0; j < T.get_n(); ++j) {
				temp_vector_[i][j][0] = boundary_temperature_;
				temp_vector_[i][j][T.get_d() - 1] = boundary_temperature_;
		}
	}
	//left right y z
	for (int j = 0; j < T.get_n(); ++j) {
		for (int k = 0; k < T.get_d()	; ++k) {
				temp_vector_[0][j][k] = boundary_temperature_;
				temp_vector_[T.get_m() - 1][j][k] = boundary_temperature_;
		}
	}
	// down top  x z
	for (int i = 0; i < T.get_m(); ++i) {
		for (int k = 0; k < T.get_d(); ++k) {
				temp_vector_[i][0][k] = boundary_temperature_;
				temp_vector_[i][T.get_n() - 1][k] = boundary_temperature_;
				}
	}
}

void DiffusionSimulator::diffuseTemperatureImplicit2D(float timeStep) {
	// solve A T = b

	// Flatten the 2D grid into a 1D vector

	// This is just an example to show how to work with the PCG solver,
	/*const int nx = 5;
	const int ny = 5;
	const int nz = 5;*/
	const int N = T.get_m() * T.get_n();

	SparseMatrix<Real> A(N);
	A.zero();
	std::vector<Real> b(N);
	for (int i = 0; i < b.size(); i++) {
		b[i] = 0;
	}
	// This is the part where you have to assemble the system matrix A and the right-hand side b!
	for (int i = 0; i < T.get_m(); ++i) {
		for (int j = 0; j < T.get_n(); ++j) {
				int index = i * T.get_n() * T.get_d() + j * T.get_d() ;
				if (i == 0 || j == 0 || (i == T.get_m() - 1) || (j == T.get_n() - 1)) {
					A.set_element(index, index, 1);
				}
				else
				{
					double center = 1 + 4 * alpha_ * timeStep;
					double left = -alpha_ * timeStep;
					double right = -alpha_ * timeStep;
					double top = -alpha_ * timeStep;
					double down = -alpha_ * timeStep;
					A.set_element(index, index, center);
					A.set_element(index - 1, index, left);
					A.set_element(index + 1, index, right);
					A.set_element(index - T.get_m(), index, top);
					A.set_element(index + T.get_m(), index, down);
				}
			}
		}
	for (int i = 0; i < T.get_m(); ++i) {
		for (int j = 0; j < T.get_n(); ++j) {
				b[i * T.get_n() * T.get_d() + j * T.get_d()] = temp_vector_[i][j][0];
		}
	}

	// perform solve
	Real pcg_target_residual = 1e-05;
	Real pcg_max_iterations = 1000;
	Real ret_pcg_residual = 1e10;
	int  ret_pcg_iterations = -1;

	SparsePCGSolver<Real> solver;
	solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97, 0.25);

	std::vector<Real> x(N, 0.0);
	for (int j = 0; j < N; ++j) { x[j] = 0.; }

	// preconditioners: 0 off, 1 diagonal, 2 incomplete cholesky
	solver.solve(A, b, x, ret_pcg_residual, ret_pcg_iterations, 0);

	// Final step is to extract the grid temperatures from the solution vector x
	// to be implemented
	int idx = 0;
	for (int i = 0; i < T.get_m(); ++i) {
		for (int j = 0; j < T.get_n(); ++j) {
				temp_vector_[i][j][0] = x[idx++];
		}
	}
	//left right y z
	for (int j = 0; j < T.get_n(); ++j) {
		for (int k = 0; k < T.get_d()	; ++k) {
				temp_vector_[0][j][0] = boundary_temperature_;
				temp_vector_[T.get_m() - 1][j][0] = boundary_temperature_;
		}
	}
	// down top  x z
	for (int i = 0; i < T.get_m(); ++i) {
		for (int k = 0; k < T.get_d(); ++k) {
				temp_vector_[i][0][0] = boundary_temperature_;
				temp_vector_[i][T.get_n() - 1][0] = boundary_temperature_;
				}
	}
}


void DiffusionSimulator::simulateTimestep(float timeStep)
{
	// update current setup for each frame
	switch (m_iTestCase)
	{
	case 0:
		// feel free to change the signature of this function
		diffuseTemperatureExplicit(timeStep);
		break;
	case 1:
		// feel free to change the signature of this function
		if(T.get_d() > 1)
		{
			diffuseTemperatureImplicit(timeStep);
		}
		else
		{
			diffuseTemperatureImplicit2D(timeStep);
		}
		break;
	}
}

void DiffusionSimulator::drawObjects()
{
	std::mt19937 eng;
	std::uniform_real_distribution<float> randCol(0.0f, 1.0f);
	std::uniform_real_distribution<float> randPos(-0.5f, 0.5f);
	for(size_t x = 0; x < T.get_m(); x++)
	{
		Vec3 pos = m_vfMovableObjectPos + Vec3(float(0.1*x), 0.0, 0.0);
		for(size_t y = 0; y < T.get_n(); y++)
		{
			Vec3 pos_y = pos + Vec3(0.0, float(0.1 * y), 0.0);
			for (size_t z = 0; z < T.get_d(); z++)
			{
				pos_y = pos_y + Vec3(0.0, 0.0, 0.1);
				Vec3 color;
				double temperature = temp_vector_[x][y][z];
				if (temperature <= 0.0) {
					color = 0.6 * Vec3(0, 0, 0);
				}
				else {
					double normalizedTemp = temperature / max_temperature_;
					color = normalizedTemp * Vec3(1,1,1);
				}
				if(T.get_d() <= 2){
					DUC->setUpLighting(Vec3(), 0.4 * Vec3(1, 1,1), 100, color);
					DUC->drawSphere(pos_y, Vec3(0.0625));
				}else
				{
					if(z!=0)
					{
						DUC->setUpLighting(Vec3(), 0.4 * Vec3(1, 1, 1), 100, color);
						DUC->drawSphere(pos_y, Vec3(0.0425));
					}
				
				}
			}
		}
	}
	// to be implemented
	//visualization
}


void DiffusionSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext)
{
	drawObjects();
}

void DiffusionSimulator::onClick(int x, int y)
{
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

void DiffusionSimulator::onMouse(int x, int y)
{
	m_oldtrackmouse.x = x;
	m_oldtrackmouse.y = y;
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

Grid::Grid() :m_(16), n_(16),d_(16)
{

}

Grid::Grid(size_t const m, size_t const n, size_t const d):m_(m),n_(n),d_(d)
{

}

size_t Grid::get_m() const
{
	return this->m_;
}

size_t Grid::get_n() const
{
	return this->n_;
}

size_t Grid::get_d() const
{
	return this->d_;
}
