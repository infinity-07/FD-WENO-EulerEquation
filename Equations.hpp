#ifndef __EQUATIONS__HPP__
#define __EQUATIONS__HPP__

#include "scTools/scTools.h"

enum RIEMANNFLUXTYPE
{
	LLF = 0
};
enum SCHEMETYPE
{
	WENO = 0,
	WENOZ = 1,
};
enum TESTCASETYPE
{
	SMOOTH = 0,
	Riemann1 = 1,
	Riemann2 = 2,
	RTI = 3
};
enum BOUNDARYTYPE
{
	PERIOD = 0,
	slip = 1,
	DIRICHLET = 2,
	special = 3,
	symmetric = 4
};

// 方程

class EulerEquation
{
public:
	double xL, xR, yL, yR;
	double outputtime;
	BOUNDARYTYPE leftBoundaryCondition, rightBoundaryCondition;
	BOUNDARYTYPE topBoundaryCondition, bottomBoundaryCondition;

	std::function<double(double, double, double)> theVarExact;
	std::function<double(Array1D<double>)> theVarUh;
	bool u_exact_exist;

public:
	// ------------- (1) 方程形式与特征值，特征矩阵 -------------//
	// (1-0) 变量
	double euler_gamma;

	// (1-1) 方程数量
	inline int getVarNum()
	{
		return 4;
	}

	// (1-2) PhyFlux
	inline void getPhyFlux(const Array1D<double> Uh, Array1D<double> &Flux, double nx, double ny)
	{
		// U
		double rho = Uh[0];
		double u = Uh[1] / Uh[0];
		double v = Uh[2] / Uh[0];
		double E = Uh[3];
		double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));

		if (nx > ny)
		{
			// Flux
			Flux[0] = rho * u;
			Flux[1] = rho * u * u + p;
			Flux[2] = rho * u * v;
			Flux[3] = (E + p) * u;
		}
		else
		{
			Flux[0] = rho * v;
			Flux[1] = rho * u * v;
			Flux[2] = rho * v * v + p;
			Flux[3] = (E + p) * v;
		}
	}

	// (1-3) eigen Values
	inline double getMaxEigenValue(const Array1D<double> &Uh, double nx, double ny)
	{
		double rho = Uh[0];
		double u = Uh[1] / Uh[0];
		double v = Uh[2] / Uh[0];
		double E = Uh[3];
		double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		double c = sqrt(euler_gamma * p / rho);
		double MaxeigValue = fabs(u * nx + v * ny) + c;

		return MaxeigValue;
	}

	// (1-4) Left Eigen Matrix
	inline void getLEigenMatrix(const Array1D<double> &Uh, double nx, double ny, Array2D<double> &eigMatrix)
	{
		double rho = Uh[0];
		double u = Uh[1] / Uh[0];
		double v = Uh[2] / Uh[0];
		double E = Uh[3];
		double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		double c = sqrt(euler_gamma * p / rho);
		double B1 = (euler_gamma - 1) / pow(c, 2);
		double B2 = B1 * (pow(u, 2) + pow(v, 2)) / 2;

		if (fabs(nx) > fabs(ny)) // x-direction
		{
			eigMatrix[0][0] = (B2 + u / c) / 2;
			eigMatrix[1][0] = v;
			eigMatrix[2][0] = (1 - B2);
			eigMatrix[3][0] = (B2 - u / c) / 2;

			eigMatrix[0][1] = -(B1 * u + 1.0 / c) / 2;
			eigMatrix[1][1] = 0;
			eigMatrix[2][1] = B1 * u;
			eigMatrix[3][1] = -(B1 * u - 1.0 / c) / 2;

			eigMatrix[0][2] = -(B1 * v) / 2;
			eigMatrix[1][2] = -1;
			eigMatrix[2][2] = B1 * v;
			eigMatrix[3][2] = -B1 * v / 2;

			eigMatrix[0][3] = B1 / 2;
			eigMatrix[1][3] = 0;
			eigMatrix[2][3] = -B1;
			eigMatrix[3][3] = B1 / 2;
		}
		else // y-direction
		{
			eigMatrix[0][0] = (B2 + v / c) / 2;
			eigMatrix[1][0] = -u;
			eigMatrix[2][0] = 1 - B2;
			eigMatrix[3][0] = (B2 - v / c) / 2;

			eigMatrix[0][1] = -(B1 * u) / 2;
			eigMatrix[1][1] = 1;
			eigMatrix[2][1] = B1 * u;
			eigMatrix[3][1] = -(B1 * u) / 2;

			eigMatrix[0][2] = -(B1 * v + 1.0 / c) / 2;
			eigMatrix[1][2] = 0;
			eigMatrix[2][2] = B1 * v;
			eigMatrix[3][2] = -(B1 * v - 1.0 / c) / 2;

			eigMatrix[0][3] = B1 / 2;
			eigMatrix[1][3] = 0;
			eigMatrix[2][3] = -B1;
			eigMatrix[3][3] = B1 / 2;
		}
	}

	// (1-5) Right Eigen Matrix
	inline void getREigenMatrix(const Array1D<double> &Uh, double nx, double ny, Array2D<double> &eigMatrix)
	{
		double rho = Uh[0];
		double u = Uh[1] / Uh[0];
		double v = Uh[2] / Uh[0];
		double E = Uh[3];
		double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		double c = sqrt(euler_gamma * p / rho);
		double H = (E + p) / rho;

		if (fabs(nx) > fabs(ny)) // x-direction
		{
			eigMatrix[0][0] = 1;
			eigMatrix[1][0] = u - c;
			eigMatrix[2][0] = v;
			eigMatrix[3][0] = H - c * u;

			eigMatrix[0][1] = 0;
			eigMatrix[1][1] = 0;
			eigMatrix[2][1] = -1;
			eigMatrix[3][1] = -v;

			eigMatrix[0][2] = 1;
			eigMatrix[1][2] = u;
			eigMatrix[2][2] = v;
			eigMatrix[3][2] = (pow(u, 2) + pow(v, 2)) / 2;

			eigMatrix[0][3] = 1;
			eigMatrix[1][3] = u + c;
			eigMatrix[2][3] = v;
			eigMatrix[3][3] = H + c * u;
		}
		else // y-direction
		{
			eigMatrix[0][0] = 1;
			eigMatrix[1][0] = u;
			eigMatrix[2][0] = v - c;
			eigMatrix[3][0] = H - c * v;

			eigMatrix[0][1] = 0;
			eigMatrix[1][1] = 1;
			eigMatrix[2][1] = 0;
			eigMatrix[3][1] = u;

			eigMatrix[0][2] = 1;
			eigMatrix[1][2] = u;
			eigMatrix[2][2] = v;
			eigMatrix[3][2] = (pow(u, 2) + pow(v, 2)) / 2;

			eigMatrix[0][3] = 1;
			eigMatrix[1][3] = u;
			eigMatrix[2][3] = v + c;
			eigMatrix[3][3] = H + c * v;
		}
	}

	// (1-6) Num flux - (a) LLF
	inline void getLLFRiemannFlux(const Array1D<double> UL, const Array1D<double> UR, Array1D<double> &Flux, double nx, double ny)
	{
		// LLF approximate Riemann flux
		int varNum = getVarNum();
		Array1D<double> FUL(varNum), FUR(varNum);
		double ws_L(0), ws_R(0), ws(0);

		getPhyFlux(UL, FUL, nx, ny);
		getPhyFlux(UR, FUR, nx, ny);
		ws_L = getMaxEigenValue(UL, nx, ny);
		ws_R = getMaxEigenValue(UR, nx, ny);
		ws = std::max(ws_L, ws_R);

		// conservative variable
		for (int r = 0; r != varNum; ++r)
			Flux[r] = 0.5 * (FUL[r] + FUR[r] - ws * (UR[r] - UL[r]));
	}

	// -------- (2) 方程的初边值条件与参数设置 ----------- //
	// (2-1) 需要用到的物理变量
	std::function<double(double, double)> u0, v0, pre0, rho0;

	// (2-2) 根据算例进行方程设置
	inline void setEquationParameters(TESTCASETYPE type)
	{
		switch (type)
		{
		case SMOOTH:
			xL = 0.0;
			xR = 2.0;
			yL = 0.0;
			yR = 2.0;
			outputtime = 2.0;
			leftBoundaryCondition = PERIOD;
			rightBoundaryCondition = PERIOD;
			topBoundaryCondition = PERIOD;
			bottomBoundaryCondition = PERIOD;
			euler_gamma = 1.4;

			rho0 = [this](double xP, double yP)
			{
				return 1 + 0.2 * sin(M_PI * (xP + yP));
			};
			u0 = [this](double x, double y)
			{
				return 1;
			};
			v0 = [this](double x, double y)
			{
				return 1;
			};
			pre0 = [this](double x, double y)
			{
				return 1;
			};

			// 函数解析解
			this->u_exact_exist = true;
			this->theVarExact = [](double x, double y, double t)
			{ return 1 + 0.2 * sin(M_PI * (x + y - 2 * t)); };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;
		case Riemann1:
			xL = 0.0;
			xR = 2.0;
			yL = 0.0;
			yR = 2.0;
			outputtime = 0.52;
			leftBoundaryCondition = slip;
			rightBoundaryCondition = slip;
			topBoundaryCondition = slip;
			bottomBoundaryCondition = slip;
			euler_gamma = 1.4;

			rho0 = [this](double xP, double yP)
			{
				if (xP < 1.0 && yP < 1.0)
					return 0.8;
				else if (xP < 1.0 && yP >= 1.0)
					return 1.0;
				else if (xP >= 1.0 && yP < 1.0)
					return 1.0;
				else
					return 0.5313;

				std::cout << "initial value error!" << std::endl;
				return 0.0;
			};
			u0 = [this](double xP, double yP)
			{
				if (xP < 1.0 && yP < 1.0)
					return 0.0;
				else if (xP < 1.0 && yP >= 1.0)
					return 0.7276;
				else if (xP >= 1.0 && yP < 1.0)
					return 0.0;
				else
					return 0.0;

				std::cout << "initial value error!" << std::endl;
				return 0.0;
			};
			v0 = [this](double xP, double yP)
			{
				if (xP < 1.0 && yP < 1.0)
					return 0.0;
				else if (xP < 1.0 && yP >= 1.0)
					return 0.0;
				else if (xP >= 1.0 && yP < 1.0)
					return 0.7276;
				else
					return 0.0;

				std::cout << "initial value error!" << std::endl;
				return 0.0;
			};
			pre0 = [this](double xP, double yP)
			{
				if (xP < 1.0 && yP < 1.0)
					return 1.0;
				else if (xP < 1.0 && yP >= 1.0)
					return 1.0;
				else if (xP >= 1.0 && yP < 1.0)
					return 1.0;
				else
					return 0.4;

				std::cout << "initial value error!" << std::endl;
				return 0.0;
			};

			// 函数解析解
			this->u_exact_exist = false;
			this->theVarExact = [](double x, double y, double t)
			{ return 0.0; };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;
		case Riemann2:
			xL = 0.0;
			xR = 1.0;
			yL = 0.0;
			yR = 1.0;
			outputtime = 0.8;
			leftBoundaryCondition = symmetric;
			rightBoundaryCondition = symmetric;
			topBoundaryCondition = symmetric;
			bottomBoundaryCondition = symmetric;
			euler_gamma = 1.4;

			rho0 = [this](double xP, double yP)
			{
				if (xP < 0.8 && yP < 0.8)
					return 0.138;
				else if (xP < 0.8 && yP >= 0.8)
					return 0.5323;
				else if (xP >= 0.8 && yP < 0.8)
					return 0.5323;
				else
					return 1.5;

				std::cout << "initial value error!" << std::endl;
				return 0.0;
			};
			u0 = [this](double xP, double yP)
			{
				if (xP < 0.8 && yP < 0.8)
					return 1.206;
				else if (xP < 0.8 && yP >= 0.8)
					return 1.206;
				else if (xP >= 0.8 && yP < 0.8)
					return 0.0;
				else
					return 0.0;

				std::cout << "initial value error!" << std::endl;
				return 0.0;
			};
			v0 = [this](double xP, double yP)
			{
				if (xP < 0.8 && yP < 0.8)
					return 1.206;
				else if (xP < 0.8 && yP >= 0.8)
					return 0.0;
				else if (xP >= 0.8 && yP < 0.8)
					return 1.206;
				else
					return 0.0;

				std::cout << "initial value error!" << std::endl;
				return 0.0;
			};
			pre0 = [this](double xP, double yP)
			{
				if (xP < 0.8 && yP < 0.8)
					return 0.029;
				else if (xP < 0.8 && yP >= 0.8)
					return 0.3;
				else if (xP >= 0.8 && yP < 0.8)
					return 0.3;
				else
					return 1.5;

				std::cout << "initial value error!" << std::endl;
				return 0.0;
			};

			// 函数解析解
			this->u_exact_exist = false;
			this->theVarExact = [](double x, double y, double t)
			{ return 0.0; };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;
		case RTI:
			xL = 0.0;
			xR = 0.25;
			yL = 0.0;
			yR = 1.0;
			outputtime = 1.95;
			leftBoundaryCondition = slip;
			rightBoundaryCondition = slip;
			topBoundaryCondition = special;
			bottomBoundaryCondition = special;
			euler_gamma = 5.0 / 3.0;

			rho0 = [this](double x, double y)
			{
				if (y <= 0.5)
					return 2.0;
				else
					return 1.0;
			};
			u0 = [this](double x, double y)
			{
				return 0;
			};
			v0 = [this](double x, double y)
			{
				double pre = pre0(x, y);
				double rho = rho0(x, y);

				double c = sqrt(euler_gamma * pre / rho);

				return -0.025 * c * cos(8 * M_PI * x);
			};
			pre0 = [this](double x, double y)
			{
				if (y <= 0.5)
					return 1 + 2 * y;
				else
					return y + 1.5;
			};
			break;
		default:
			throw std::invalid_argument("Invalid example number");
		}
	}

	// (2-3) 将物理变量转换成双曲型方程 Ut+f(U)x = 0 的 U
	inline void getU0(double xP, double yP, Array1D<double> &U)
	{
		double rho = rho0(xP, yP);
		double u = u0(xP, yP);
		double v = v0(xP, yP);
		double pre = pre0(xP, yP);

		double E = pre / (euler_gamma - 1.0) + 0.5 * rho * (u * u + v * v);

		// Flux
		U[0] = rho;
		U[1] = rho * u;
		U[2] = rho * v;
		U[3] = E;
	}

	// ---------- (3) 输出 --------------//
	// (3-1) 输出变量的数量
	inline int getVitalVarNum()
	{
		return 6;
	}

	// (3-2) 输出变量的名字
	inline void getVitalVarName(Array1D<std::string> &VitalVarName)
	{
		VitalVarName[0] = "rho";
		VitalVarName[1] = "u";
		VitalVarName[2] = "v";
		VitalVarName[3] = "pre";
		VitalVarName[4] = "c";
		VitalVarName[5] = "E";
	}

	// (3-3) 输出变量的值
	inline void getVitalVarVal(const Array1D<double> Uh, Array1D<double> &VitalVar)
	{
		double rho = Uh[0];
		double u = Uh[1] / Uh[0];
		double v = Uh[2] / Uh[0];
		double E = Uh[3];
		double pre = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		double c = sqrt(euler_gamma * pre / rho);

		VitalVar[0] = rho;
		VitalVar[1] = u;
		VitalVar[2] = v;
		VitalVar[3] = pre;
		VitalVar[4] = c;
		VitalVar[5] = E;
	}
};

#endif