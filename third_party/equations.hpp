#ifndef __EQUATIONS__HPP__
#define __EQUATIONS__HPP__

enum RIEMANNFLUXTYPE
{
	LLF = 0
};
enum SCHEMETYPE
{
	WENOJS = 0,
	WENOZ = 1,
	WENOZP = 2,
	WENOZPI = 3
};
// 使用 inline 定义字符串数组
inline const std::string SCHEMETYPE_STRINGS[] = {
	"WENOJS", // SCHEMETYPE::WENO
	"WENOZ",  // SCHEMETYPE::WENOZ
	"WENOZP",
	"WENOZPI",
};

enum TESTCASETYPE
{
	SMOOTH = 0,	  // 光滑算例
	Riemann1 = 1, // 黎曼算例1
	Riemann2 = 2, // 黎曼算例2
	RTI = 3		  // 瑞丽泰勒不稳定算例
};
// 使用 inline 定义字符串数组
inline const std::string TESTCASETYPE_STRINGS[] = {
	"SMOOTH",
	"Riemann1",
	"Riemann2",
	"RTI",
};

enum BOUNDARYTYPE
{
	PERIOD = 0,	  // 周期边界条件
	SLIP = 1,	  // 无滑移边界条件
	SPECIAL = 2,  // 流入边界条件
	SYMMETRIC = 3 // 流出边界条件
};

// 方程

class EulerEquation
{
public:
	double xL, xR, yL, yR;
	double outputtime;
	BOUNDARYTYPE leftBoundaryCondition, rightBoundaryCondition;
	BOUNDARYTYPE topBoundaryCondition, bottomBoundaryCondition;

	std::function<double(const double, const double, const double)> theVarExact;
	std::function<double(const Array1D<double>)> theVarUh;
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
	inline void getPhyFlux(const Array1D<double> &Uh, Array1D<double> &Flux, double nx, double ny)
	{
		// U
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));

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

	inline void getPhyFluxX(const Array1D<double> &Uh, Array1D<double> &Flux)
	{
		// U
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));

		// Flux
		Flux[0] = rho * u;
		Flux[1] = rho * u * u + p;
		Flux[2] = rho * u * v;
		Flux[3] = (E + p) * u;
	}

	inline void getPhyFluxY(const Array1D<double> &Uh, Array1D<double> &Flux)
	{
		// U
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));

		Flux[0] = rho * v;
		Flux[1] = rho * u * v;
		Flux[2] = rho * v * v + p;
		Flux[3] = (E + p) * v;
	}

	// (1-3) eigen Values
	inline double getMaxEigenValue(const Array1D<double> &Uh, double nx, double ny)
	{
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		const double c = sqrt(euler_gamma * p / rho);

		return fabs(u * nx + v * ny) + c; // MaxeigValue
	}

	inline double getMaxEigenValueX(const Array1D<double> &Uh)
	{
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		const double c = sqrt(euler_gamma * p / rho);

		return fabs(u) + c; // MaxeigValue
	}

	inline double getMaxEigenValueY(const Array1D<double> &Uh)
	{
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		const double c = sqrt(euler_gamma * p / rho);

		return fabs(v) + c; // MaxeigValue
	}

	// (1-4) Left Eigen Matrix
	inline void getLEigenMatrixX(const Array1D<double> &Uh, Array2D<double> &eigMatrix)
	{
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		const double c = sqrt(euler_gamma * p / rho);
		const double B1 = (euler_gamma - 1) / (c * c);
		const double B2 = B1 * (u * u + v * v) / 2;

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

	inline void getLEigenMatrixY(const Array1D<double> &Uh, Array2D<double> &eigMatrix)
	{
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		const double c = sqrt(euler_gamma * p / rho);
		const double B1 = (euler_gamma - 1) / (c * c);
		const double B2 = B1 * (u * u + v * v) / 2;

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

	// (1-5) Right Eigen Matrix
	inline void getREigenMatrixX(const Array1D<double> &Uh, Array2D<double> &eigMatrix)
	{
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		const double c = sqrt(euler_gamma * p / rho);
		const double H = (E + p) / rho;

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
		eigMatrix[3][2] = (u * u + v * v) / 2;

		eigMatrix[0][3] = 1;
		eigMatrix[1][3] = u + c;
		eigMatrix[2][3] = v;
		eigMatrix[3][3] = H + c * u;
	}

	inline void getREigenMatrixY(const Array1D<double> &Uh, Array2D<double> &eigMatrix)
	{
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double p = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		const double c = sqrt(euler_gamma * p / rho);
		const double H = (E + p) / rho;

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
		eigMatrix[3][2] = (u * u + v * v) / 2;

		eigMatrix[0][3] = 1;
		eigMatrix[1][3] = u;
		eigMatrix[2][3] = v + c;
		eigMatrix[3][3] = H + c * v;
	}

	// (1-6) Num flux - (a) LLF
	inline void getLLFRiemannFlux(const Array1D<double> &UL, const Array1D<double> &UR, Array1D<double> &Flux, double nx, double ny)
	{
		// LLF approximate Riemann flux
		const int varNum = getVarNum();
		const double ws_L = getMaxEigenValue(UL, nx, ny);
		const double ws_R = getMaxEigenValue(UR, nx, ny);
		const double ws = std::max(ws_L, ws_R);

		Array1D<double> FUL(varNum), FUR(varNum);
		getPhyFlux(UL, FUL, nx, ny);
		getPhyFlux(UR, FUR, nx, ny);

		// conservative variable
		for (int r = 0; r != varNum; ++r)
		{
			Flux[0] = 0.5 * (FUL[0] + FUR[0] - ws * (UR[0] - UL[0]));
			Flux[1] = 0.5 * (FUL[1] + FUR[1] - ws * (UR[1] - UL[1]));
			Flux[2] = 0.5 * (FUL[2] + FUR[2] - ws * (UR[2] - UL[2]));
			Flux[3] = 0.5 * (FUL[3] + FUR[3] - ws * (UR[3] - UL[3]));
		}
	}

	inline void getLLFRiemannFluxX(const Array1D<double> &UL, const Array1D<double> &UR, Array1D<double> &Flux)
	{
		// LLF approximate Riemann flux
		const int varNum = getVarNum();
		const double ws_L = getMaxEigenValueX(UL);
		const double ws_R = getMaxEigenValueX(UR);
		const double ws = std::max(ws_L, ws_R);

		Array1D<double> FUL(varNum), FUR(varNum);
		getPhyFluxX(UL, FUL);
		getPhyFluxX(UR, FUR);

		// conservative variable
		for (int r = 0; r != varNum; ++r)
		{
			Flux[0] = 0.5 * (FUL[0] + FUR[0] - ws * (UR[0] - UL[0]));
			Flux[1] = 0.5 * (FUL[1] + FUR[1] - ws * (UR[1] - UL[1]));
			Flux[2] = 0.5 * (FUL[2] + FUR[2] - ws * (UR[2] - UL[2]));
			Flux[3] = 0.5 * (FUL[3] + FUR[3] - ws * (UR[3] - UL[3]));
		}
	}

	inline void getLLFRiemannFluxY(const Array1D<double> &UL, const Array1D<double> &UR, Array1D<double> &Flux)
	{
		// LLF approximate Riemann flux
		const int varNum = getVarNum();
		const double ws_L = getMaxEigenValueY(UL);
		const double ws_R = getMaxEigenValueY(UR);
		const double ws = std::max(ws_L, ws_R);

		Array1D<double> FUL(varNum), FUR(varNum);
		getPhyFluxY(UL, FUL);
		getPhyFluxY(UR, FUR);

		// conservative variable
		for (int r = 0; r != varNum; ++r)
		{
			Flux[0] = 0.5 * (FUL[0] + FUR[0] - ws * (UR[0] - UL[0]));
			Flux[1] = 0.5 * (FUL[1] + FUR[1] - ws * (UR[1] - UL[1]));
			Flux[2] = 0.5 * (FUL[2] + FUR[2] - ws * (UR[2] - UL[2]));
			Flux[3] = 0.5 * (FUL[3] + FUR[3] - ws * (UR[3] - UL[3]));
		}
	}

	inline void fluxSplit(const Array1D<double> &uh, double nx, double ny, Array1D<double> &flux_minus, Array1D<double> &flux_plus)
	{
		// Split the flux into positive and negative parts
		Array1D<double> F(4);
		getPhyFlux(uh, F, nx, ny);

		const double ws = getMaxEigenValue(uh, nx, ny);

		for (int r = 0; r != 4; ++r)
		{
			flux_minus[r] = 0.5 * (F[r] - ws * uh[r]);
			flux_plus[r] = 0.5 * (F[r] + ws * uh[r]);
		}
	}

	inline void fluxSplitX(const Array1D<double> &uh, Array1D<double> &flux_minus, Array1D<double> &flux_plus)
	{
		// Split the flux into positive and negative parts
		Array1D<double> F(4);
		getPhyFluxX(uh, F);

		const double ws = getMaxEigenValueX(uh);

		for (int r = 0; r != 4; ++r)
		{
			flux_minus[r] = 0.5 * (F[r] - ws * uh[r]);
			flux_plus[r] = 0.5 * (F[r] + ws * uh[r]);
		}
	}

	inline void fluxSplitY(const Array1D<double> &uh, Array1D<double> &flux_minus, Array1D<double> &flux_plus)
	{
		// Split the flux into positive and negative parts
		Array1D<double> F(4);
		getPhyFluxY(uh, F);

		const double ws = getMaxEigenValueY(uh);

		for (int r = 0; r != 4; ++r)
		{
			flux_minus[r] = 0.5 * (F[r] - ws * uh[r]);
			flux_plus[r] = 0.5 * (F[r] + ws * uh[r]);
		}
	}

	// -------- (2) 方程的初边值条件与参数设置 ----------- //
	// (2-1) 需要用到的物理变量
	std::function<double(const double, const double)> u0, v0, pre0, rho0;

	// (2-2) 根据算例进行方程设置
	inline void setEquationParameters(TESTCASETYPE type)
	{
		switch (type)
		{
		case SMOOTH:
			// Qiu, J., & Shu, C.-W. (2005).
			// Runge--Kutta Discontinuous Galerkin Method Using WENO Limiters.
			// SIAM Journal on Scientific Computing, 26(3), 907-929. doi:10.1137/s1064827503425298
			xL = 0.0;
			xR = 2.0;
			yL = 0.0;
			yR = 2.0;
			outputtime = 1.0;
			leftBoundaryCondition = PERIOD;
			rightBoundaryCondition = PERIOD;
			topBoundaryCondition = PERIOD;
			bottomBoundaryCondition = PERIOD;
			euler_gamma = 1.4;

			rho0 = [this](const double xP, const double yP)
			{
				return 1 + 0.2 * sin(M_PI * (xP + yP));
			};
			u0 = [this](const double x, const double y)
			{
				return 0.7;
			};
			v0 = [this](const double x, const double y)
			{
				return 0.3;
			};
			pre0 = [this](const double x, const double y)
			{
				return 1;
			};

			// 函数解析解
			this->u_exact_exist = true;
			this->theVarExact = [](const double x, const double y, const double t)
			{ return 1 + 0.2 * sin(M_PI * (x + y - t)); };
			this->theVarUh = [](Array1D<double> Uh)
			{ return Uh[0]; };
			break;
		case Riemann1:
			// Fan, C., & Wu, K. (2024).
			// High-order oscillation-eliminating Hermite WENO method for hyperbolic conservation laws.
			// Journal of Computational Physics, 519. doi:10.1016/j.jcp.2024.113435
			xL = 0.0;
			xR = 2.0;
			yL = 0.0;
			yR = 2.0;
			outputtime = 0.52;
			leftBoundaryCondition = SYMMETRIC;
			rightBoundaryCondition = SYMMETRIC;
			topBoundaryCondition = SYMMETRIC;
			bottomBoundaryCondition = SYMMETRIC;
			euler_gamma = 1.4;

			rho0 = [this](const double xP, const double yP)
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
			u0 = [this](const double xP, const double yP)
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
			v0 = [this](const double xP, const double yP)
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
			this->theVarExact = [](const double x, const double y, const double t)
			{ return 0.0; };
			this->theVarUh = [](const Array1D<double> &Uh)
			{ return Uh[0]; };
			break;
		case Riemann2:
			// Luo, X., &Wu, S.- p.(2021).
			// Improvement of the WENO - Z + scheme.
			// Computers &Fluids, 218. doi : 10.1016 / j.compfluid.2021.104855
			xL = 0.0;
			xR = 1.0;
			yL = 0.0;
			yR = 1.0;
			outputtime = 0.8;
			leftBoundaryCondition = SYMMETRIC;
			rightBoundaryCondition = SYMMETRIC;
			topBoundaryCondition = SYMMETRIC;
			bottomBoundaryCondition = SYMMETRIC;
			euler_gamma = 1.4;

			rho0 = [this](const double xP, const double yP)
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
			u0 = [this](const double xP, const double yP)
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
			v0 = [this](const double xP, const double yP)
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
			pre0 = [this](const double xP, const double yP)
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
			this->theVarExact = [](const double x, const double y, const double t)
			{ return 0.0; };
			this->theVarUh = [](const Array1D<double> &Uh)
			{ return Uh[0]; };
			break;
		case RTI:
			// Fleischmann, N., Adami, S., & Adams, N. A. (2019).
			// Numerical symmetry-preserving techniques for low-dissipation shock-capturing schemes.
			// Computers & Fluids, 189, 94-107. doi:10.1016/j.compfluid.2019.04.004
			xL = 0.0;
			xR = 0.25;
			yL = 0.0;
			yR = 1.0;
			outputtime = 1.95;
			leftBoundaryCondition = SLIP;
			rightBoundaryCondition = SLIP;
			topBoundaryCondition = SPECIAL;
			bottomBoundaryCondition = SPECIAL;
			euler_gamma = 5.0 / 3.0;

			rho0 = [this](const double x, const double y)
			{
				if (y <= 0.5)
					return 2.0;
				else
					return 1.0;
			};
			u0 = [this](const double x, const double y)
			{
				return 0;
			};
			v0 = [this](const double x, const double y)
			{
				const double pre = pre0(x, y);
				const double rho = rho0(x, y);
				const double c = sqrt(euler_gamma * pre / rho);

				return -0.025 * c * cos(8 * M_PI * x);
			};
			pre0 = [this](const double x, const double y)
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
	inline void getU0(const double xP, const double yP, Array1D<double> &U)
	{
		const double rho = rho0(xP, yP);
		const double u = u0(xP, yP);
		const double v = v0(xP, yP);
		const double pre = pre0(xP, yP);
		const double E = pre / (euler_gamma - 1.0) + 0.5 * rho * (u * u + v * v);

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
	inline void getVitalVarVal(const Array1D<double> &Uh, Array1D<double> &VitalVar)
	{
		const double rho = Uh[0];
		const double u = Uh[1] / Uh[0];
		const double v = Uh[2] / Uh[0];
		const double E = Uh[3];
		const double pre = (euler_gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
		const double c = sqrt(euler_gamma * pre / rho);

		VitalVar[0] = rho;
		VitalVar[1] = u;
		VitalVar[2] = v;
		VitalVar[3] = pre;
		VitalVar[4] = c;
		VitalVar[5] = E;
	}
};

#endif