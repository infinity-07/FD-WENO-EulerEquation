#ifndef __WENOFV__HPP__
#define __WENOFV__HPP__

#include <iostream>
#include <fstream>
#include "scTools/scTools.h"
#include "Equations.hpp"
#include <map>

class Cgrid
{
public:
    // 成员变量
    double m_xCenter;
    double m_yCenter;

    // 默认构造函数
    Cgrid()
        : m_xCenter(0), m_yCenter(0)
    {
    }

    // 带参数的构造函数
    Cgrid(double xCenter, double yCenter)
        : m_xCenter(xCenter), m_yCenter(yCenter)
    {
    }
};

class CgridFlux
{
public:
    // 成员变量
    Array1D<double> rightFlux;
    Array1D<double> leftFlux;
    Array1D<double> topFlux;
    Array1D<double> bottomFlux;
};

class Cvector
{
public:
    // 成员变量
    Array1D<double> vector;
};

class minusAndPlusFlux
{
public:
    Array1D<double> flux_minus;
    Array1D<double> flux_plus;
};

class CWENOFD
{
public:
    CWENOFD() {}
    ~CWENOFD();

    CWENOFD(EulerEquation *equation)
    {
        this->equation = equation;
    }

    EulerEquation *equation;

public:
    int m_procIndexX;
    int m_rank, m_size;
    int m_worldPointNumX, m_worldPointNumY;

    int m_ghostCellNum;
    int m_startPointX, m_endPointX;
    int m_startPointY, m_endPointY;
    Array2D<double> m_worldTotalTheta;

    int m_numPointsX, m_numPointsY;

    int m_totalPointNumX, m_totalPointNumY;

    int m_worldStartPointX, m_worldStartPointY;
    int m_worldEndPointX, m_worldEndPointY;
    int m_varNum;

    double m_globalXL;
    double m_globalXR;
    double m_globalYL;
    double m_globalYR;

    double m_cfl;
    double m_deltaX;
    double m_deltaY;
    double m_outputTime;
    double m_restartTime;
    double m_now;

    Array2D<double> m_totalTheta;

public:
    Array2D<Cvector> m_Uh;
    Array2D<Cvector> m_worldUh;
    Array2D<Cgrid> m_grids;
    Array2D<Cgrid> m_worldGrids;
    Array1D<minusAndPlusFlux> m_flux;
    Array2D<CgridFlux> m_gridFlux;
    Array2D<CgridFlux> m_worldGridFlux;
    double m_start_time;

public:
    Timer m_timer; // Timer
    Timer m_solverTimer;
    Timer m_mainTimer;
    Timer m_wenoTimer;
    Timer m_FDTimer;
    Timer m_outputTimer;
    Timer m_initialiTimer;
    Timer m_fluxSplitTimer;
    Timer m_MPITimer;

    void initializeSolver(std::map<std::string, std::string> option);
    void initializeAve();
    void outputAve(std::string filename);
    void outputError(std::string filename);
    void outputAccuracy(std::string prefix);
    double calculateDeltaT();

    void fluxSplit(Array1D<double> uh, double nx, double ny, Array1D<double> &flux_minus, Array1D<double> &flux_plus);

    void getFlux();
    void run(std::map<std::string, std::string> m_option);
    void assembleRHS();
    void assembleSourceTerm();
    // void getWorldUh();
    // void worldToProcUh();

    TESTCASETYPE m_testcase;
    SCHEMETYPE m_scheme;

    enum RKMETHOD
    {
        RK1 = 0,
        RK2 = 1,
        RK3 = 2,
    };

    RKMETHOD m_rkMeth;

    // void setBoundaryGPs(void);
    void setBoundary(void);
    void MPICommunication(void);

    Array2D<Cvector> m_rhs;
    Array2D<Cvector> m_Un;
    void RunRK1(double deltaT);
    void RunRK2(double deltaT);
    void RunRK3(double deltaT);

    void GatherAllUhToRank0();

    double useWENO(double uavemm, double uavem, double uave, double uavep, double uavepp, int gp);

    Array2D<double> m_worldTotalPr, m_totalPr;
    // void copyConfig();
    std::map<std::string, std::string> m_option;
};

///////////////////////////////////////////////////////////////////////
#endif