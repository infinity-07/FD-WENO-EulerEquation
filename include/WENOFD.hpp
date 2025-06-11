#ifndef __WENOFV__HPP__
#define __WENOFV__HPP__

#include "../third_party/scTools.h"

enum RKMETHOD
{
    RK1 = 0,
    RK2 = 1,
    RK3 = 2,
};

struct Cgrid
{
    double m_xCenter, m_yCenter;
};

struct Cvector
{
    Array1D<double> vector;
};

struct minusAndPlusFlux
{
    Array1D<double> flux_minus, flux_plus;
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
    int m_varNum;

public:
    int m_procIndexX;
    int m_rank, m_size;

    int m_worldPointNumX, m_worldPointNumY;
    int m_ghostCellNum;
    int m_startPointX, m_endPointX;
    int m_startPointY, m_endPointY;
    int m_numPointsX, m_numPointsY;
    int m_totalPointNumX, m_totalPointNumY;
    int m_worldStartPointX, m_worldStartPointY;
    int m_worldEndPointX, m_worldEndPointY;

    double m_globalXL, m_globalXR, m_globalYL, m_globalYR;
    double m_deltaX, m_deltaY;

    double m_outputTime;
    double m_now;
    double m_cfl;
    double m_start_time;

    TESTCASETYPE m_testcase;
    SCHEMETYPE m_scheme;
    RKMETHOD m_rkMeth;

    std::string m_outputDir;
    void backupProjectFiles();

public:
    Array2D<Cvector> m_rhs;
    Array2D<Cvector> m_Un;
    Array2D<Cvector> m_Uh;
    Array2D<Cgrid> m_grids;
    Array1D<minusAndPlusFlux> m_flux;

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

public:
    void run(std::map<std::string, std::string> m_option);

    void initializeSolver(std::map<std::string, std::string> option);
    void initializeAve();

    void assembleRHS();
    void assembleSourceTerm();
    void assembleBoundaryTerm();

    double calculateDeltaT();

    double useWENO(double &uavemm, double &uavem, double &uave, double &uavep, double &uavepp);

    // void fluxSplit(Array1D<double> uh, double nx, double ny, Array1D<double> &flux_minus, Array1D<double> &flux_plus);

    void outputAve(std::string prefix);
    void outputError(std::string prefix);
    void outputAccuracy(std::string prefix);

    void setBoundary(void);
    void exchangeGhostCellsValue(void);

public:
    void RunRK1(double deltaT);
    void RunRK2(double deltaT);
    void RunRK3(double deltaT);

public:
    void GatherAllUhToRank0();
    Array2D<Cgrid> m_worldGrids;
    Array2D<Cvector> m_worldUh;
};

///////////////////////////////////////////////////////////////////////
#endif