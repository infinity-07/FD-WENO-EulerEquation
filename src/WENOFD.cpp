#include <mpi.h>

#include "../include/WENOFD.hpp"
#include "../include/Utility.hpp"

CWENOFD::~CWENOFD()
{
}

void CWENOFD::initializeSolver(std::map<std::string, std::string> option)
{
    // Note: Variables containing 'world' refer to the entire computational domain
    m_start_time = MPI_Wtime();

    // Get thread rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank); // 线程编号
    MPI_Comm_size(MPI_COMM_WORLD, &m_size); // 线程总数

    m_procIndexX = m_rank;

    if (m_rank == 0)
        std::cout << "Initializing variables..." << std::endl;

    m_testcase = (TESTCASETYPE)std::stoi(option["TESTCASE"]);
    m_rkMeth = (RKMETHOD)std::stoi(option["RKMETH"]);

    // Set equation parameters based on the testcase
    equation->setEquationParameters(m_testcase);
    m_outputTime = equation->outputtime;
    m_varNum = equation->getVarNum();

    // Initialize scheme, PPlimiter, and CFL number from options
    m_scheme = (SCHEMETYPE)std::stoi(option["SCHEME"]);
    m_cfl = std::stod(option["CFL"]);

    m_ghostCellNum = 3;

    // Initialize number of elements in X and Y directions
    m_worldPointNumX = std::stoi(option["ELEMNUMX"]);
    m_worldPointNumY = std::stoi(option["ELEMNUMY"]);

    m_worldStartPointX = m_ghostCellNum;
    m_worldStartPointY = m_ghostCellNum;
    m_worldEndPointX = m_worldStartPointX + m_worldPointNumX;
    m_worldEndPointY = m_worldStartPointY + m_worldPointNumY;

    const int basePointNumX = m_worldPointNumX / m_size; // 每个线程基本分配元素数
    const int remainder = m_worldPointNumX % m_size;     // 剩余元素数

    // 将余数分散到前多个线程
    m_numPointsX = basePointNumX + (m_rank < remainder ? 1 : 0);

    // 给线程rank分配localElemNumX个元素
    m_numPointsY = m_worldPointNumY;

    m_startPointX = m_ghostCellNum;
    m_startPointY = m_ghostCellNum;
    m_endPointX = m_startPointX + m_numPointsX;
    m_endPointY = m_startPointY + m_numPointsY;
    m_totalPointNumX = m_numPointsX + 2 * m_ghostCellNum;
    m_totalPointNumY = m_numPointsY + 2 * m_ghostCellNum;

    // Calculate element size in X and Y directions
    m_deltaX = (equation->xR - equation->xL) / (m_worldPointNumX - 1); // Structured grid
    m_deltaY = (equation->yR - equation->yL) / (m_worldPointNumY - 1); // Structured grid

    const double procXLeft = equation->xL + m_procIndexX * basePointNumX * m_deltaX;
    const double procYLeft = equation->yL;

    const double procGhostXLeft = procXLeft - m_deltaX * m_ghostCellNum;
    const double procGhostYLeft = procYLeft - m_deltaY * m_ghostCellNum;

    // Set the output directory based on the current timestamp and test case types
    m_outputDir = "./output/";
    m_outputDir += getTimestamp() + "_" + TESTCASETYPE_STRINGS[m_testcase] + "_";
    m_outputDir += SCHEMETYPE_STRINGS[m_scheme] + "_";
    m_outputDir += std::to_string(m_worldPointNumX) + "x" + std::to_string(m_worldPointNumY) + "/";

    if (m_rank == 0)
        backupProjectFiles();
    // if (m_rank == 0)
    //     std::cout << "Generating grid..." << std::endl;

    // Generate computational grid
    double loc_xCenter,
        loc_yCenter;
    m_grids.Resize(m_totalPointNumX, m_totalPointNumY);
    for (int ei = 0; ei < m_totalPointNumX; ++ei)
    {
        for (int ej = 0; ej < m_totalPointNumY; ++ej)
        {
            loc_xCenter = procGhostXLeft + ei * m_deltaX;
            loc_yCenter = procGhostYLeft + ej * m_deltaY;

            m_grids[ei][ej].m_xCenter = loc_xCenter;
            m_grids[ei][ej].m_yCenter = loc_yCenter;
        }
    }

    if (m_rank == 0)
    {
        m_worldGrids.Resize(m_worldPointNumX + 2 * m_ghostCellNum, m_worldPointNumY + 2 * m_ghostCellNum);
        for (int ei = 0; ei < m_worldPointNumX + 2 * m_ghostCellNum; ++ei)
        {
            for (int ej = 0; ej < m_worldPointNumY + 2 * m_ghostCellNum; ++ej)
            {
                loc_xCenter = equation->xL - m_ghostCellNum * m_deltaX + ei * m_deltaX;
                loc_yCenter = equation->yL - m_ghostCellNum * m_deltaY + ej * m_deltaY;

                m_worldGrids[ei][ej].m_xCenter = loc_xCenter;
                m_worldGrids[ei][ej].m_yCenter = loc_yCenter;
            }
        }
    }

    // Resize solution matrix and initialize variables
    m_Uh.Resize(m_totalPointNumX, m_totalPointNumY);
    for (int ei = 0; ei != m_totalPointNumX; ++ei)
        for (int ej = 0; ej != m_totalPointNumY; ++ej)
            m_Uh[ei][ej].vector.Resize(m_varNum);

    if (m_rank == 0)
    {
        // Initialize arrays to store global data
        m_worldUh.Resize(m_worldPointNumX + 2 * m_ghostCellNum, m_worldPointNumY + 2 * m_ghostCellNum);
        for (int ei = 0; ei != m_worldPointNumX + 2 * m_ghostCellNum; ei++)
            for (int ej = 0; ej != m_worldPointNumY + 2 * m_ghostCellNum; ej++)
                m_worldUh[ei][ej].vector.Resize(m_varNum);
    }

    // Resize and initialize Un matrix for storing solutions
    m_Un.Resize(m_totalPointNumX, m_totalPointNumY);
    for (int ei = 0; ei != m_totalPointNumX; ++ei)
        for (int ej = 0; ej != m_totalPointNumY; ++ej)
            m_Un[ei][ej].vector.Resize(m_varNum);

    // Resize and initialize RHS matrix
    m_rhs.Resize(m_totalPointNumX, m_totalPointNumY);
    for (int ei = 0; ei != m_totalPointNumX; ++ei)
        for (int ej = 0; ej != m_totalPointNumY; ++ej)
            m_rhs[ei][ej].vector.Resize(m_varNum);

    m_solverTimer.reset("solver");

    m_mainTimer.reset("main", &m_solverTimer);

    // m_wenoTimer.reset("weno", &m_mainTimer);
    // m_FDTimer.reset("FD", &m_mainTimer);
    // m_fluxSplitTimer.reset("flux split", &m_mainTimer);

    m_MPITimer.reset("mpi", &m_solverTimer);
    // m_outputTimer.reset("output", &m_solverTimer);
    m_initialiTimer.reset("initialization", &m_solverTimer);

    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFD::initializeAve(void)
{
    m_initialiTimer.start();
    Array1D<double> Conserved_var(m_varNum);

    if (m_rank == 0)
        std::cout << "Initializing cell averages..." << std::endl;

    for (int ei = m_startPointX; ei != m_endPointX; ++ei)
    {
        for (int ej = m_startPointY; ej != m_endPointY; ++ej)
        {
            m_Uh[ei][ej].vector.setZero();

            // Get the initial conserved variables at the Gauss point
            equation->getU0(m_grids[ei][ej].m_xCenter, m_grids[ei][ej].m_yCenter, Conserved_var);

            // Normalize by the cell area to get the average value
            for (int r = 0; r != m_varNum; ++r)
                m_Uh[ei][ej].vector[r] = Conserved_var[r];
        }
    }
    if (m_rank == 0)
        std::cout << "Cell averages initialization completed..." << std::endl;

    m_initialiTimer.pause();
}
void CWENOFD::run(std::map<std::string, std::string> option)
{
    initializeSolver(option);

    // Initialize current time
    m_now = 0;

    int count(0);
    double deltaT(0);
    double end_time;

    // Runge-Kutta iteration
    if (m_rank == 0)
        std::cout << "Start iteration..." << std::endl;

    m_solverTimer.start();

    initializeAve();

    outputAve("initial");

    m_mainTimer.start();

    while (fabs(m_now - m_outputTime) > 1e-9)
    {
        // Calculate time step
        deltaT = calculateDeltaT();

        // deltaT = 1e-3;

        // Perform RK3 method
        switch (m_rkMeth)
        {
        case RK1:
            RunRK1(deltaT);
            break;
        case RK2:
            RunRK2(deltaT);
            break;
        case RK3:
            RunRK3(deltaT);
            break;
        default:
            std::cout << "Error: Invalid Runge-Kutta method" << std::endl;
            std::cin.get();
            exit(1);
        }

        // Update current time and count iterations
        m_now += deltaT;
        count++;

        // Output progress every 100 iterations
        if (m_rank == 0)
        {
            end_time = MPI_Wtime();
            // showProgressBar(equation->outputtime, m_now, end_time - m_start_time);

            if (count % 20 == 1)
                std::cout << std::fixed << std::setprecision(2)
                          << "Iteration count: " << count
                          << std::scientific << std::setprecision(2)
                          << " DeltaT: " << deltaT
                          << std::fixed << std::setprecision(2)
                          << " Elapsed time: " << end_time - m_start_time
                          << " seconds" << std::endl;
        }

        // Output average values every 1000 iterations
        // if (count % 100 == 1)
        // {
        //     m_mainTimer.pause();
        //     outputAve("intermediate");
        //     m_mainTimer.start();
        // }
    }

    m_mainTimer.pause();

    outputAve("final");

    if (equation->u_exact_exist)
    {
        outputAccuracy("final");
        outputError("final");
    }

    if (m_rank == 0)
        m_solverTimer.printHierarchy();
}

void CWENOFD::RunRK1(double deltaT)
{
    // Assemble the right-hand side (RHS) for the equations
    assembleRHS();

    // Update solution vector using the RK2 coefficients
    for (int ei = m_startPointX; ei != m_endPointX; ei++)
        for (int ej = m_startPointY; ej != m_endPointY; ej++)
            for (int r = 0; r < m_varNum; r++)
                m_Uh[ei][ej].vector[r] += deltaT * m_rhs[ei][ej].vector[r];
}
void CWENOFD::RunRK2(double deltaT)
{
    // Copy current solution to temporary storage
    for (int ei = m_startPointX; ei != m_endPointX; ei++)
        for (int ej = m_startPointY; ej != m_endPointY; ej++)
            for (int r = 0; r != m_varNum; ++r)
                m_Un[ei][ej].vector[r] = m_Uh[ei][ej].vector[r];

    for (int step = 0; step != 2; step++)
    {
        // Assemble the right-hand side (RHS) for the equations
        assembleRHS();

        double a, b, c;
        // Set coefficients for the second-order Runge-Kutta scheme based on the step
        switch (step)
        {
        case 0:
            a = 1.0;
            b = 0.0;
            c = 1.0;
            break;
        case 1:
            a = 0.5;
            b = 0.5;
            c = 0.5;
            break;
        default:
            std::cout << "Error: Invalid RK2 step" << std::endl;
            std::cin.get();
            exit(1);
            break;
        }

        // Update solution vector using the RK2 coefficients
        for (int ei = m_startPointX; ei != m_endPointX; ei++)
            for (int ej = m_startPointY; ej != m_endPointY; ej++)
                for (int r = 0; r < m_varNum; r++)
                    m_Uh[ei][ej].vector[r] = a * m_Un[ei][ej].vector[r] + b * m_Uh[ei][ej].vector[r] + c * deltaT * m_rhs[ei][ej].vector[r];
    }
}
void CWENOFD::RunRK3(double deltaT)
{
    // Copy current solution to temporary storage
    for (int ei = m_startPointX; ei != m_endPointX; ei++)
        for (int ej = m_startPointY; ej != m_endPointY; ej++)
            for (int r = 0; r != m_varNum; ++r)
                m_Un[ei][ej].vector[r] = m_Uh[ei][ej].vector[r];

    for (int step = 0; step != 3; step++)
    {
        // Assemble the right-hand side (RHS) for the equations
        assembleRHS();

        double a, b, c;
        // Set coefficients for the third-order Runge-Kutta scheme based on the step
        switch (step)
        {
        case 0:
            a = 1.0;
            b = 0.0;
            c = 1.0;
            break;
        case 1:
            a = 3.0 / 4.0;
            b = 1.0 / 4.0;
            c = 1.0 / 4.0;
            break;
        case 2:
            a = 1.0 / 3.0;
            b = 2.0 / 3.0;
            c = 2.0 / 3.0;
            break;
        default:
            std::cout << "Error: Invalid RK3 step" << std::endl;
            std::cin.get();
            exit(1);
            break;
        }

        // Update solution vector using the RK3 coefficients
        for (int ei = m_startPointX; ei != m_endPointX; ei++)
            for (int ej = m_startPointY; ej != m_endPointY; ej++)
                for (int r = 0; r < m_varNum; r++)
                    m_Uh[ei][ej].vector[r] = a * m_Un[ei][ej].vector[r] + b * m_Uh[ei][ej].vector[r] + c * deltaT * m_rhs[ei][ej].vector[r];
    }
}

double CWENOFD::calculateDeltaT()
{
    // Initialize timestep to a large value
    double timestep(1.0);

    // Temporary variables for calculations
    double tmp(0);
    Array1D<double> loc_Uh(m_varNum);
    double eigenvalueX(0), eigenvalueY(0);

    // calculate the time step
    loc_Uh.setZero();
    for (int ei = m_startPointX; ei != m_endPointX; ei++)
    {
        for (int ej = m_startPointY; ej != m_endPointY; ej++)
        {
            // Extract average solution values for the cell
            for (int r = 0; r != m_varNum; ++r)
                loc_Uh[r] = m_Uh[ei][ej].vector[r];

            // Get maximum eigenvalues in X and Y directions
            eigenvalueX = equation->getMaxEigenValue(loc_Uh, 1, 0);
            eigenvalueY = equation->getMaxEigenValue(loc_Uh, 0, 1);

            // Calculate local timestep based on CFL condition
            tmp = m_cfl / (eigenvalueX / m_deltaX + eigenvalueY / m_deltaY);

            // Update timestep if smaller timestep is found
            if (tmp > timestep)
                continue;
            else
                timestep = tmp;
        }
    }

    // Adjust timestep if remaining time is smaller
    if (m_outputTime - m_now < timestep)
        timestep = m_outputTime - m_now;

    // Check for NaN timestep
    if (timestep != timestep)
    {
        std::cout << "Error: Timestep is not a number..." << std::endl;
        std::cin.get();
        exit(1);
    }

    // Obtain the minimum timestep across all processes using MPI_Reduce
    double localTimestep = timestep;
    double globalTimestep;
    MPI_Reduce(&localTimestep, &globalTimestep, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    // Broadcast the global minimum timestep to all processes
    MPI_Bcast(&globalTimestep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return globalTimestep;
}

void CWENOFD::assembleRHS(void)
{
    for (int ei = m_startPointX; ei != m_endPointX; ei++)
        for (int ej = m_startPointY; ej != m_endPointY; ej++)
            m_rhs[ei][ej].vector.setZero();

    // Calculate numerical fluxes using WENO scheme
    assembleBoundaryTerm();

    // Assemble source term for specific test cases
    assembleSourceTerm();
}
// void CWENOFD::fluxSplit(Array1D<double> uh, double nx, double ny, Array1D<double> &flux_minus, Array1D<double> &flux_plus)
// {
//     // Split the flux into positive and negative parts
//     Array1D<double> F(m_varNum);
//     equation->getPhyFlux(uh, F, nx, ny);

//     const double ws = equation->getMaxEigenValue(uh, nx, ny);

//     for (int r = 0; r != m_varNum; ++r)
//     {
//         flux_minus[r] = 0.5 * (F[r] - ws * uh[r]);
//         flux_plus[r] = 0.5 * (F[r] + ws * uh[r]);
//     }
// }

void CWENOFD::assembleBoundaryTerm(void)
{
    m_mainTimer.pause();
    exchangeGhostCellsValue();
    m_mainTimer.start();

    setBoundary();

    Array2D<double> fluxFaceX(m_totalPointNumX + 1, m_varNum);
    Array2D<double> flux_minusX(m_totalPointNumX, m_varNum);
    Array2D<double> flux_plusX(m_totalPointNumX, m_varNum);

    Array2D<double> fluxFaceY(m_totalPointNumY + 1, m_varNum);
    Array2D<double> flux_minusY(m_totalPointNumY, m_varNum);
    Array2D<double> flux_plusY(m_totalPointNumY, m_varNum);

    Array1D<double> loc_Uh(m_varNum);
    Array1D<double> loc_flux_minus(m_varNum);
    Array1D<double> loc_flux_plus(m_varNum);

    Array1D<double> aveFace(m_varNum);
    Array2D<double> eigenMatrixL(m_varNum, m_varNum);
    Array2D<double> eigenMatrixR(m_varNum, m_varNum);

    Array2D<double> v_plus(5, m_varNum);
    Array2D<double> v_minus(5, m_varNum);
    Array1D<double> v_half_plus(m_varNum);
    Array1D<double> v_half_minus(m_varNum);
    Array1D<double> v_half(m_varNum);
    Array1D<double> flux_half(m_varNum);

    // x-direction
    for (int ej = m_startPointY; ej != m_endPointY; ej++)
    {
        for (int ei = 0; ei != m_totalPointNumX; ei++)
        {
            equation->fluxSplitX(m_Uh[ei][ej].vector, loc_flux_minus, loc_flux_plus);
            Array1D<double> loc_flux_minus1(m_varNum), loc_flux_plus1(m_varNum);
            equation->fluxSplit(m_Uh[ei][ej].vector, 1, 0, loc_flux_minus1, loc_flux_plus1);

            double diff = 0;
            for (int r = 0; r != m_varNum; r++)
                diff += abs(loc_flux_minus1[r] - loc_flux_minus[r]) + abs(loc_flux_plus[r] - loc_flux_plus1[r]);
            if (diff > 1e-6)
            {
                std::cout << "WARNING: fluxSplitX 和 fluxSplit 结果不一致 (总差异: " << diff << ")\n";

                // 输出每个分量的详细差异
                std::cout << "Index\tfluxSplitX_minus\tfluxSplit_minus\t差值\t|\tfluxSplitX_plus\tfluxSplit_plus\t差值\n";
                for (int r = 0; r != m_varNum; r++)
                {
                    double minus_diff = std::abs(loc_flux_minus1[r] - loc_flux_minus[r]);
                    double plus_diff = std::abs(loc_flux_plus[r] - loc_flux_plus1[r]);

                    std::cout << r << "\t"
                              << loc_flux_minus[r] << "\t" << loc_flux_minus1[r] << "\t" << minus_diff << "\t|\t"
                              << loc_flux_plus[r] << "\t" << loc_flux_plus1[r] << "\t" << plus_diff << "\n";
                }

                // 可以在这里添加调试断点或异常处理
                // throw std::runtime_error("fluxSplit 结果不一致");
                std::cin.get();
            }

            for (int r = 0; r != m_varNum; ++r)
            {
                flux_minusX[ei][r] = loc_flux_minus[r];
                flux_plusX[ei][r] = loc_flux_plus[r];
            }
        }

        for (int faceI = m_startPointX; faceI != m_endPointX + 1; faceI++)
        {
            // 1. 计算界面平均量用于特征分解
            for (int r = 0; r != m_varNum; r++)
                aveFace[r] = 0.5 * (m_Uh[faceI - 1][ej].vector[r] + m_Uh[faceI][ej].vector[r]);

            equation->getLEigenMatrixX(aveFace, eigenMatrixL);
            equation->getREigenMatrixX(aveFace, eigenMatrixR);

            // 2. 将 flux_plus 和 flux_minus 投影到特征空间
            for (int s = -2; s <= 2; ++s) // 共5点模板
            {
                // 方法1
                // for (int i = 0; i < m_varNum; ++i) // 第 i 个特征变量
                // {
                //     v_plus[s + 2][i] = 0.0;
                //     v_minus[s + 2][i] = 0.0;
                //     for (int j = 0; j < m_varNum; ++j)
                //     {
                //         v_plus[s + 2][i] += eigenMatrixL[i][j] * flux_plusX[faceI + s - 1][j];
                //         v_minus[s + 2][i] += eigenMatrixL[i][j] * flux_minusX[faceI + s][j];
                //     }
                // }

                // 方法2
                // i = 0
                v_plus[s + 2][0] = eigenMatrixL[0][0] * flux_plusX[faceI + s - 1][0];
                v_plus[s + 2][0] += eigenMatrixL[0][1] * flux_plusX[faceI + s - 1][1];
                v_plus[s + 2][0] += eigenMatrixL[0][2] * flux_plusX[faceI + s - 1][2];
                v_plus[s + 2][0] += eigenMatrixL[0][3] * flux_plusX[faceI + s - 1][3];

                v_minus[s + 2][0] = eigenMatrixL[0][0] * flux_minusX[faceI + s][0];
                v_minus[s + 2][0] += eigenMatrixL[0][1] * flux_minusX[faceI + s][1];
                v_minus[s + 2][0] += eigenMatrixL[0][2] * flux_minusX[faceI + s][2];
                v_minus[s + 2][0] += eigenMatrixL[0][3] * flux_minusX[faceI + s][3];

                // i = 1
                v_plus[s + 2][1] = eigenMatrixL[1][0] * flux_plusX[faceI + s - 1][0];
                v_plus[s + 2][1] += eigenMatrixL[1][1] * flux_plusX[faceI + s - 1][1];
                v_plus[s + 2][1] += eigenMatrixL[1][2] * flux_plusX[faceI + s - 1][2];
                v_plus[s + 2][1] += eigenMatrixL[1][3] * flux_plusX[faceI + s - 1][3];

                v_minus[s + 2][1] = eigenMatrixL[1][0] * flux_minusX[faceI + s][0];
                v_minus[s + 2][1] += eigenMatrixL[1][1] * flux_minusX[faceI + s][1];
                v_minus[s + 2][1] += eigenMatrixL[1][2] * flux_minusX[faceI + s][2];
                v_minus[s + 2][1] += eigenMatrixL[1][3] * flux_minusX[faceI + s][3];

                // i = 2
                v_plus[s + 2][2] = eigenMatrixL[2][0] * flux_plusX[faceI + s - 1][0];
                v_plus[s + 2][2] += eigenMatrixL[2][1] * flux_plusX[faceI + s - 1][1];
                v_plus[s + 2][2] += eigenMatrixL[2][2] * flux_plusX[faceI + s - 1][2];
                v_plus[s + 2][2] += eigenMatrixL[2][3] * flux_plusX[faceI + s - 1][3];

                v_minus[s + 2][2] = eigenMatrixL[2][0] * flux_minusX[faceI + s][0];
                v_minus[s + 2][2] += eigenMatrixL[2][1] * flux_minusX[faceI + s][1];
                v_minus[s + 2][2] += eigenMatrixL[2][2] * flux_minusX[faceI + s][2];
                v_minus[s + 2][2] += eigenMatrixL[2][3] * flux_minusX[faceI + s][3];

                // i = 3
                v_plus[s + 2][3] = eigenMatrixL[3][0] * flux_plusX[faceI + s - 1][0];
                v_plus[s + 2][3] += eigenMatrixL[3][1] * flux_plusX[faceI + s - 1][1];
                v_plus[s + 2][3] += eigenMatrixL[3][2] * flux_plusX[faceI + s - 1][2];
                v_plus[s + 2][3] += eigenMatrixL[3][3] * flux_plusX[faceI + s - 1][3];

                v_minus[s + 2][3] = eigenMatrixL[3][0] * flux_minusX[faceI + s][0];
                v_minus[s + 2][3] += eigenMatrixL[3][1] * flux_minusX[faceI + s][1];
                v_minus[s + 2][3] += eigenMatrixL[3][2] * flux_minusX[faceI + s][2];
                v_minus[s + 2][3] += eigenMatrixL[3][3] * flux_minusX[faceI + s][3];
            }

            // 3. 对每个特征变量进行 WENO 重构
            for (int k = 0; k < m_varNum; ++k)
            {
                // 左偏重构：
                v_half_plus[k] = useWENO(
                    v_plus[0][k], v_plus[1][k], v_plus[2][k], v_plus[3][k], v_plus[4][k]);

                // 右偏重构：
                v_half_minus[k] = useWENO(
                    v_minus[4][k], v_minus[3][k], v_minus[2][k], v_minus[1][k], v_minus[0][k]);
            }

            // 4. 投影回物理空间：f_{i+1/2} = R * (v^+ + v^-)
            for (int k = 0; k < m_varNum; ++k)
                v_half[k] = v_half_plus[k] + v_half_minus[k];

            // 方法1
            // for (int i = 0; i < m_varNum; ++i)
            // {
            //     flux_half[i] = 0.0;
            //     for (int j = 0; j < m_varNum; ++j)
            //         flux_half[i] += eigenMatrixR[i][j] * v_half[j];
            // }
            // // 写入界面通量
            // for (int r = 0; r < m_varNum; ++r)
            //     fluxFaceX[faceI][r] = flux_half[r];

            // 方法2
            // 展开 flux_half 计算 (m_varNum = 4)
            fluxFaceX[faceI][0] = eigenMatrixR[0][0] * v_half[0];
            fluxFaceX[faceI][0] += eigenMatrixR[0][1] * v_half[1];
            fluxFaceX[faceI][0] += eigenMatrixR[0][2] * v_half[2];
            fluxFaceX[faceI][0] += eigenMatrixR[0][3] * v_half[3];

            fluxFaceX[faceI][1] = eigenMatrixR[1][0] * v_half[0];
            fluxFaceX[faceI][1] += eigenMatrixR[1][1] * v_half[1];
            fluxFaceX[faceI][1] += eigenMatrixR[1][2] * v_half[2];
            fluxFaceX[faceI][1] += eigenMatrixR[1][3] * v_half[3];

            fluxFaceX[faceI][2] = eigenMatrixR[2][0] * v_half[0];
            fluxFaceX[faceI][2] += eigenMatrixR[2][1] * v_half[1];
            fluxFaceX[faceI][2] += eigenMatrixR[2][2] * v_half[2];
            fluxFaceX[faceI][2] += eigenMatrixR[2][3] * v_half[3];

            fluxFaceX[faceI][3] = eigenMatrixR[3][0] * v_half[0];
            fluxFaceX[faceI][3] += eigenMatrixR[3][1] * v_half[1];
            fluxFaceX[faceI][3] += eigenMatrixR[3][2] * v_half[2];
            fluxFaceX[faceI][3] += eigenMatrixR[3][3] * v_half[3];
        }

        for (int ei = m_startPointX; ei != m_endPointX; ei++)
            for (int r = 0; r != m_varNum; ++r)
                m_rhs[ei][ej].vector[r] += -(fluxFaceX[ei + 1][r] - fluxFaceX[ei][r]) / m_deltaX;
    }

    // y-direction
    for (int ei = m_startPointX; ei != m_endPointX; ei++)
    {
        fluxFaceY.setZero();

        for (int ej = 0; ej != m_totalPointNumY; ej++)
        {
            // y 方向 unit normal: (0,1)
            equation->fluxSplit(m_Uh[ei][ej].vector, 0, 1, loc_flux_minus, loc_flux_plus);

            for (int r = 0; r != m_varNum; ++r)
            {
                flux_minusY[ej][r] = loc_flux_minus[r];
                flux_plusY[ej][r] = loc_flux_plus[r];
            }
        }

        for (int faceJ = m_startPointY; faceJ != m_endPointY + 1; ++faceJ)
        {
            // 1. 计算界面平均量用于特征分解（用 Uh 而不是 flux）
            for (int r = 0; r < m_varNum; ++r)
                aveFace[r] = 0.5 * (m_Uh[ei][faceJ - 1].vector[r] + m_Uh[ei][faceJ].vector[r]);

            // 获取特征矩阵 L 和 R（法向量 (0,1)）
            equation->getLEigenMatrixY(aveFace, eigenMatrixL);
            equation->getREigenMatrixY(aveFace, eigenMatrixR);

            // 2. 投影 flux_plus/flux_minus 到特征空间

            // 方法1
            // for (int s = -2; s <= 2; ++s)
            // {
            //     for (int i = 0; i < m_varNum; ++i)
            //     {
            //         v_plus[s + 2][i] = 0.0;
            //         v_minus[s + 2][i] = 0.0;
            //         for (int j = 0; j < m_varNum; ++j)
            //         {
            //             // 正向通量左偏：索引 faceJ+s-1
            //             v_plus[s + 2][i] += eigenMatrixL[i][j] * flux_plusY[faceJ + s - 1][j];

            //             // 负向通量右偏：索引 faceJ+s
            //             v_minus[s + 2][i] += eigenMatrixL[i][j] * flux_minusY[faceJ + s][j];
            //         }
            //     }
            // }

            // 方法2
            for (int s = -2; s <= 2; ++s)
            {
                // i = 0
                v_plus[s + 2][0] = eigenMatrixL[0][0] * flux_plusY[faceJ + s - 1][0];
                v_plus[s + 2][0] += eigenMatrixL[0][1] * flux_plusY[faceJ + s - 1][1];
                v_plus[s + 2][0] += eigenMatrixL[0][2] * flux_plusY[faceJ + s - 1][2];
                v_plus[s + 2][0] += eigenMatrixL[0][3] * flux_plusY[faceJ + s - 1][3];

                v_minus[s + 2][0] = eigenMatrixL[0][0] * flux_minusY[faceJ + s][0];
                v_minus[s + 2][0] += eigenMatrixL[0][1] * flux_minusY[faceJ + s][1];
                v_minus[s + 2][0] += eigenMatrixL[0][2] * flux_minusY[faceJ + s][2];
                v_minus[s + 2][0] += eigenMatrixL[0][3] * flux_minusY[faceJ + s][3];

                // i = 1
                v_plus[s + 2][1] = eigenMatrixL[1][0] * flux_plusY[faceJ + s - 1][0];
                v_plus[s + 2][1] += eigenMatrixL[1][1] * flux_plusY[faceJ + s - 1][1];
                v_plus[s + 2][1] += eigenMatrixL[1][2] * flux_plusY[faceJ + s - 1][2];
                v_plus[s + 2][1] += eigenMatrixL[1][3] * flux_plusY[faceJ + s - 1][3];

                v_minus[s + 2][1] = eigenMatrixL[1][0] * flux_minusY[faceJ + s][0];
                v_minus[s + 2][1] += eigenMatrixL[1][1] * flux_minusY[faceJ + s][1];
                v_minus[s + 2][1] += eigenMatrixL[1][2] * flux_minusY[faceJ + s][2];
                v_minus[s + 2][1] += eigenMatrixL[1][3] * flux_minusY[faceJ + s][3];

                // i = 2
                v_plus[s + 2][2] = eigenMatrixL[2][0] * flux_plusY[faceJ + s - 1][0];
                v_plus[s + 2][2] += eigenMatrixL[2][1] * flux_plusY[faceJ + s - 1][1];
                v_plus[s + 2][2] += eigenMatrixL[2][2] * flux_plusY[faceJ + s - 1][2];
                v_plus[s + 2][2] += eigenMatrixL[2][3] * flux_plusY[faceJ + s - 1][3];

                v_minus[s + 2][2] = eigenMatrixL[2][0] * flux_minusY[faceJ + s][0];
                v_minus[s + 2][2] += eigenMatrixL[2][1] * flux_minusY[faceJ + s][1];
                v_minus[s + 2][2] += eigenMatrixL[2][2] * flux_minusY[faceJ + s][2];
                v_minus[s + 2][2] += eigenMatrixL[2][3] * flux_minusY[faceJ + s][3];

                // i = 3
                v_plus[s + 2][3] = eigenMatrixL[3][0] * flux_plusY[faceJ + s - 1][0];
                v_plus[s + 2][3] += eigenMatrixL[3][1] * flux_plusY[faceJ + s - 1][1];
                v_plus[s + 2][3] += eigenMatrixL[3][2] * flux_plusY[faceJ + s - 1][2];
                v_plus[s + 2][3] += eigenMatrixL[3][3] * flux_plusY[faceJ + s - 1][3];

                v_minus[s + 2][3] = eigenMatrixL[3][0] * flux_minusY[faceJ + s][0];
                v_minus[s + 2][3] += eigenMatrixL[3][1] * flux_minusY[faceJ + s][1];
                v_minus[s + 2][3] += eigenMatrixL[3][2] * flux_minusY[faceJ + s][2];
                v_minus[s + 2][3] += eigenMatrixL[3][3] * flux_minusY[faceJ + s][3];
            }

            // 3. 对每个特征分量做 WENO 重构
            for (int k = 0; k < m_varNum; ++k)
            {
                // 左偏重构
                v_half_plus[k] = useWENO(v_plus[0][k], v_plus[1][k], v_plus[2][k], v_plus[3][k], v_plus[4][k]);

                // 右偏重构
                v_half_minus[k] = useWENO(v_minus[4][k], v_minus[3][k], v_minus[2][k], v_minus[1][k], v_minus[0][k]);
            }

            // 4. 投影回物理空间： f_{j+1/2} = R * (v^+ + v^-)
            for (int k = 0; k < m_varNum; ++k)
                v_half[k] = v_half_plus[k] + v_half_minus[k];

            // 方法1
            // for (int i = 0; i < m_varNum; ++i)
            // {
            //     flux_half[i] = 0.0;
            //     for (int j = 0; j < m_varNum; ++j)
            //         flux_half[i] += eigenMatrixR[i][j] * v_half[j];
            // }

            // // 写入 y 方向界面通量
            // for (int r = 0; r < m_varNum; ++r)
            //     fluxFaceY[faceJ][r] = flux_half[r];

            // 方法2：手动展开计算并直接写入 y 方向界面通量 (faceJ)
            fluxFaceY[faceJ][0] = eigenMatrixR[0][0] * v_half[0];
            fluxFaceY[faceJ][0] += eigenMatrixR[0][1] * v_half[1];
            fluxFaceY[faceJ][0] += eigenMatrixR[0][2] * v_half[2];
            fluxFaceY[faceJ][0] += eigenMatrixR[0][3] * v_half[3];

            fluxFaceY[faceJ][1] = eigenMatrixR[1][0] * v_half[0];
            fluxFaceY[faceJ][1] += eigenMatrixR[1][1] * v_half[1];
            fluxFaceY[faceJ][1] += eigenMatrixR[1][2] * v_half[2];
            fluxFaceY[faceJ][1] += eigenMatrixR[1][3] * v_half[3];

            fluxFaceY[faceJ][2] = eigenMatrixR[2][0] * v_half[0];
            fluxFaceY[faceJ][2] += eigenMatrixR[2][1] * v_half[1];
            fluxFaceY[faceJ][2] += eigenMatrixR[2][2] * v_half[2];
            fluxFaceY[faceJ][2] += eigenMatrixR[2][3] * v_half[3];

            fluxFaceY[faceJ][3] = eigenMatrixR[3][0] * v_half[0];
            fluxFaceY[faceJ][3] += eigenMatrixR[3][1] * v_half[1];
            fluxFaceY[faceJ][3] += eigenMatrixR[3][2] * v_half[2];
            fluxFaceY[faceJ][3] += eigenMatrixR[3][3] * v_half[3];
        }

        for (int ej = m_startPointY; ej != m_endPointY; ej++)
            for (int r = 0; r != m_varNum; ++r)
                m_rhs[ei][ej].vector[r] += -(fluxFaceY[ej + 1][r] - fluxFaceY[ej][r]) / m_deltaY;
    }
}

void CWENOFD::assembleSourceTerm(void)
{
    Array1D<double> Uh(m_varNum);
    switch (m_testcase)
    {
    case RTI:
    {
        const double g_g0 = 1.0;
        for (int ei = m_startPointX; ei != m_endPointX; ei++)
        {
            for (int ej = m_startPointY; ej != m_endPointY; ej++)
            {
                m_rhs[ei][ej].vector[3] += m_Uh[ei][ej].vector[2] * g_g0; // energy
                m_rhs[ei][ej].vector[2] += m_Uh[ei][ej].vector[0] * g_g0; // Moment y
            }
        }
        break;
    }
    default:
        break;
    }
}

double CWENOFD::useWENO(double &uavemm, double &uavem, double &uave, double &uavep, double &uavepp)
{
    // m_wenoTimer.start();
    double u_hat(0);
    switch (m_scheme)
    {
    case WENO:
        u_hat = WENO5threconstruction(uavemm, uavem, uave, uavep, uavepp);
        break;
    case WENOZ:
        u_hat = WENO5Zthreconstruction(uavemm, uavem, uave, uavep, uavepp);
        break;
    case WENOZP:
        u_hat = WENO5ZPtheconstruction(uavemm, uavem, uave, uavep, uavepp, m_deltaX);
        break;
    case WENOZPI:
        u_hat = WENO5ZPItheconstruction(uavemm, uavem, uave, uavep, uavepp);
        break;
    default:
        std::cout << "the scheme is not supported" << std::endl;
        std::cin.get();
        break;
    }
    // m_wenoTimer.pause();

    return u_hat;
}

void CWENOFD::exchangeGhostCellsValue(void)
{
    m_MPITimer.start();
    // 沿X方向与相邻进程交换ghost单元数据，用于后续WENO计算
    MPI_Barrier(MPI_COMM_WORLD);
    const int MLENGTH = m_varNum * m_ghostCellNum * m_numPointsY;
    std::vector<double> sendBufLeft(MLENGTH), sendBufRight(MLENGTH);
    std::vector<double> recvBufLeft(MLENGTH), recvBufRight(MLENGTH);

    // Fill sendBufLeft
    int index = 0;
    for (int i = 0; i < m_ghostCellNum; ++i)
    {
        int ei(0);
        if (m_rank == 0)
            ei = m_startPointX + i + 1;
        else
            ei = m_startPointX + i;

        for (int ej = m_startPointY; ej < m_endPointY; ++ej)
            for (int r = 0; r < m_varNum; ++r)
                sendBufLeft[index++] = m_Uh[ei][ej].vector[r];
    }

    // Fill sendBufRight
    index = 0;
    for (int i = 0; i < m_ghostCellNum; ++i)
    {
        int ei(0);
        if (m_rank == m_size - 1)
            ei = m_endPointX - m_ghostCellNum + i - 1;
        else
            ei = m_endPointX - m_ghostCellNum + i;

        for (int ej = m_startPointY; ej < m_endPointY; ++ej)
            for (int r = 0; r < m_varNum; ++r)
                sendBufRight[index++] = m_Uh[ei][ej].vector[r];
    }

    MPI_Request requests[4];
    MPI_Status statuses[4];

    // Define neighbors
    int leftRank = (m_rank == 0) ? m_size - 1 : m_rank - 1;
    int rightRank = (m_rank == m_size - 1) ? 0 : m_rank + 1;

    // Post receives
    MPI_Irecv(recvBufLeft.data(), MLENGTH, MPI_DOUBLE, leftRank, 1, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(recvBufRight.data(), MLENGTH, MPI_DOUBLE, rightRank, 2, MPI_COMM_WORLD, &requests[1]);

    // Send buffers
    MPI_Isend(sendBufRight.data(), MLENGTH, MPI_DOUBLE, rightRank, 1, MPI_COMM_WORLD, &requests[2]);
    MPI_Isend(sendBufLeft.data(), MLENGTH, MPI_DOUBLE, leftRank, 2, MPI_COMM_WORLD, &requests[3]);

    // Wait for completion
    MPI_Waitall(4, requests, statuses);

    // Unpack recvBufLeft into m_Uh[ei < m_startPointX]
    index = 0;
    for (int i = 0; i < m_ghostCellNum; ++i)
    {
        int ei = m_startPointX - m_ghostCellNum + i;
        for (int ej = m_startPointY; ej < m_endPointY; ++ej)
            for (int r = 0; r < m_varNum; ++r)
                m_Uh[ei][ej].vector[r] = recvBufLeft[index++];
    }

    // Unpack recvBufRight into m_Uh[ei >= m_endPointX]
    index = 0;
    for (int i = 0; i < m_ghostCellNum; ++i)
    {
        int ei = m_endPointX + i;
        for (int ej = m_startPointY; ej < m_endPointY; ++ej)
            for (int r = 0; r < m_varNum; ++r)
                m_Uh[ei][ej].vector[r] = recvBufRight[index++];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    m_MPITimer.pause();
}
void CWENOFD::GatherAllUhToRank0()
{
    // 收集所有线程的数据到线程0，用于输出
    MPI_Barrier(MPI_COMM_WORLD);

    // Step 1: 收集所有线程的 m_numPointsX 到线程0
    std::vector<int> m_numPointsX_per_rank(m_size, 0);
    MPI_Gather(&m_numPointsX, 1, MPI_INT, m_numPointsX_per_rank.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Step 2: 计算每个线程的 startPointX（线程0完成）
    std::vector<int> m_startPointX_per_rank(m_size, 0);
    if (m_rank == 0)
    {
        m_startPointX_per_rank[0] = 0;
        for (int i = 1; i < m_size; ++i)
            m_startPointX_per_rank[i] = m_startPointX_per_rank[i - 1] + m_numPointsX_per_rank[i - 1];
    }

    // Step 3: 广播 startPointX_per_rank 给所有线程（也可以只线程0用，但这样通用）
    MPI_Bcast(m_startPointX_per_rank.data(), m_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Step 4: 数据收集和发送
    if (m_rank == 0)
    {
        // 收集其他线程的数据
        for (int src_rank = 1; src_rank < m_size; ++src_rank)
        {
            int recvNX = m_numPointsX_per_rank[src_rank];
            std::vector<double> recvBuffer(recvNX * m_numPointsY * m_varNum);
            MPI_Status status;
            MPI_Recv(recvBuffer.data(), recvBuffer.size(), MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &status);

            int offsetX = m_startPointX_per_rank[src_rank];

            for (int ei = 0; ei < recvNX; ++ei)
            {
                for (int ej = 0; ej < m_numPointsY; ++ej)
                {
                    for (int r = 0; r < m_varNum; ++r)
                    {
                        int index = (ei * m_numPointsY * m_varNum) + (ej * m_varNum) + r;
                        m_worldUh[ei + offsetX + m_ghostCellNum][ej + m_ghostCellNum].vector[r] = recvBuffer[index];
                    }
                }
            }
        }

        // 拷贝自身数据
        int offsetX = m_startPointX_per_rank[0];
        for (int ei = 0; ei < m_numPointsX; ++ei)
        {
            for (int ej = 0; ej < m_numPointsY; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    m_worldUh[ei + offsetX + m_ghostCellNum][ej + m_ghostCellNum].vector[r] =
                        m_Uh[ei + m_ghostCellNum][ej + m_ghostCellNum].vector[r];
                }
            }
        }
    }
    else
    {
        // 序列化并发送 m_Uh 给线程0
        std::vector<double> sendBuffer(m_numPointsX * m_numPointsY * m_varNum);
        for (int ei = 0; ei < m_numPointsX; ++ei)
        {
            for (int ej = 0; ej < m_numPointsY; ++ej)
            {
                for (int r = 0; r < m_varNum; ++r)
                {
                    int index = (ei * m_numPointsY * m_varNum) + (ej * m_varNum) + r;
                    sendBuffer[index] = m_Uh[ei + m_ghostCellNum][ej + m_ghostCellNum].vector[r];
                }
            }
        }

        MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
void CWENOFD::setBoundary(void)
{
    // 给整个计算区域的 ghost cell 赋值 （使用边界条件）

    Array1D<double> loc_Uh(m_varNum);

    // y 方向
    for (int ei = 0; ei != m_totalPointNumX; ei++)
    {
        for (int e = 0; e != m_ghostCellNum; ++e)
        {
            switch (equation->bottomBoundaryCondition)
            {
            case PERIOD:
            {
                const int e1 = e + 1;
                for (int r = 0; r != m_varNum; ++r)
                {
                    m_Uh[ei][m_startPointY - e1].vector[r] = m_Uh[ei][m_endPointY - 1 - e1].vector[r]; // bottom
                    m_Uh[ei][m_endPointY - 1 + e1].vector[r] = m_Uh[ei][m_startPointY + e1].vector[r]; // top
                }
                break;
            }
            case SYMMETRIC:
            {
                const int e1 = e + 1;
                for (int r = 0; r != m_varNum; ++r)
                {
                    m_Uh[ei][m_startPointY - e1].vector[r] = m_Uh[ei][m_startPointY + e1].vector[r];     // bottom
                    m_Uh[ei][m_endPointY - 1 + e1].vector[r] = m_Uh[ei][m_endPointY - 1 - e1].vector[r]; // top
                }
                break;
            }
            case SLIP:
            {
                const int e1 = e + 1;
                for (int r = 0; r != m_varNum; ++r)
                {
                    m_Uh[ei][m_startPointY - e1].vector[r] = m_Uh[ei][m_startPointY + e1].vector[r];     // bottom
                    m_Uh[ei][m_endPointY - 1 + e1].vector[r] = m_Uh[ei][m_endPointY - 1 - e1].vector[r]; // top
                }

                m_Uh[ei][m_startPointY - e1].vector[2] *= -1;   // y direction
                m_Uh[ei][m_endPointY - 1 + e1].vector[2] *= -1; // y direction

                break;
            }
            case SPECIAL:
            {
                const int e1 = e + 1;

                equation->getU0(m_grids[ei][m_startPointY].m_xCenter, m_grids[ei][m_startPointY].m_yCenter, loc_Uh); // bottom
                for (int r = 0; r != m_varNum; ++r)
                    m_Uh[ei][m_startPointY - e1].vector[r] = loc_Uh[r]; // bottom

                equation->getU0(m_grids[ei][m_endPointY - 1].m_xCenter, m_grids[ei][m_endPointY - 1].m_yCenter, loc_Uh); // top
                for (int r = 0; r != m_varNum; ++r)
                    m_Uh[ei][m_endPointY - 1 + e1].vector[r] = loc_Uh[r]; // top

                break;
            }
            default:
                std::cout << "Error: Invalid boundary condition for y direction" << std::endl;
                break;
            }
        }
    }

    if (m_rank == 0)
    {
        // left
        for (int ej = 0; ej != m_totalPointNumY; ej++)
        {
            switch (equation->leftBoundaryCondition)
            {
            case SYMMETRIC:
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    const int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                        m_Uh[m_startPointX - e1][ej].vector[r] = m_Uh[m_startPointX + e1][ej].vector[r];
                }
                break;

            case SLIP:
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    const int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                        m_Uh[m_startPointX - e1][ej].vector[r] = m_Uh[m_startPointX + e1][ej].vector[r];

                    m_Uh[m_startPointX - e1][ej].vector[1] *= -1; // x direction | left
                }
                break;
            case SPECIAL:
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    const int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        equation->getU0(m_grids[m_startPointX][ej].m_xCenter, m_grids[m_startPointX][ej].m_yCenter, loc_Uh);
                        m_Uh[m_startPointX - e1][ej].vector[r] = loc_Uh[r]; // left
                    }
                }
                break;
            case PERIOD:
                //  在 exchangeGhostCellsValue 函数里已经顺手处理过了
                break;
            default:
                std::cout << "Error: Invalid boundary condition for x direction" << std::endl;
                break;
            }
        }
    }
    if (m_rank == m_size - 1) // 不能使用else if，因为可能 m_rank == 0 == m_size - 1
    {
        // right
        for (int ej = 0; ej != m_totalPointNumY; ej++)
        {
            switch (equation->rightBoundaryCondition)
            {
            case SYMMETRIC:
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    const int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                        m_Uh[m_endPointX - 1 + e1][ej].vector[r] = m_Uh[m_endPointX - 1 - e1][ej].vector[r];
                }
                break;
            case SLIP:
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    const int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                        m_Uh[m_endPointX - 1 + e1][ej].vector[r] = m_Uh[m_endPointX - 1 - e1][ej].vector[r];

                    m_Uh[m_endPointX - 1 + e1][ej].vector[1] *= -1; // x direction | right
                }
                break;

            case SPECIAL:
                for (int e = 0; e != m_ghostCellNum; ++e)
                {
                    const int e1 = e + 1;
                    for (int r = 0; r != m_varNum; ++r)
                    {
                        equation->getU0(m_grids[m_endPointX][ej].m_xCenter, m_grids[m_endPointX][ej].m_yCenter, loc_Uh);

                        m_Uh[m_endPointX - 1 + e1][ej].vector[r] = loc_Uh[r]; // left
                    }
                }
                break;
            case PERIOD:
                //  在 exchangeGhostCellsValue 函数里已经顺手处理过了
                break;
            default:
                std::cout << "Error: Invalid boundary condition for x direction" << std::endl;
                break;
            }
        }
    }
}

void CWENOFD::outputAve(std::string prefix)
{
    // 收集数据到线程0，使用线程0输出数据到 plt 文件里
    // m_MPITimer.start();
    exchangeGhostCellsValue();
    // m_MPITimer.pause();

    // m_outputTimer.start();
    setBoundary();
    GatherAllUhToRank0();

    // output
    if (m_rank == 0)
    {
        logMessage("outputing average solution...");

        int VitalVarNum = equation->getVitalVarNum();
        Array1D<double> VitalVar(VitalVarNum);
        Array1D<double> VitalVarAve(VitalVarNum);

        // Names of variables to be output
        Array1D<std::string> VitalVarName(VitalVarNum);
        equation->getVitalVarName(VitalVarName);

        // Output file name
        std::string filename = m_outputDir + "average_" + prefix + ".plt";
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << "TITLE=DGSolution" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "X ";
            outputFile << "Y ";
            for (int k = 0; k != VitalVarNum; k++)
                outputFile << VitalVarName[k] << " ";

            outputFile << std::endl;
            outputFile << "ZONE T=TA, ";
            outputFile << "I="; // Y - direction
            outputFile << m_worldPointNumY << ", ";
            outputFile << "J="; // X - direction
            outputFile << m_worldPointNumX << ", ";
            outputFile << "DATAPACKING = POINT" << std::endl;

            for (int ei = m_ghostCellNum; ei != m_worldPointNumX + m_ghostCellNum; ei++)
            {
                for (int ej = m_ghostCellNum; ej != m_worldPointNumY + m_ghostCellNum; ej++)
                {
                    Array1D<double> Uhh(m_varNum);
                    for (int r = 0; r != m_varNum; r++)
                        Uhh[r] = m_worldUh[ei][ej].vector[r];

                    // Output the values of the required variables one by one
                    equation->getVitalVarVal(Uhh, VitalVar);

                    outputFile << m_worldGrids[ei][ej].m_xCenter << " ";
                    outputFile << m_worldGrids[ei][ej].m_yCenter << " ";

                    for (int k = 0; k != VitalVarNum; k++)
                        outputFile << VitalVar[k] << " ";
                    outputFile << std::endl;
                }
            }

            outputFile.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
        logMessage("The file " + filename + " has been output successfully...");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // m_outputTimer.pause();
}
void CWENOFD::outputError(std::string prefix)
{
    // m_MPITimer.start();
    exchangeGhostCellsValue();
    // m_MPITimer.pause();

    GatherAllUhToRank0();

    // output
    if (m_rank == 0)
    {
        logMessage("outputing average solution...");

        int VitalVarNum = equation->getVitalVarNum();
        Array1D<double> VitalVar(VitalVarNum);
        Array1D<double> VitalVarAve(VitalVarNum);

        // Names of variables to be output
        Array1D<std::string> VitalVarName(VitalVarNum);
        equation->getVitalVarName(VitalVarName);

        // Output file name
        std::string filename = m_outputDir + "error_" + prefix + ".plt";
        std::ofstream outputFile(filename);

        if (outputFile.is_open())
        {
            outputFile << "TITLE=DGSolution" << std::endl;
            outputFile << "VARIABLES=";
            outputFile << "X ";
            outputFile << "Y ";
            outputFile << "error";
            outputFile << std::endl;
            outputFile << "ZONE T=TA, ";
            outputFile << "I="; // Y - direction
            outputFile << m_worldPointNumY << ", ";
            outputFile << "J="; // X - direction
            outputFile << m_worldPointNumX << ", ";
            outputFile << "DATAPACKING = POINT" << std::endl;

            for (int ei = m_ghostCellNum; ei != m_worldPointNumX + m_ghostCellNum; ei++)
            {
                for (int ej = m_ghostCellNum; ej != m_worldPointNumY + m_ghostCellNum; ej++)
                {
                    Array1D<double> Uhh(m_varNum);
                    for (int r = 0; r != m_varNum; r++)
                        Uhh[r] = m_worldUh[ei][ej].vector[r];

                    double Uref = equation->theVarExact(m_worldGrids[ei][ej].m_xCenter, m_worldGrids[ei][ej].m_yCenter, m_now);
                    double Ucal = equation->theVarUh(Uhh);

                    outputFile << m_worldGrids[ei][ej].m_xCenter << " ";
                    outputFile << m_worldGrids[ei][ej].m_yCenter << " ";
                    outputFile << fabs(Uref - Ucal) << " " << std::endl;
                }
            }

            outputFile.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
        logMessage("The file " + filename + " has been output successfully...");
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// void CWENOFD::copyConfig()
// {
//     // Input file path and name
//     std::filesystem::path inputFileName = "./input/config.cfg";
//     // Output file path and name
//     std::filesystem::path outputFileName = "./output/config_copy.cfg";
//     try
//     {
//         // Use std::filesystem::copy to copy the file
//         std::filesystem::copy(inputFileName, outputFileName, std::filesystem::copy_options::overwrite_existing);
//         std::cout << "File copied successfully." << std::endl;
//     }
//     catch (const std::filesystem::filesystem_error &e)
//     {
//         std::cerr << "Error: " << e.what() << std::endl;
//     }
// }
void CWENOFD::outputAccuracy(std::string prefix)
{
    if (m_rank == 0)
        std::cout << "Verifying accuracy..." << std::endl;

    // Arrays to store exact solution and error
    Array2D<double> Uexact(m_totalPointNumX + 2 * m_ghostCellNum, m_totalPointNumY + 2 * m_ghostCellNum);
    Array2D<double> Uerr(m_totalPointNumX + 2 * m_ghostCellNum, m_totalPointNumY + 2 * m_ghostCellNum);

    double Uref(0);
    double err_1(0), err_2(0), err_inf(0);

    // Compute exact solution averages over each cell
    Uexact.setZero();
    for (int ei = m_startPointX; ei != m_endPointX; ei++)
    {
        for (int ej = m_startPointY; ej != m_endPointY; ej++)
        {
            Uref = equation->theVarExact(m_grids[ei][ej].m_xCenter, m_grids[ei][ej].m_yCenter, m_now);
            Uexact[ei][ej] = Uref;
        }
    }

    // Compute error norms
    Array1D<double> Uave(m_varNum);
    for (int ei = m_startPointX; ei != m_endPointX; ei++)
    {
        for (int ej = m_startPointY; ej != m_endPointY; ej++)
        {
            for (int r = 0; r != m_varNum; ++r)
                Uave[r] = m_Uh[ei][ej].vector[r];

            // Calculate error at each cell
            Uerr[ei][ej] = fabs(Uexact[ei][ej] - equation->theVarUh(Uave));

            // Update error norms
            err_1 += Uerr[ei][ej];
            err_2 += pow(Uerr[ei][ej], 2);
            if (Uerr[ei][ej] > err_inf)
                err_inf = Uerr[ei][ej];
        }
    }

    double err1Sum, err2Sum, errinfMax;

    MPI_Allreduce(&err_1, &err1Sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&err_2, &err2Sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&err_inf, &errinfMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Calculate average error norms
    err1Sum = err1Sum / m_worldPointNumX / m_worldPointNumY;
    err2Sum = sqrt(err2Sum / m_worldPointNumX / m_worldPointNumY);

    if (m_rank == 0)
    {
        // Output to file
        std::fstream fileout;
        std::string filename = m_outputDir + "accuracy_" + prefix + ".csv";

        // Open file for writing
        std::ifstream fileExists(filename.c_str());
        if (fileExists)
            // File exists, append data
            fileout.open(filename.c_str(), std::ios::out | std::ios::app);
        else
        {
            // File doesn't exist, create a new file and write header
            fileout.open(filename.c_str(), std::ios::out);
            fileout << "Accuracy file does not exist, created new file..." << std::endl;
            fileout << "elemNum, Linf-norm (rho), L1-norm (rho), L2-norm (rho)" << std::endl;
        }

        fileout << m_worldPointNumX << ", ";
        fileout << std::setprecision(15) << std::setw(20) << std::setiosflags(std::ios::scientific) << errinfMax << ", ";
        fileout << std::setprecision(15) << std::setw(20) << std::setiosflags(std::ios::scientific) << err1Sum << ", ";
        fileout << std::setprecision(15) << std::setw(20) << std::setiosflags(std::ios::scientific) << err2Sum;

        fileout << std::endl;
        fileout.close();
        std::cout << "Accuracy verification completed..." << std::endl;
    }
}

inline void copyDirectory(const std::filesystem::path &src, const std::filesystem::path &dst)
{
    // copy the directory from src to dst
    try
    {
        std::filesystem::copy(src, dst,
                              std::filesystem::copy_options::recursive |
                                  std::filesystem::copy_options::overwrite_existing);
        // std::filesystem::copy_options::recursive 当复制源是目录时，递归复制其所有内容
        // std::filesystem::copy_options::overwrite_existing 目标文件已存在时直接覆盖
        std::cout << "复制成功: " << src << " -> " << dst << std::endl;
    }
    catch (const std::filesystem::filesystem_error &e)
    {
        std::cerr << "错误: " << e.what() << std::endl;
    }
}

void CWENOFD::backupProjectFiles()
{
    // 创建时间戳目录
    std::filesystem::path codeOutputDir = m_outputDir + "sourceCode/";
    std::filesystem::create_directories(codeOutputDir);
    std::filesystem::create_directories(codeOutputDir / "output/");

    // 需要复制的目录列表
    const std::vector<std::pair<std::filesystem::path, std::filesystem::path>> dirsToCopy = {
        {"./bin", codeOutputDir / "bin"},
        {"./build", codeOutputDir / "build"},
        {"./docs", codeOutputDir / "docs"},
        {"./include", codeOutputDir / "include"},
        {"./src", codeOutputDir / "src"},
        {"./makefile", codeOutputDir / "makefile"},
        {"./third_party", codeOutputDir / "third_party"}};

    // 执行复制操作
    for (const auto &[src, dst] : dirsToCopy)
    {
        if (std::filesystem::exists(src))
        {
            copyDirectory(src, dst);
        }
        else
        {
            std::cerr << "警告: 源目录不存在 " << src << std::endl;
        }
    }
}
