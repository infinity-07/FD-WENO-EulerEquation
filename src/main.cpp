#include "../third_party/scTools.h"
#include "../include/WENOFD.hpp"
#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // ------------------- Euler Equation  ----------------//
    EulerEquation *equation = new EulerEquation();

    std::string inputfile = "./input/config.cfg";

    // ---------------- Reading Configuration File - --------------//
    if (world_rank == 0)
        std::cout << "Starting to read configuration file..." << std::endl;

    CConfig config(inputfile);
    std::map<std::string, std::string> option;
    config.read(option);

    // ------------ Pass Equation and Initial Boundary Conditions to Solver -------------- //
    CWENOFD solver(equation);

    solver.run(option);

    MPI_Finalize();

    return 0;
}
