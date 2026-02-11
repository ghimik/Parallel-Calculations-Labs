#include "../include/matrix.hpp"
#include "../include/timer.hpp"
#include <iostream>
#include <mpi.h>

Matrix inverse_sequential(const Matrix& A);

Matrix inverse_openmp(const Matrix& A);


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    int n = 512;

    Matrix A = Matrix::random(n);

    Timer timer;

    timer.start();
    Matrix inv1 = inverse_sequential(A);
    double t1 = timer.stop();

    timer.start();
    Matrix inv2 = inverse_openmp(A);
    double t2 = timer.stop();

    Timer timer;
    if(rank == 0) std::cout << "MPI inversion start...\n";

    timer.start();
    Matrix inv = inverse_mpi(A);
    double t3 = timer.stop();

    if(rank == 0) std::cout << "MPI inversion done: " << 3 << " sec\n";

    MPI_Finalize();

    timer.start();
    Matrix inv = inverse_cuda(A);
    double t4 = timer.stop();

    

    std::cout << "Sequential: " << t1 << " sec\n";
    std::cout << "OpenMP:     " << t2 << " sec\n";
    std::cout << "CUDA: " << t4 << " sec\n";


    return 0;
}
