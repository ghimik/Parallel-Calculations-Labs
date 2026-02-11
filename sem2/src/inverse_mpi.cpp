#include "../include/matrix.hpp"
#include "../include/timer.hpp"
#include <mpi.h>
#include <cmath>
#include <stdexcept>
#include <iostream>

Matrix inverse_mpi(const Matrix& A) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = A.n;
    Matrix aug = A;
    Matrix result = Matrix::identity(n);

    for (int k = 0; k < n; k++) {

        double localMax = -1.0;
        int localRow = -1;

        for (int i = k + rank; i < n; i += size) {
            double val = std::abs(aug(i, k));
            if (val > localMax) {
                localMax = val;
                localRow = i;
            }
        }

        struct { double value; int index; } localData{localMax, localRow}, globalData;
        MPI_Allreduce(&localData, &globalData, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        int maxRow = globalData.index;
        double maxVal = globalData.value;

        if (maxVal < 1e-12) {
            if (rank == 0) throw std::runtime_error("Matrix is singular");
            else return Matrix(n);
        }

        if (maxRow != k) {
            for (int j = 0; j < n; j++) {
                std::swap(aug(k, j), aug(maxRow, j));
                std::swap(result(k, j), result(maxRow, j));
            }
        }

        double pivot = aug(k, k);
        for (int j = 0; j < n; j++) {
            aug(k, j) /= pivot;
            result(k, j) /= pivot;
        }

        MPI_Bcast(&aug(k, 0), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&result(k, 0), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int i = k + rank; i < n; i += size) {
            if (i == k) continue;
            double factor = aug(i, k);
            for (int j = 0; j < n; j++) {
                aug(i, j) -= factor * aug(k, j);
                result(i, j) -= factor * result(k, j);
            }
        }

        for (int i = 0; i < n; i += size) {
            int owner = i % size;
            MPI_Bcast(&aug(i, 0), n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
            MPI_Bcast(&result(i, 0), n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
        }
    }

    return result;
}