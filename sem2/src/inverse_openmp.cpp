#include "../include/matrix.hpp"
#include <cmath>
#include <stdexcept>
#include <omp.h>

Matrix inverse_openmp(const Matrix& A) {
    int n = A.n;

    Matrix aug = A;
    Matrix result = Matrix::identity(n);

    for (int k = 0; k < n; k++) {

        // 🔹 Partial pivoting (поиск max строки)
        int maxRow = k;
        double maxVal = std::abs(aug(k, k));

        #pragma omp parallel
        {
            int localRow = maxRow;
            double localMax = maxVal;

            #pragma omp for nowait
            for (int i = k + 1; i < n; i++) {
                double val = std::abs(aug(i, k));
                if (val > localMax) {
                    localMax = val;
                    localRow = i;
                }
            }

            #pragma omp critical
            {
                if (localMax > maxVal) {
                    maxVal = localMax;
                    maxRow = localRow;
                }
            }
        }

        if (maxVal < 1e-12)
            throw std::runtime_error("Matrix is singular");

        // Swap rows (лучше не параллелить)
        if (maxRow != k) {
            for (int j = 0; j < n; j++) {
                std::swap(aug(k, j), aug(maxRow, j));
                std::swap(result(k, j), result(maxRow, j));
            }
        }

        // Нормализация pivot строки
        double pivot = aug(k, k);

        #pragma omp parallel for
        for (int j = 0; j < n; j++) {
            aug(k, j) /= pivot;
            result(k, j) /= pivot;
        }

        // Обнуление остальных строк
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {

            if (i == k) continue;

            double factor = aug(i, k);

            for (int j = 0; j < n; j++) {
                aug(i, j) -= factor * aug(k, j);
                result(i, j) -= factor * result(k, j);
            }
        }
    }

    return result;
}
