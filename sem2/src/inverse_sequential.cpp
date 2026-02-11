#include "../include/matrix.hpp"
#include <cmath>
#include <stdexcept>

Matrix inverse_sequential(const Matrix& A) {
    int n = A.n;

    Matrix aug(n);
    Matrix result = Matrix::identity(n);

    // Копируем A
    aug = A;

    for (int k = 0; k < n; k++) {

        // Partial pivoting
        int maxRow = k;
        double maxVal = std::abs(aug(k, k));

        for (int i = k + 1; i < n; i++) {
            if (std::abs(aug(i, k)) > maxVal) {
                maxVal = std::abs(aug(i, k));
                maxRow = i;
            }
        }

        if (maxVal < 1e-12)
            throw std::runtime_error("Matrix is singular");

        // Swap rows
        if (maxRow != k) {
            for (int j = 0; j < n; j++) {
                std::swap(aug(k, j), aug(maxRow, j));
                std::swap(result(k, j), result(maxRow, j));
            }
        }

        // Normalize pivot row
        double pivot = aug(k, k);
        for (int j = 0; j < n; j++) {
            aug(k, j) /= pivot;
            result(k, j) /= pivot;
        }

        // Eliminate other rows
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
