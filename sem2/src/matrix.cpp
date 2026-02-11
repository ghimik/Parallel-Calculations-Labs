#include "../include/matrix.hpp"
#include <random>
#include <cmath>

Matrix::Matrix(int size) : n(size), data(size * size, 0.0) {}

double& Matrix::operator()(int i, int j) {
    return data[i * n + j];
}

double Matrix::operator()(int i, int j) const {
    return data[i * n + j];
}

Matrix Matrix::identity(int n) {
    Matrix I(n);
    for (int i = 0; i < n; i++)
        I(i, i) = 1.0;
    return I;
}

Matrix Matrix::random(int n) {
    Matrix A(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dist(-10.0, 10.0);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A(i, j) = dist(gen);

    return A;
}

void Matrix::print() const {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            std::cout << (*this)(i, j) << " ";
        std::cout << "\n";
    }
}
