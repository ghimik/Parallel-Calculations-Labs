#pragma once
#include <vector>
#include <iostream>

class Matrix {
public:
    int n;
    std::vector<double> data;

    Matrix(int size);

    double& operator()(int i, int j);
    double operator()(int i, int j) const;

    static Matrix identity(int n);
    static Matrix random(int n);

    void print() const;
};
