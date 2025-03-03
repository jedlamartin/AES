#ifndef KMETHOD_H
#define KMETHOD_H

#include <eigen3/Eigen/Dense>
#include <vector>

using namespace Eigen;

class EulerSimplified {
public:
    EulerSimplified(int, double, bool);
    void process(std::vector<double>& input, std::vector<double>& output);

private:
    int N;

    static const double Res;
    static const double Cap;
    static const double eps;
    static const double IS0;
    static const double UT;

    Matrix<double, Dynamic, 1> wBuffer;

    Matrix<double, Dynamic, Dynamic> D;
    Matrix<double, Dynamic, 1> E;
    Matrix<double, Dynamic, Dynamic> H;
    Matrix<double, Dynamic, Dynamic> G;
    Matrix<double, Dynamic, 1> J;
};
#endif