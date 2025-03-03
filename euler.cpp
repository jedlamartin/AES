#include "euler.h"

#include <iostream>

const double Euler::Res = 5000.f;

const double Euler::Cap = 10e-6;

const double Euler::eps = 10e-6;

const double Euler::IS0 = 10e-9;

const double Euler::UT = 26e-3;

Euler::Euler(int N, double fs, bool increase) : N(N) {
    Matrix<double, Dynamic, Dynamic> A;
    Matrix<double, Dynamic, 1> B;
    Matrix<double, Dynamic, Dynamic> C;
    Matrix<double, Dynamic, Dynamic> F;
    if(increase) {
        A.setZero(N, N);
        B.setZero(N, 1);
        C.setZero(N, N);
        this->D.setZero(N, N);
        this->E.setZero(N, 1);
        F.setZero(N, N);
        C(N - 1, N - 1) = 1 / Euler::Cap;
        this->D(0, 0) = -1;
        this->E(0, 0) = 1;
        for(int i = 0; i < N - 1; i++) {
            C(i, i) = 1 / Euler::Cap;
            C(i, i + 1) = -1 / Euler::Cap;
            this->D(i + 1, i) = 1;
            this->D(i + 1, i + 1) = -1;
        }
    } else {
        A.setZero(N, N);
        B.setZero(N, 1);
        C.setZero(N, 1);
        this->D.setZero(1, N);
        this->E.setZero(1, 1);
        F.setZero(1, 1);

        for(int i = 0; i < N; i++) {
            if(i == 0) {
                // First row
                A(i, 0) = -1.0 / (Euler::Cap * Euler::Res);
                A(i, 1) = 1.0 / (Euler::Cap * Euler::Res);
            } else if(i == N - 1) {
                // Last row
                A(i, N - 2) = 1.0 / (Euler::Cap * Euler::Res);
                A(i, N - 1) = -1.0 / (Euler::Cap * Euler::Res);
            } else {
                // Intermediate rows
                A(i, i - 1) = 1.0 / (Euler::Cap * Euler::Res);
                A(i, i) = -(1.0 / (Euler::Cap * Euler::Res) +
                            1.0 / (Euler::Cap * Euler::Res));
                A(i, i + 1) = 1.0 / (Euler::Cap * Euler::Res);
            }

            C(0, 0) = 1 / Euler::Cap;
            this->D(0, 0) = -1;
            this->E(0, 0) = 1;
        }
    }

    Matrix<double, Dynamic, Dynamic> I;
    I.setIdentity(N, N);
    this->H = fs * (fs * I - A).inverse();
    this->G = (fs * I - A).inverse() * C;
    this->J = (fs * I - A).inverse() * B;
}

void Euler::process(std::vector<double>& input, std::vector<double>& output) {
    Matrix<double, Dynamic, 1> w;
    wBuffer.setZero(N, 1);
    Matrix<double, Dynamic, Dynamic> f;
    Matrix<double, Dynamic, Dynamic> J;
    Matrix<double, Dynamic, Dynamic> I;
    Matrix<double, Dynamic, 1> ones;
    I.setIdentity(N, N);
    ones.setOnes(this->E.rows(), 1);
    for(int i = 0; i < input.size(); i++) {
        w = wBuffer;
        while(true) {
            f = -w + H * wBuffer +
                G * IS0 *
                    (((this->D * w + E * input[i]) / Euler::UT)
                         .unaryExpr([](double x) { return std::exp(x); }) -
                     ones);
            if((f.array().abs() < Euler::eps).all()) {
                break;
            }
            J = -I + G * Euler::IS0 / Euler::UT *
                         ((((D * w + E * input[i]) / Euler::UT)
                               .unaryExpr([](double x) { return std::exp(x); }))
                              .asDiagonal()) *
                         D;
            w = w - J.inverse() * f;
        }
        output[i] = w(0, 0);
        wBuffer = w;
    }
}
