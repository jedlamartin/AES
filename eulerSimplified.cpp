#include "eulerSimplified.h"

const double EulerSimplified::Res = 5000.f;

const double EulerSimplified::Cap = 10e-6;

const double EulerSimplified::eps = 10e-6;

const double EulerSimplified::IS0 = 10e-9;

const double EulerSimplified::UT = 26e-3;

EulerSimplified::EulerSimplified(int N, double fs, bool increase) : N(N) {
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
        C(N - 1, N - 1) = 1 / EulerSimplified::Cap;
        this->D(0, 0) = -1;
        this->E(0, 0) = 1;
        for(int i = 0; i < N - 1; i++) {
            C(i, i) = 1 / EulerSimplified::Cap;
            C(i, i + 1) = -1 / EulerSimplified::Cap;
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
                A(i, 0) = -1.0 / (EulerSimplified::Cap * EulerSimplified::Res);
                A(i, 1) = 1.0 / (EulerSimplified::Cap * EulerSimplified::Res);
            } else if(i == N - 1) {
                // Last row
                A(i, N - 2) =
                    1.0 / (EulerSimplified::Cap * EulerSimplified::Res);
                A(i, N - 1) =
                    -1.0 / (EulerSimplified::Cap * EulerSimplified::Res);
            } else {
                // Intermediate rows
                A(i, i - 1) =
                    1.0 / (EulerSimplified::Cap * EulerSimplified::Res);
                A(i, i) =
                    -(1.0 / (EulerSimplified::Cap * EulerSimplified::Res) +
                      1.0 / (EulerSimplified::Cap * EulerSimplified::Res));
                A(i, i + 1) =
                    1.0 / (EulerSimplified::Cap * EulerSimplified::Res);
            }

            C(0, 0) = 1 / EulerSimplified::Cap;
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

void EulerSimplified::process(std::vector<double>& input,
                              std::vector<double>& output) {
                                  Matrix<double, Dynamic, 1> w;
                                  wBuffer.setZero(N, 1);
                                  int maxSize=N>1?1:0;
    Matrix<double, Dynamic, Dynamic> f;
    Matrix<double, Dynamic, Dynamic> J;
    Matrix<double, Dynamic, Dynamic> I;
    Matrix<double, Dynamic, 1> ones;
    I.setIdentity(maxSize+1, maxSize+1);
    ones.setOnes(maxSize+1, 1);

    
    const Matrix<double, Dynamic, Dynamic> Dtmp(this->D.block(0,0,maxSize,0));
    const Matrix<double, Dynamic, Dynamic> Gtmp(G.block(0,0,maxSize,G.cols()-1));
    const Matrix<double, Dynamic, Dynamic> Htmp(H.block(0,0,maxSize,H.cols()-1));
    Matrix<double, Dynamic, 1> wTmp;
    Matrix<double, Dynamic, 1> wMul;
    

    for(int i = 0; i < input.size(); i++) {
        w = wBuffer;
        w.block(maxSize+1,0,w.rows()-1,0)=this->H.block(maxSize+1,0,this->H.rows()-1,this->H.cols()-1)*wBuffer;
        while(true) {
            wTmp=w.block(0,0,maxSize,0);
            wMul.block(0,0,maxSize,0)=wBuffer.block(0,0,maxSize,0);
            wMul.block(maxSize+1,0,w.rows()-1,0)=w.block(maxSize+1,0,w.rows()-1,0)
            f = -wTmp + Htmp * wMul +
                Gtmp * IS0 *
                    (((Dtmp * wTmp + E * input[i]) / EulerSimplified::UT)
                         .unaryExpr([](double x) { return std::exp(x); }) -
                     ones);
            if((f.array().abs() < EulerSimplified::eps).all()) {
                break;
            }
            J = -I + Gtmp * EulerSimplified::IS0 / EulerSimplified::UT *
                         ((((Dtmp * wTmp + E * input[i]) / EulerSimplified::UT)
                               .unaryExpr([](double x) { return std::exp(x); }))
                              .asDiagonal()) *
                         Dtmp;
            
            w.block(0,0,maxSize,0) = w.block(0,0,maxSize,0) - J.inverse() * f;
        }
        output[i] = w(0, 0);
        wBuffer = w;
    }
}
