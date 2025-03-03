#include "K.h"

const double Kmethod::Res = 5000.f;

const double Kmethod::Cap = 10e-6;

const double Kmethod::IS0 = 10e-9;

const double Kmethod::UT = 26e-3;

Kmethod::Kmethod(int N, double fs, double eps, bool increase) : N(N) {
    int numOfTables = increase ? N : 1;
    this->table.resize(
        numOfTables,
        std::make_pair(std::vector<double>(), std::vector<double>()));

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
        C(N - 1, N - 1) = 1 / Kmethod::Cap;
        this->D(0, 0) = -1;
        this->E(0, 0) = 1;
        for(int i = 0; i < N - 1; i++) {
            C(i, i) = 1 / Kmethod::Cap;
            C(i, i + 1) = -1 / Kmethod::Cap;
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
                A(i, 0) = -1.0 / (Kmethod::Cap * Kmethod::Res);
                A(i, 1) = 1.0 / (Kmethod::Cap * Kmethod::Res);
            } else if(i == N - 1) {
                // Last row
                A(i, N - 2) = 1.0 / (Kmethod::Cap * Kmethod::Res);
                A(i, N - 1) = -1.0 / (Kmethod::Cap * Kmethod::Res);
            } else {
                // Intermediate rows
                A(i, i - 1) = 1.0 / (Kmethod::Cap * Kmethod::Res);
                A(i, i) = -(1.0 / (Kmethod::Cap * Kmethod::Res) +
                            1.0 / (Kmethod::Cap * Kmethod::Res));
                A(i, i + 1) = 1.0 / (Kmethod::Cap * Kmethod::Res);
            }

            C(0, 0) = 1 / Kmethod::Cap;
            this->D(0, 0) = -1;
            this->E(0, 0) = 1;
        }
    }
    double h = 2 * fs;
    Matrix<double, Dynamic, Dynamic> I;
    I.setIdentity(N, N);
    this->H = (h * I - A).inverse() * (h * I + A);
    this->G = (h * I - A).inverse() * C;
    this->J = (h * I - A).inverse() * B;
    this->K = D * this->G + F;

    // generate lookup tables
    std::pair<std::vector<double>, std::vector<double>> characteristic;
    for(double u = -1.f; u <= 1.f; u += eps) {
        double i = Kmethod::IS0 * (std::exp(u / Kmethod::UT) - 1);
        characteristic.second.push_back(i);
        characteristic.first.push_back(u);
    }
    for(int i = 0; i < numOfTables; i++) {
        int size = characteristic.first.size();
        std::vector<double> p(size, 0.f);
        for(int j = 0; j < size; j++) {
            double tmp = 0;
            for(int k = 0; k < numOfTables; k++) {
                tmp += K(i, k) * characteristic.second[j];
            }
            p[j] = characteristic.first[j] - tmp;
        }
        this->table[i] = std::make_pair(p, characteristic.second);
    }
}

void Kmethod::process(std::vector<double>& input, std::vector<double>& output) {
    Matrix<double, Dynamic, Dynamic> w;
    w.setZero(N, 1);

    Matrix<double, Dynamic, Dynamic> y;
    y.setZero(this->E.rows(), 1);

    Matrix<double, Dynamic, Dynamic> pk;
    Matrix<double, Dynamic, 1> p;

    double uBuffer = 0.f;

    for(int i = 0; i < input.size(); i++) {
        pk = H * w + J * (input[i] + uBuffer) + G * y;
        p = D * pk + E * input[i];
        y = binarySearch(p);
        w = G * y + pk;
        output[i] = w(0, 0);
        uBuffer = input[i];
    }
}

Matrix<double, Dynamic, 1> Kmethod::binarySearch(
    Matrix<double, Dynamic, 1>& p) const {
    Matrix<double, Dynamic, 1> y;
    y.setZero(p.rows(), 1);

    for(int index = 0; index < p.rows(); index++) {
        int size = this->table[index].first.size();
        // std::cout << p(index, 0) << std::endl;
        if(p(index, 0) <= this->table[index].first[0]) {
            y(index, 0) = this->table[index].second[0];
        } else if(p(index, 0) >= this->table[index].first[size - 1]) {
            y(index, 0) = this->table[index].second[size - 1];
        } else {
            int i = 0, j = this->table[index].first.size(), mid;
            while(i < j) {
                mid = (i + j) / 2;
                if(this->table[index].first[mid] == p(index, 0)) {
                    y(index, 0) = table[index].second[mid];
                    break;
                }

                if(p(index, 0) < this->table[index].first[mid]) {
                    if(mid > 0 &&
                       p(index, 0) > this->table[index].first[mid - 1]) {
                        mid = (p(index, 0) - this->table[index].first[mid - 1] <
                               this->table[index].first[mid] - p(index, 0)) ?
                                  mid - 1 :
                                  mid;
                        break;
                    }
                    j = mid;
                } else if(p(index, 0) > this->table[index].first[mid]) {
                    if(mid < size - 1 &&
                       p(index, 0) < this->table[index].first[mid + 1]) {
                        mid =
                            (p(index, 0) - this->table[index].first[mid] <
                             this->table[index].first[mid + 1] - p(index, 0)) ?
                                mid :
                                mid + 1;
                        break;
                    }
                    i = mid + 1;
                }
            }
            y(index, 0) = this->table[index].second[mid];
        }
    }
    // std::cout << std::endl;
    return y;
}
