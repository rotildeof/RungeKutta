#ifndef _RungeKuttaSimul_
#define _RungeKuttaSimul_

#include <valarray>
#include <array>
#include <vector>
#include <valarray>
#include <functional>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using func_detail = std::function<double(double, double*)>;

template <std::size_t N>
class RungeKuttaSimul {
public:
    RungeKuttaSimul();
    ~RungeKuttaSimul() {};
    void SetInitValues(double x, double* y);
    void SetStep(double step) { h = step; };
    void SetMaximumX(double x_max) { x_max_ = x_max;  is_max_set = true; };
    void AssignFunction(int i_func, std::function<double(double, double*)> func);
    void Solve(double x_max);
    void WriteFile(const char* filename);

    int64_t GetSize() { return x_.size(); };
    double GetValueX(double i) { return x_[i]; };
    double GetValueY(double i, double k) { return y_[i][k]; };

protected:
    std::array<func_detail, N> funcs; // funcs[N]
    std::array<std::valarray<double>, 4> k_; // k_[4][N]
    std::vector<std::valarray<double> > y_; // y_[][N]
    std::vector<double> x_; // x_[]
    double h = 0.001;
    double x_min_;
    double x_max_;
    bool is_max_set = false;
    bool is_min_set = false;
    bool is_first_solve = true;
};

template <std::size_t N>
RungeKuttaSimul<N>::RungeKuttaSimul() {
    std::for_each(std::begin(k_), std::end(k_), [](std::valarray<double>& k) {k.resize(N);});
    std::for_each(std::begin(y_), std::end(y_), [](std::valarray<double>& y) {y.resize(N);});
}

template <std::size_t N>
void RungeKuttaSimul<N>::SetInitValues(double x, double* y) {
    x_.push_back(x);
    x_min_ = x;
    is_min_set = true;
    y_.emplace_back(y, N);
}

template <std::size_t N>
void RungeKuttaSimul<N>::AssignFunction(int i_func, func_detail func) {
    funcs[i_func] = std::move(func);
}

template <std::size_t N>
void RungeKuttaSimul<N>::Solve(double x_max) {
    if (!is_min_set) {
        std::cout << "Error in Function \"Solve()\" : initial value is not given. " << std::endl;
        return;
    }
    int i_now = x_.size();
    int64_t rep = (x_max - x_[i_now - 1]) / h;
    if (is_max_set && is_first_solve) {
        int64_t size_capacity = (x_max_ - x_min_) / h + 10;// + 10 for buffer
        y_.reserve(size_capacity);
        x_.reserve(size_capacity);
    }
    is_first_solve = false;

    std::cout << "==== Performing Runge Kutta ====" << std::endl;
    std::cout << "init value x = " << x_[i_now - 1] << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "init value y[" << i << "] = " << y_[i_now - 1][i] << std::endl;
    }
    std::cout << "Last x value = " << x_max << std::endl;
    std::cout << "Number of repetition = " << rep << std::endl;
    for (int64_t i = i_now; i <= rep + i_now - 1; ++i) {
        std::transform(std::begin(funcs), std::end(funcs), std::begin(k_[0]),
                       [&](func_detail f) {
                           return h * f(x_[i - 1], &y_[i - 1][0]);
                       });
        std::transform(std::begin(funcs), std::end(funcs), std::begin(k_[1]),
                       [&](func_detail f) {
                           std::valarray<double> y_plus_k0_over2 = y_[i - 1] + k_[0] / 2;
                           return h * f(x_[i - 1] + h / 2, &y_plus_k0_over2[0]);
                       });
        std::transform(std::begin(funcs), std::end(funcs), std::begin(k_[2]),
                       [&](func_detail f) {
                           std::valarray<double> y_plus_k1_over2 = y_[i - 1] + k_[1] / 2;
                           return h * f(x_[i - 1] + h / 2, &y_plus_k1_over2[0]);
                       });
        std::transform(std::begin(funcs), std::end(funcs), std::begin(k_[3]),
                       [&](func_detail f) {
                           std::valarray<double> y_plus_k2 = y_[i - 1] + k_[2];
                           return h * f(x_[i - 1] + h, &y_plus_k2[0]);
                       });
        x_.emplace_back(x_[i - 1] + h);
        y_.emplace_back(y_[i - 1] + (k_[0] + 2 * k_[1] + 2 * k_[2] + k_[3]) / 6.);
    }
    // std::cout << "==== RungeKutta End ====" << std::endl;
}

template <std::size_t N>
void RungeKuttaSimul<N>::WriteFile(const char* filename) {
    std::ofstream ofile(filename);
    int64_t size = this->GetSize();
    ofile << std::setprecision(14) << std::scientific << std::fixed;
    for (int64_t i = 0; i < size; ++i) {
        ofile << x_[i] << " ";
        for (int j = 0; j < N; j++) {
            ofile << y_[i][j] << " ";
        }
        ofile << "\n";
    }
    ofile << std::flush;
}
#endif
