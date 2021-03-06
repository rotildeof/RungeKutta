#ifndef _RungeKutta_
#define _RungeKutta_

#include "RungeKuttaSimul.hh"

template <std::size_t N>
class RungeKutta : public RungeKuttaSimul<N> {
    using func_detail = std::function<double(double, const double*)>;
public:
    RungeKutta();
    ~RungeKutta() {};
    void AssignFunction(func_detail func);
};

template <std::size_t N>
RungeKutta<N>::RungeKutta() {}

template <std::size_t N>
void RungeKutta<N>::AssignFunction(func_detail func) {
    for (int i = 0; i < N - 1; ++i) {
        this->funcs[i] = [=](double x, const double* y) {
            return y[i + 1];
        };
    }
    this->funcs[N - 1] = std::move(func);
}

#endif
