#include "numerics/integrate/solver.hpp"

#include <iostream>

int main()
{
    // create(Type, ODE::Function f, double t_0, const Eigen::VectorXd& y_0);

    /// std::exp(std::sin(t))
    auto f = [](double t, const Eigen::VectorXd& y)->Eigen::VectorXd
    {
        Eigen::VectorXd y_dot(1);
        y_dot[0] = y[0] * std::cos(t);
        return y_dot;
    };

    Eigen::VectorXd y0(1);
    y0[0] = 1.;
    double t = 0;
    auto rk45 = wfem::ODE::Solver::create(wfem::ODE::Solver::Type::RK45, f, t, y0);
    auto rk23 = wfem::ODE::Solver::create(wfem::ODE::Solver::Type::RK23, f, t, y0);
    auto dop853 = wfem::ODE::Solver::create(wfem::ODE::Solver::Type::DOP853, f, t, y0);
    auto bdf = wfem::ODE::Solver::create(wfem::ODE::Solver::Type::BDF, f, t, y0);

    const double h_sample = 1;
    while (t <= 100)
    {
        t += h_sample;
        rk45->step(t);
        rk23->step(t);
        dop853->step(t);
        bdf->step(t);

        double target = std::exp(std::sin(t));
        std::cout << t << ", " << target << ", " << rk45->y()[0] - target << ", " << rk23->y()[0] - target << ", " << dop853->y()[0] - target << ", " << bdf->y()[0] - target << std::endl;
    }
    return 0;
}