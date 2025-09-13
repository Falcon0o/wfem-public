#include "wfem/numerics/ode/solver.hpp"

#include "wfem/simulink/model.hpp"
#include <iostream>

#include <cmath>

int main(int argc, char** argv)
{
    wfem::Simulink::Model sim_model;
    sim_model.deserialize_json("/home/wangxinyu/workspace/myself/wfem/sample/20250707/sim_model.json");
    return 0;
    sim_model.show_direct_graph_puml("/home/wangxinyu/workspace/myself/wfem/sample/20250707/new.puml");
    auto func = [](double* y, double t) {
        y[0] = std::exp(std::sin(t));
    };

    auto solver = wfem::ODE::Solver::create();
    solver->set_relative_error(1e-9);

    int n_cont_states = sim_model.n_continuous_states();
    std::vector<double> y_end(n_cont_states);
    std::vector<double> tmp(n_cont_states);
    std::vector<double> y_0(n_cont_states);

    double t = 0;
    double h = 1;
    sim_model.initialize_conditions(y_0.data());

    std::cout << "时间, 数值结, 解析解\n";
    std::cout << t << ", " << y_0[0] << ", " << y_0[0] << std::endl;;
    while (t < 100)
    {
        double t_end = t + h;
        solver->run(1, y_end.data(), y_0.data(), t, t_end,
            [&sim_model](double* y_dot, double t, const double* y) {
                sim_model.evaluate(y_dot, t, y);
            }
        );

        y_0 = y_end; 
        t = t_end;

        func(tmp.data(), t);        
        std::cout << t << ", " << y_0[0] << ", " << tmp[0] << std::endl;;
    }
    return 0;
}

