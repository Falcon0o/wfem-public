#pragma once

#include <memory>
#include <vector>

#include <functional>

#include "Eigen/Dense"

namespace wfem {
namespace ODE {



using Function = std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>;

class Solver
{
public:
    enum ErrCode {
        OK = 0,
        TOO_SMALL_STEP
    };
    virtual ErrCode step(double t_next) = 0;
    enum class Type {
        ODE45,      RK45,
        ODE23,      RK23,
        DOP853,
        ODE15s,     BDF
    };
    static std::shared_ptr<Solver> create(Type, ODE::Function f, double t_0, const Eigen::VectorXd& y_0);

    const Eigen::VectorXd& y() const { return m_y; }
    Eigen::VectorXd& y() { return m_y; }

    const Eigen::VectorXd& f() const { return m_f; }
    Eigen::VectorXd& f() { return m_f; }


    virtual int order() const = 0;

protected:

    Solver(ODE::Function f, double t_0, const Eigen::VectorXd& y_0, int order);

    Eigen::VectorXd compute_y_dot(double t, const Eigen::VectorXd& y) const {
        ++m_evaluate_times;
        return m_ode_function(t, y);
    }

    double select_initial_step(double t_0, const Eigen::VectorXd& y_0, int order) const;

    static double RMS_norm(const Eigen::VectorXd& x)
    {
        return std::sqrt(x.squaredNorm() / x.size());
    }

    
    const double& t() const { return m_t; }
    double& t() { return m_t; }

    const double& h() const { return m_h; }
    double& h() { return m_h; }

    const double& h_prev() const { return m_h_prev; }
    double& h_prev() { return m_h_prev; }

    const Eigen::VectorXd& atol() const { return m_atol; }
    const double& rtol() const { return m_rtol; }

    const double& max_step() const { return m_max_step; }

private:

    double          m_max_step;
    Eigen::VectorXd m_atol;
    double          m_rtol;
    
    ODE::Function   m_ode_function;
    mutable int     m_evaluate_times;

    double          m_t;

    Eigen::VectorXd m_y;
    Eigen::VectorXd m_f;

    double m_h;
    double m_h_prev;
};


} // namespace ODE
} // namespace wfem
