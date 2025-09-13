#include "numerics/integrate/solver.hpp"

#include "numerics/integrate/rk.hpp"
#include "numerics/integrate/bdf.hpp"

namespace wfem {

std::shared_ptr<ODE::Solver> ODE::Solver::create(Type tp, ODE::Function f, double t_0, const Eigen::VectorXd& y_0)
{
    switch (tp)
    {
    case Type::ODE45:
    case Type::RK45:

        return std::make_shared<RK45>(f, t_0, y_0);
    
    case Type::ODE23:
    case Type::RK23:

        return std::make_shared<RK23>(f, t_0, y_0);
    
    case Type::DOP853:

        return std::make_shared<DOP853>(f, t_0, y_0);

    case Type::BDF:
    case Type::ODE15s:

        return std::make_shared<BDF>(f, t_0, y_0);

    default:
        break;
    }
    return nullptr;
}

ODE::Solver::Solver(ODE::Function f, double t_0, const Eigen::VectorXd& y_0, int order)
:   m_max_step(std::numeric_limits<double>::infinity()),
    m_atol(Eigen::VectorXd::Constant(y_0.size(), 1.e-6)),
    m_rtol(1.e-3),
    m_t(t_0),
    m_y(y_0),
    m_evaluate_times(0),
    m_ode_function(f)
{
    m_f = this->compute_y_dot(t_0, y_0);
    m_h = this->select_initial_step(m_t, m_y, order);
}

double ODE::Solver::select_initial_step(double t_0, const Eigen::VectorXd& y_0, int order) const
{
    Eigen::VectorXd scale = m_atol + y_0.cwiseAbs() * m_rtol;
    Eigen::VectorXd scale_inv = scale.cwiseInverse();

    Eigen::VectorXd f_0 = compute_y_dot(t_0, y_0);
    double d_0 = RMS_norm(y_0.cwiseProduct(scale_inv));
    double d_1 = RMS_norm(f_0.cwiseProduct(scale_inv));

    double h_0 = (d_0 < 1e-5 || d_1 < 1e-5) ? 1e-6 : 
        (0.01 * d_0 / d_1);

    Eigen::VectorXd y_1 = y_0 + h_0 * f_0;
    Eigen::VectorXd f_1 = compute_y_dot(t_0 + h_0, y_1);
    double d_2 = RMS_norm((f_1 - f_0).cwiseProduct(scale_inv)) / h_0;

    double h_1 = (d_1 <= 1e-15 && d_2 <= 1e-15) ? std::max(1e-6, h_0 * 1e-3) : std::pow(0.01 / std::max(d_1, d_2), 1. / ((double)order + 1.));

    return std::min(100. * h_0, h_1);
}

} // namespace wfem