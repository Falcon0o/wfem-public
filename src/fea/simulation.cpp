
#include "fea/simulation_impl.hpp"

#include "numerics/linear_solvers/amgx.hpp"

namespace wfem {
namespace fea {

Simulation::Simulation()
:   m_impl(nullptr)
{
}

Simulation::~Simulation()
{
}

Simulation& Simulation::initialize(const Parameters& parameters, int* argc, char*** argv)
{
    LinearSolver::AMGX::initialize();
    m_impl = new Simulation::Impl;
    return *this;
}

Simulation& Simulation::finalize()
{
    if (m_impl)
    {
        delete m_impl;
        m_impl = nullptr;
    }

    LinearSolver::AMGX::finalize();
    return *this;
}






Simulation& Simulation::load_inp_file(const char* inp_file)
{
    m_impl->load_inp_file(inp_file);
    return *this;
}

Simulation& Simulation::run()
{
    m_impl->run();
    return *this;
}

} // namespace fea
} // namespace wfem