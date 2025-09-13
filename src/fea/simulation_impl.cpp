#include "fea/simulation_impl.hpp"

#include "parsers/abaqus/inp_parser.hpp"
#include "parsers/abaqus/inp_report.hpp"

#include "fea/materials/material.hpp"

#include "fea/element.hpp"

#include <iostream>

#include "base/assert.h"

#include <omp.h>
#include <iostream>
#include <chrono>



#include "parallel/mpi.hpp"
#include "parallel/mpi/comm.hpp"
#include "base/logger.hpp"

#include "numerics/sparse_matrices/elimination_tree.hpp"

#include "numerics/linear_solver.hpp"

#include <fstream>

namespace wfem {
namespace fea {

Simulation::Impl::Impl()
:   m_pardiso_solver(nullptr),
    m_cudss_solver(nullptr),
    m_amgx_solver(nullptr),
    m_num_equations(0),
    m_num_variables(0)
{
}

Simulation::Impl::~Impl()
{
    if (m_pardiso_solver) {
        delete m_pardiso_solver;
        m_pardiso_solver = nullptr;
    }

    if (m_cudss_solver) {
        delete m_cudss_solver;
        m_cudss_solver = nullptr;
    }

    if (m_amgx_solver) {
        delete m_amgx_solver;
        m_amgx_solver = nullptr;
    }
}
void Simulation::Impl::load_inp_file(const char* inp_file)
{
    wfem::abaqus::InpParser inp_parser;
    inp_parser.parse("/home/wangxinyu/workspace/myself/wfem/data/conrod.inp");
    auto report = inp_parser.get_report();

    m_x0.assign(report->node_coords(), report->node_coords() + report->n_nodes() * 3);
    m_u.assign(report->n_nodes() * 3, 0);

    // std::cout << report->n_nodes() * 3 << std::endl;
    const auto& elements = report->elements();
    for (auto iter : elements)
    {
        m_element_factory.create(m_elements, iter.first, iter.second.ids.size(), iter.second.nodes.data());
    }
}

void Simulation::Impl::run()
{
    m_num_equations = m_num_variables = m_u.size();

    m_unpacked_lhs.set_nrows(m_num_equations);
    m_unpacked_lhs.set_ncols(m_num_variables);


    struct TimePoint {
        TimePoint(const char* lb)
        :   label(lb),
            time_point(std::chrono::high_resolution_clock::now())
        {}
        
        std::chrono::high_resolution_clock::time_point time_point;
        const char* label;
    };

    std::vector<TimePoint> time_points;
    time_points.emplace_back("初始化");

    {
        int nnz_lhs = 0;
        int nnz_rhs = 0;

        for (Element* element : m_elements)
        {
            int dof = element->dof(); 
            nnz_lhs += dof * dof;
            nnz_rhs += dof;
        }
        m_unpacked_lhs.reserve(nnz_lhs);
        m_unpacked_rhs.reserve(nnz_rhs);
    }
    
    for (Element* element : m_elements)
    {
        element->evaluate(m_unpacked_lhs.rows(), m_unpacked_lhs.cols(), m_unpacked_rhs.rows);
    }

    time_points.emplace_back("计算稀疏矩阵拓扑结构");

    // rho 7.850 10^3 kg / m^3 = 7.85 e-9 t / mm^3  
    // K 2.1e11 kg / (m * s^2) = 2.1e5 t / (mm * s^2)
    Material<LinearElastic, Isotropy> mat(7.85e-09, 210000., 0.3);
    
    Element::Context ctx{
        .x0 = m_x0.data(),
        .u = m_u.data(),
        .v = nullptr,
        .material = &mat
    };

    m_unpacked_lhs.set_zero();
    m_unpacked_rhs.set_zero();

    ::omp_set_num_threads(8);
    #pragma omp parallel for schedule(static)
    for (int ii = 0; ii < m_elements.size(); ++ii)
    {
        Element* element = m_elements[ii];
        element->evaluate(ctx, m_unpacked_lhs.values().data(), m_unpacked_rhs.values.data());
    }

    time_points.emplace_back("完成线性方程组左端项和右端项计算");

#if 0
    /// 从coo映射待csr
    std::vector<int> matrix_map(unpacked_lhs.size());
    for (int ii = 0; ii < unpacked_lhs.size(); ++ii)
    {
        int jj = pardiso_lhs.find(unpacked_lhs.row(ii), unpacked_lhs.col(ii));
        WFEM_ASSERT(jj >= 0 && jj < pardiso_lhs.nnz());
        matrix_map[ii] = jj;
    }
#endif
    DenseVector<double> rhs(m_num_equations);
    rhs.set_zero();
    // time_points.emplace_back("更新CSR稀疏矩阵");

    for (int ii = 0; ii < m_unpacked_rhs.size(); ++ii)
    {
        rhs[m_unpacked_rhs.rows[ii]] += m_unpacked_rhs.values[ii];
    }
    for (int ii = 0; ii < rhs.nrows(); ++ii)
    {
        rhs[ii] = 1e-5;
    }

    DenseVector<double> pardiso_sol(m_num_equations);
    this->solve_pardiso(pardiso_sol, rhs);

    DenseVector<double> cudss_sol(m_num_equations);
    this->solve_cudss(cudss_sol, rhs);

    DenseVector<double> amgx_sol(m_num_equations);
    this->solve_amgx(amgx_sol, rhs);

    std::ofstream res_csv("zzz.csv", std::ios_base::trunc);

    for (int ii = 0; ii < m_num_variables; ++ii)
    {
        res_csv << pardiso_sol[ii] << ", "
                << cudss_sol[ii] << ", "
                << amgx_sol[ii] << ", "
                << "\n";
    }
#if 0


    for (int ii = 1; ii < time_points.size(); ++ii)
    {
        float time_used = float(std::chrono::duration_cast<std::chrono::microseconds>(time_points[ii].time_point - time_points[ii - 1].time_point).count()) / 1000;
        std::cout << ii << ": [" << time_points[ii].label << "] " << time_used << "ms" << std::endl; 
    }

#endif
}

void Simulation::Impl::solve_pardiso(DenseVector<double>& x, const DenseVector<double>& b)
{
    if (m_pardiso_solver == nullptr)
    {
        m_pardiso_solver = new LinearSolver::Pardiso;
    }

    m_pardiso_lhs = m_unpacked_lhs;

    m_pardiso_solver->analyze(m_pardiso_lhs);

    m_pardiso_solver->factorize(m_pardiso_lhs);

    m_pardiso_solver->solve(x, b);
}

void Simulation::Impl::solve_cudss(DenseVector<double>& x, const DenseVector<double>& b)
{
    if (m_cudss_solver == nullptr)
    {
        m_cudss_solver = new LinearSolver::CudssOnHost;
    }

    m_cudss_lhs = m_unpacked_lhs;

    m_cudss_solver->analyze(m_cudss_lhs);

    m_cudss_solver->factorize(m_cudss_lhs);

    m_cudss_solver->solve(x, b);
}

void Simulation::Impl::solve_amgx(DenseVector<double>& x, const DenseVector<double>& b)
{
    if (m_amgx_solver == nullptr)
    {
        m_amgx_solver = new LinearSolver::AMGX;
    }

    m_amgx_lhs = m_unpacked_lhs;

    m_amgx_solver->setup(m_amgx_lhs, false);
    m_amgx_solver->solve_with_0_initial_guess(x, b);
}
} // namespace fea
} // namespace wfem