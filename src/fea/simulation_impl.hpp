#pragma once

#include "fea/simulation.hpp"

#include "fea/element_factory.hpp"

#include "numerics/sparse_matrix.hpp"
#include "numerics/linear_solver.hpp"


namespace wfem {

namespace fea {

class Simulation::Impl 
{
public:

    using UnpackedMatrix = SparseMatrix::COO<int, double, SparseMatrix::BaseZero>;

    struct UnpackedVector {
        std::vector<int> rows;
        std::vector<double> values;

        int size() const { return rows.size(); }
        void set_zero() {
            values.assign(rows.size(), 0);
        }

        void reserve(int s) { 
            values.reserve(s);
            rows.reserve(s);
        }
    };

public:

    Impl();
    ~Impl();

    void load_inp_file(const char* inp_file);
    void run();

private:

    void solve_pardiso(DenseVector<double>& x, const DenseVector<double>& b);
    void solve_cudss(DenseVector<double>& x, const DenseVector<double>& b);
    void solve_amgx(DenseVector<double>& x, const DenseVector<double>& b);

private:

    ElementFactory m_element_factory;

    std::vector<Element*> m_elements;

    std::vector<double> m_x0;
    std::vector<double> m_u;

    UnpackedMatrix  m_unpacked_lhs;
    UnpackedVector  m_unpacked_rhs;

    int m_num_equations;
    int m_num_variables;

    /// 这里用指针是为了便于统计时间
    LinearSolver::Pardiso*                  m_pardiso_solver;
    LinearSolver::Pardiso::MATRIX           m_pardiso_lhs;

    LinearSolver::CudssOnHost*              m_cudss_solver;
    LinearSolver::CudssOnHost::MATRIX       m_cudss_lhs;
    
    LinearSolver::AMGX*                     m_amgx_solver;
    LinearSolver::AMGX::MATRIX              m_amgx_lhs;
};
} // namespace fea
} // namespace wfem