
#ifndef SPECTRA_JD_SYM_EIGS_BASE_H
#define SPECTRA_JD_SYM_EIGS_BASE_H

#include <Eigen/Core>
#include <vector>     // std::vector
#include <cmath>      // std::abs, std::pow
#include <algorithm>  // std::min
#include <stdexcept>  // std::invalid_argument
#include <iostream>

#include "Util/SelectionRule.h"
#include "Util/CompInfo.h"
#include "LinAlg/SearchSpace.h"
#include "LinAlg/RitzPairs.h"

namespace Spectra {
    template <typename Derived, typename OpType>
    class GDSymEigsBase
    {
    protected:
        using Index = Eigen::Index;
        using Scalar = typename OpType::Scalar;
        using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        const OpType& m_matrix_operator;  // object to conduct matrix operation,
        // e.g. matrix-vector product
        Index niter_ = 0;
        const Index m_number_eigenvalues;  // number of eigenvalues requested
        Index m_max_search_space_size;
        Index m_initial_search_space_size;
        Index m_correction_size;             // how many correction vectors are added in each iteration
        RitzPairs<Scalar> m_ritz_pairs;      // Ritz eigen pair structure
        SearchSpace<Scalar> m_search_space;  // search space

    private:
        CompInfo m_info = CompInfo::NotComputed;  // status of the computation

        void check_argument() const
        {
            if (m_number_eigenvalues < 1 || m_number_eigenvalues > m_matrix_operator.cols() - 1)
                throw std::invalid_argument("nev must satisfy 1 <= nev <= n - 1, n is the size of matrix");
        }

        void initialize()
        {
            // TODO better input validation and checks
            if (m_matrix_operator.cols() < m_max_search_space_size)
            {
                m_max_search_space_size = m_matrix_operator.cols();
            }
            if (m_matrix_operator.cols() < m_initial_search_space_size + m_correction_size)
            {
                m_initial_search_space_size = m_matrix_operator.cols() / 3;
                m_correction_size = m_matrix_operator.cols() / 3;
            }
        }

    public:
        //成员初始化列表
        GDSymEigsBase(OpType& op, OpType& M, Index nev, Index nvec_init, Index nvec_max) :
            m_matrix_operator(op),
			M_matrix_operator(M),//新加的M矩阵
            m_number_eigenvalues(nev), //要搜索的特征值个数 设定搜索空间和最大搜索空间，然后直接返回前nev个特征值
            m_max_search_space_size(nvec_max < op.rows() ? nvec_max : 10 * m_number_eigenvalues),
            m_initial_search_space_size(nvec_init < op.rows() ? nvec_init : 2 * m_number_eigenvalues),
            m_correction_size(m_number_eigenvalues) //!!
        {
            check_argument();
            initialize();
        }
        //委托构造函数
        GDSymEigsBase(OpType& op, OpType& M, Index nev) :
            GDSymEigsBase(op, M, nev, 2 * nev, 10 * nev) {
        }

        ///
        /// Sets the Maxmium SearchspaceSize after which is deflated
        ///
        void set_max_search_space_size(Index max_search_space_size)
        {
            m_max_search_space_size = max_search_space_size;
        }
        ///
        /// Sets how many correction vectors are added in each iteration
        ///
        void set_correction_size(Index correction_size)
        {
            m_correction_size = correction_size;
        }

        ///
        /// Sets the Initial SearchspaceSize for Ritz values
        ///
        void set_initial_search_space_size(Index initial_search_space_size)
        {
            m_initial_search_space_size = initial_search_space_size;
        }

        ///
        /// Virtual destructor
        ///
        virtual ~GDSymEigsBase() {}

        ///
        /// Returns the status of the computation.
        /// The full list of enumeration values can be found in \ref Enumerations.
        ///
        CompInfo info() const { return m_info; }

        ///
        /// Returns the number of iterations used in the computation.
        ///
        Index num_iterations() const { return niter_; }

        //在进行guess的时候，已经计算了特征值和特征向量
        Vector eigenvalues() const { return m_ritz_pairs.ritz_values().head(m_number_eigenvalues); }

        Matrix eigenvectors() const { return m_ritz_pairs.ritz_vectors().leftCols(m_number_eigenvalues); }

        Index compute(SortRule selection = SortRule::LargestMagn, Index maxit = 100,
            Scalar tol = 100 * Eigen::NumTraits<Scalar>::dummy_precision())
        {
            Derived& derived = static_cast<Derived&>(*this); //在基类中调用派生类的成员函数或访问派生类的成员变量
            Matrix intial_space = derived.setup_initial_search_space(selection); //在GD类里
            return compute_with_guess(intial_space, selection, maxit, tol);
        }

        Index compute_with_guess(const Eigen::Ref<const Matrix>& initial_space,
            SortRule selection = SortRule::LargestMagn,
            Index maxit = 100,
            Scalar tol = 100 * Eigen::NumTraits<Scalar>::dummy_precision())

        {
            m_search_space.initialize_search_space(initial_space);
            niter_ = 0;
            for (niter_ = 0; niter_ < maxit; niter_++)
            {
                bool do_restart = (m_search_space.size() > m_max_search_space_size);

                if (do_restart)
                {
                    m_search_space.restart(m_ritz_pairs, m_initial_search_space_size);
                }

                m_search_space.update_operator_basis_product(m_matrix_operator, M_matrix_operator);

                //在子空间中求解小型特征值问题  要改
                Eigen::ComputationInfo small_problem_info = m_ritz_pairs.compute_eigen_pairs(m_search_space); //在这里进行子空间里面的求解
                if (small_problem_info != Eigen::ComputationInfo::Success)
                {
                    m_info = CompInfo::NumericalIssue;
                    break;
                }
                m_ritz_pairs.sort(selection);  //要排序

                bool converged = m_ritz_pairs.check_convergence(tol, m_number_eigenvalues);
                if (converged)
                {
                    m_info = CompInfo::Successful;
                    break;
                }
                else if (niter_ == maxit - 1)
                {
                    m_info = CompInfo::NotConverging;
                    break;
                }
                //计算修正向量
                Derived& derived = static_cast<Derived&>(*this);
                Matrix corr_vect = derived.calculate_correction_vector();

                //加入新的搜索空间
                m_search_space.extend_basis(corr_vect);
            }
            return (m_ritz_pairs.converged_eigenvalues()).template cast<Index>().head(m_number_eigenvalues).sum();
        }
    };

}

#endif  
