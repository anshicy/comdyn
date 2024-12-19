
#ifndef SPECTRA_DAVIDSON_SYM_EIGS_SOLVER_H
#define SPECTRA_DAVIDSON_SYM_EIGS_SOLVER_H

#include <Eigen/Core>

#include "JDSymEigsBase.h"
#include "Util/SelectionRule.h"

namespace Spectra {
    ///
    template <typename OpType>
    class GDavidsonSymEigsSolver : public GDSymEigsBase<DavidsonSymEigsSolver<OpType>, OpType>
    {
    private:
        using Index = Eigen::Index;
        using Scalar = typename OpType::Scalar;
        using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        Vector m_1_diagonal;
        Vector m_diagonal;

    public:
        GDavidsonSymEigsSolver(OpType& op, OpType& M, Index nev, Index nvec_init, Index nvec_max) :
            GDSymEigsBase<DavidsonSymEigsSolver<OpType>, OpType>(op, M,nev, nvec_init, nvec_max)
        {
            m_diagonal.resize(this->m_matrix_operator.rows());
			m_1_diagonal.resize(this->m_matrix_operator.rows());
            for (Index i = 0; i < op.rows(); i++)
            {
                m_diagonal(i) = op(i, i);
				m_1_diagonal(i) = M(i, i);
				//这里开两个数组，存刚度阵和质量阵的对角线元素
            }
        }

        GDavidsonSymEigsSolver(OpType& op, OpType& M, Index nev) :
            GDavidsonSymEigsSolver(op, M, nev, 2 * nev, 10 * nev) {
        }
        Matrix setup_initial_search_space(SortRule selection) const
        {
            std::vector<Eigen::Index> indices_sorted = argsort(selection, m_diagonal);

            Matrix initial_basis = Matrix::Zero(this->m_matrix_operator.rows(), this->m_initial_search_space_size);

            for (Index k = 0; k < this->m_initial_search_space_size; k++)
            {
                Index row = indices_sorted[k];
                initial_basis(row, k) = 1.0;
            }
            return initial_basis;
        }

        /// 用了DPR修正向量 适用于广义特征值问题!!!
        Matrix calculate_correction_vector() const
        {
            const Matrix& residues = this->m_ritz_pairs.residues();
            const Vector& eigvals = this->m_ritz_pairs.ritz_values(); //记得看看是什么
            Matrix correction = Matrix::Zero(this->m_matrix_operator.rows(), this->m_correction_size);
            for (Index k = 0; k < this->m_correction_size; k++)
            {
                Vector tmp = eigvals(k) * m_1_diagonal.array() - m_diagonal.array();
                correction.col(k) = residues.col(k).array() / tmp.array();
            }
            return correction;
        }
    };

}

#endif  