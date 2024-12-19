
#ifndef SPECTRA_DAVIDSON_SYM_EIGS_SOLVER_H
#define SPECTRA_DAVIDSON_SYM_EIGS_SOLVER_H

#include <Eigen/Core>

#include "JDSymEigsBase.h"
#include "Util/SelectionRule.h"

namespace Spectra {
///
template <typename OpType>
class DavidsonSymEigsSolver : public JDSymEigsBase<DavidsonSymEigsSolver<OpType>, OpType>
{
private:
    using Index = Eigen::Index;
    using Scalar = typename OpType::Scalar;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Vector m_diagonal;

public:
    DavidsonSymEigsSolver(OpType& op, Index nev, Index nvec_init, Index nvec_max) :
        JDSymEigsBase<DavidsonSymEigsSolver<OpType>, OpType>(op, nev, nvec_init, nvec_max)
    {
        m_diagonal.resize(this->m_matrix_operator.rows());
        for (Index i = 0; i < op.rows(); i++)
        {
            m_diagonal(i) = op(i, i);
        }
    }

    DavidsonSymEigsSolver(OpType& op, Index nev) :
        DavidsonSymEigsSolver(op, nev, 2 * nev, 10 * nev) {}
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

	/// 用了DPR修正向量 
    Matrix calculate_correction_vector() const
    {
        const Matrix& residues = this->m_ritz_pairs.residues();
        const Vector& eigvals = this->m_ritz_pairs.ritz_values();
        Matrix correction = Matrix::Zero(this->m_matrix_operator.rows(), this->m_correction_size);
        for (Index k = 0; k < this->m_correction_size; k++)
        {
            Vector tmp = eigvals(k) - m_diagonal.array();
            correction.col(k) = residues.col(k).array() / tmp.array();
        }
        return correction;
    }
};

}  

#endif  