
#ifndef SPECTRA_SYM_GEIGS_SOLVER_H
#define SPECTRA_SYM_GEIGS_SOLVER_H

#include "SymEigsBase.h"
#include "Util/GEigsMode.h"
#include "MatOp/internal/SymGEigsCholeskyOp.h"
#include "MatOp/internal/SymGEigsRegInvOp.h"

namespace Spectra {

// Empty class template
template <typename OpType, typename BOpType, GEigsMode Mode>
class SymGEigsSolver
{};
// Partial specialization for mode = GEigsMode::Cholesky
template <typename OpType, typename BOpType>
class SymGEigsSolver<OpType, BOpType, GEigsMode::Cholesky> :
    public SymEigsBase<SymGEigsCholeskyOp<OpType, BOpType>, IdentityBOp>
{
private:
    using Scalar = typename OpType::Scalar;
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    using ModeMatOp = SymGEigsCholeskyOp<OpType, BOpType>;
    using Base = SymEigsBase<ModeMatOp, IdentityBOp>;

    const BOpType& m_Bop;

public:
    SymGEigsSolver(OpType& op, BOpType& Bop, Index nev, Index ncv) :
        Base(ModeMatOp(op, Bop), IdentityBOp(), nev, ncv),
        m_Bop(Bop)
    {}

    /// \cond

    Matrix eigenvectors(Index nvec) const override
    {
        Matrix res = Base::eigenvectors(nvec);
        Vector tmp(res.rows());
        const Index nconv = res.cols();
        for (Index i = 0; i < nconv; i++)
        {
            m_Bop.upper_triangular_solve(&res(0, i), tmp.data());
            res.col(i).noalias() = tmp;
        }

        return res;
    }

    Matrix eigenvectors() const override
    {
        return SymGEigsSolver<OpType, BOpType, GEigsMode::Cholesky>::eigenvectors(this->m_nev);
    }

    /// \endcond
};
template <typename OpType, typename BOpType>
class SymGEigsSolver<OpType, BOpType, GEigsMode::RegularInverse> :
    public SymEigsBase<SymGEigsRegInvOp<OpType, BOpType>, BOpType>
{
private:
    using Index = Eigen::Index;

    using ModeMatOp = SymGEigsRegInvOp<OpType, BOpType>;
    using Base = SymEigsBase<ModeMatOp, BOpType>;

public:
    SymGEigsSolver(OpType& op, BOpType& Bop, Index nev, Index ncv) :
        Base(ModeMatOp(op, Bop), Bop, nev, ncv)
    {}
};

}  // namespace Spectra

#endif  
