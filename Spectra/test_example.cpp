#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <iostream>

using namespace Spectra;

// M = diag(1, 2, ..., 10)
class MyDiagonalTen
{
public:
    using Scalar = double;  // A typedef named "Scalar" is required
    int rows() const { return 10; }
    int cols() const { return 10; }
    // y_out = M * x_in
    void perform_op(const double* x_in, double* y_out) const
    {
        for (int i = 0; i < rows(); i++)
        {
            y_out[i] = x_in[i] * (i + 1);
        }
    }
};

int main()
{
    MyDiagonalTen op;
    DavidsonSymEigsSolver<MyDiagonalTen> eigs(op, 3, 6);
    eigs.init();
    eigs.compute(SortRule::LargestAlge);
    if (eigs.info() == CompInfo::Successful)
    {
        Eigen::VectorXd evalues = eigs.eigenvalues();
        std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    }

    return 0;
}