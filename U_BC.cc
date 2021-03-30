#include "U_BC.h"

void U_BC::non_slip_bc(Eigen::ArrayXXd &phi, const size_t &num_gc) const
{
    Eigen::Index rows = phi.rows();
    Eigen::Index cols = phi.cols();
    for (size_t i = 0; i != num_gc; ++i)
    {
        phi.row(i) = 2.0 - phi.row(2 * num_gc - i - 1);
        phi.row(rows - i - 1) = -phi.row(rows - 2 * num_gc + i);
        phi.col(i) = phi.col(2 * num_gc - i - 2);
        phi.col(cols - i - 1) = phi.col(cols - 2 * num_gc + i + 1);
    }
    phi.col(num_gc - 1) = 0.0;
    phi.col(cols - num_gc) = 0.0;
}