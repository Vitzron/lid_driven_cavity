#include "P_BC.h"

void P_BC::bc(Eigen::ArrayXXd &phi) const
{
    Eigen::Index rows = phi.rows();
    Eigen::Index cols = phi.cols();

    phi.row(0) = phi.row(1);
    phi.row(rows - 1) = phi.row(rows - 2);
    phi.col(0) = phi.col(1);
    phi.col(cols - 1) = phi.col(cols - 2);
}

void P_BC::bc(Eigen::ArrayXXd &phi, Eigen::ArrayXXd &phi_x, Eigen::ArrayXXd &phi_y) const
{
    Eigen::Index rows = phi.rows();
    Eigen::Index cols = phi.cols();

    // phi
    phi.row(0) = phi.row(1);
    phi.row(rows - 1) = phi.row(rows - 2);
    phi.col(0) = phi.col(1);
    phi.col(cols - 1) = phi.col(cols - 2);

    // phi_x
    phi_x.row(0) = phi_x.row(1);
    phi_x.row(rows - 1) = phi_x.row(rows - 2);
    phi_x.col(0) = 0.0;
    phi_x.col(cols - 1) = 0.0;

    // phi_y
    phi_y.row(0) = 0.0;
    phi_y.row(rows - 1) = 0.0;
    phi_y.col(0) = phi_y.col(1);
    phi_y.col(cols - 1) = phi_y.col(cols - 2);
}