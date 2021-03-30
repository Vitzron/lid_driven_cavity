#include <iostream>
#include "WENO.h"
#include "Eigen_Utilities.h"
using namespace Eigen;
using namespace std;

constexpr double SMALL = 1.0e-99;
constexpr double small = 1.0e-6;

void WENO::update(Eigen::ArrayXXd &phi,
        const Eigen::ArrayXXd &velo_x, const Eigen::ArrayXXd &velo_y,
        const double &dt, const double &dx, const double &dy,
        const size_t &num_gc) const
{
    Eigen::Index cols = phi.cols() - 2 * num_gc;
    Eigen::Index rows = phi.rows() - 2 * num_gc;
    if (cols != velo_x.cols() || rows != velo_x.rows() ||
        cols != velo_y.cols() || rows != velo_y.rows())
    {
        cerr << "INCOMPATIBLE shapes in WENO!" << endl;
        exit(EXIT_FAILURE);
    }

    Eigen_Utilities eu = Eigen_Utilities();
    ArrayXXd v, phi_x_m, phi_x_p, phi_y_m, phi_y_p;

    // phi_x_minus
    v = (phi.block(num_gc, 1, rows, phi.cols() - 2) - phi.block(num_gc, 0, rows, phi.cols() - 2)) / dx;
    set_phi_s_d(phi_x_m, v.middleCols(0, cols), v.middleCols(1, cols), v.middleCols(2, cols), v.middleCols(3, cols), v.middleCols(4, cols));
    // phi_x_plus
    v = (phi.block(num_gc, 2, rows, phi.cols() - 2) - phi.block(num_gc, 1, rows, phi.cols() - 2)) / dx;
    set_phi_s_d(phi_x_p, v.middleCols(4, cols), v.middleCols(3, cols), v.middleCols(2, cols), v.middleCols(1, cols), v.middleCols(0, cols));
    // phi_y_minus
    v = (phi.block(1, num_gc, phi.rows() - 2, cols) - phi.block(2, num_gc, phi.rows() - 2, cols)) / dy;
    set_phi_s_d(phi_y_m, v.middleRows(4, rows), v.middleRows(3, rows), v.middleRows(2, rows), v.middleRows(1, rows), v.middleRows(0, rows));
    // phi_y_plus
    v = (phi.block(0, num_gc, phi.rows() - 2, cols) - phi.block(1, num_gc, phi.rows() - 2, cols)) / dy;
    set_phi_s_d(phi_y_p, v.middleRows(0, rows), v.middleRows(1, rows), v.middleRows(2, rows), v.middleRows(3, rows), v.middleRows(4, rows));

    phi.block(num_gc, num_gc, rows, cols) -= dt * (eu.max_zero_x(velo_x) * phi_x_m + eu.min_zero_x(velo_x) * phi_x_p
        + eu.max_zero_x(velo_y) * phi_y_m + eu.min_zero_x(velo_y) * phi_y_p);
}

void WENO::set_phi_s_d(ArrayXXd &phi_s_d, const Eigen::ArrayXXd &v1,
        const Eigen::ArrayXXd &v2, const Eigen::ArrayXXd &v3,
        const Eigen::ArrayXXd &v4, const Eigen::ArrayXXd &v5) const
{
    ArrayXXd phi_x_1 = v1 / 3.0 - v2 * 7.0 / 6.0 + v3 * 11.0 / 6.0;
    ArrayXXd phi_x_2 = -v2 / 6.0 + v3 * 5.0 / 6.0 + v4 / 3.0;
    ArrayXXd phi_x_3 = v3 / 3.0 + v4 * 5.0 / 6.0 - v5 / 6.0;

    ArrayXXd S1 = 13.0 / 12.0 * (v1 - 2.0 * v2 + v3).square() 
        + 1.0 / 4.0 * (v1 - 4.0 * v2 + 3.0 * v3).square();
    ArrayXXd S2 = 13.0 / 12.0 * (v2 - 2.0 * v3 + v4).square()
        + 1.0 / 4.0 * (v2 - v4).square();
    ArrayXXd S3 = 13.0 / 12.0 * (v3 - 2.0 * v4 + v5).square()
        + 1.0 / 4.0 * (3.0 * v3 - 4.0 * v4 + v5).square();
    
    Eigen_Utilities eu = Eigen_Utilities();
    ArrayXXd epsilon = small * eu.max(eu.max(eu.max(eu.max(v1.square(), v2.square()), v3.square()), v4.square()), v5.square()) + SMALL;

    ArrayXXd alpha1 = 0.1 / (S1 + epsilon).square();
    ArrayXXd alpha2 = 0.6 / (S2 + epsilon).square();
    ArrayXXd alpha3 = 0.3 / (S3 + epsilon).square();
    ArrayXXd alpha = alpha1 + alpha2 + alpha3;
    phi_s_d = alpha1 / alpha * phi_x_1 + alpha2 / alpha * phi_x_2 + alpha3 / alpha * phi_x_3;
}