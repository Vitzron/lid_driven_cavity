#include <iostream>
#include <vector>
#include "Pressure.h"
#include "Eigen_Utilities.h"
using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

void Pressure::projection(Eigen::ArrayXXd &p, const Eigen::ArrayXXd &velo_x, const Eigen::ArrayXXd &velo_y,
        const double &dt, const double &dx, const double &dy) const
{
    if (p.rows() != velo_x.rows() || p.cols() != velo_y.cols() ||
        p.cols() - 1 != velo_x.cols() || p.rows() - 1 != velo_y.rows())
    {
        cerr << "INCOMPATIBLE shapes in Pressure!" << endl;
        exit(EXIT_FAILURE);
    }

    Index rows = p.rows() - 2;
    Index cols = p.cols() - 2;

    // right hand side
    VectorXd sv;
    rhs(sv, velo_x.block(1, 0, rows, cols + 1), velo_y.block(0, 1, rows + 1, cols), dt, dx, dy);

    // solver
    VectorXd vect_p = VectorXd::Zero(sv.size() + 1);
    vect_p.segment(0, sv.size()) = solver.solve(sv);
    Eigen_Utilities eu = Eigen_Utilities();
    p.block(1, 1, rows, cols) = eu.vect2C_array(vect_p, rows, cols);
}

void Pressure::set_solver(const Eigen::Index &rows, const Eigen::Index &cols, 
        const double &dx, const double &dy)
{
    std::vector<T> tripletList;
    tripletList.reserve(rows * cols * 5);
    const double X = 1.0 / dx / dx;
    const double Y = 1.0 / dy / dy;

    SpMat A = SpMat(rows * cols - 1, rows * cols - 1);
    Index i, j;
    //***********************first column*************************
    j = 0;
    // first row
    i = 0;
    tripletList.push_back(T(j * rows + i, j * rows + i, -X - Y));
    tripletList.push_back(T(j * rows + i, j * rows + i + 1, Y));
    tripletList.push_back(T(j * rows + i, j * rows + i + rows, X));
    // middle rows
    for (i = 1; i != rows - 1; ++i)
    {
        tripletList.push_back(T(j * rows + i, j * rows + i - 1, Y));
        tripletList.push_back(T(j * rows + i, j * rows + i, -X - Y - Y));
        tripletList.push_back(T(j * rows + i, j * rows + i + 1, Y));
        tripletList.push_back(T(j * rows + i, j * rows + i + rows, X));
    }
    // last row
    i = rows - 1;
    tripletList.push_back(T(j * rows + i, j * rows + i - 1, Y));
    tripletList.push_back(T(j * rows + i, j * rows + i, -X - Y));
    tripletList.push_back(T(j * rows + i, j * rows + i + rows, X));

    //*********************middle column************************
    for (j = 1; j != cols - 1; ++j)
    {
        for (i = 0; i != rows; ++i)
        {
            if (i == 0) // first row
            {
                tripletList.push_back(T(j * rows + i, j * rows + i - rows, X));
                tripletList.push_back(T(j * rows + i, j * rows + i, -X - X - Y));
                tripletList.push_back(T(j * rows + i, j * rows + i + 1, Y));
                tripletList.push_back(T(j * rows + i, j * rows + i + rows, X));
            }
            else if (i == rows - 1) // last rows
            {
                if (j != cols - 2)
                {
                    tripletList.push_back(T(j * rows + i, j * rows + i - rows, X));
                    tripletList.push_back(T(j * rows + i, j * rows + i - 1, Y));
                    tripletList.push_back(T(j * rows + i, j * rows + i, -X - X - Y));
                    tripletList.push_back(T(j * rows + i, j * rows + i + rows, X));
                }
                else // last second cols
                {
                    tripletList.push_back(T(j * rows + i, j * rows + i - rows, X));
                    tripletList.push_back(T(j * rows + i, j * rows + i - 1, Y));
                    tripletList.push_back(T(j * rows + i, j * rows + i, -X - X - Y));
                    // tripletList.push_back(T(j * rows + i, j * rows + i + rows, X));
                }
            }
            else // middle row
            {
                tripletList.push_back(T(j * rows + i, j * rows + i - rows, X));
                tripletList.push_back(T(j * rows + i, j * rows + i - 1,Y));
                tripletList.push_back(T(j * rows + i, j * rows + i, -X - X - Y - Y));
                tripletList.push_back(T(j * rows + i, j * rows + i + 1, Y));
                tripletList.push_back(T(j * rows + i, j * rows + i + rows, X));
            }
        }
    }

    //*********************last column************************
    j = cols - 1;
    // first row
    i = 0;
    tripletList.push_back(T(j * rows + i, j * rows + i - rows, X));
    tripletList.push_back(T(j * rows + i, j * rows + i, -X - Y));
    tripletList.push_back(T(j * rows + i, j * rows + i + 1, Y));
    // middle rows
    for (i = 1; i != rows - 2; ++i)
    {
        tripletList.push_back(T(j * rows + i, j * rows + i - rows, X));
        tripletList.push_back(T(j * rows + i, j * rows + i - 1, Y));
        tripletList.push_back(T(j * rows + i, j * rows + i, -X - Y - Y));
        tripletList.push_back(T(j * rows + i, j * rows + i + 1, Y));
    }
    // last second row
    i = rows - 2;
    tripletList.push_back(T(j * rows + i, j * rows + i - rows, X));
    tripletList.push_back(T(j * rows + i, j * rows + i - 1, Y));
    tripletList.push_back(T(j * rows + i, j * rows + i, -X - Y - Y));
    // tripletList.push_back(T(j * rows + i, j * rows + i + 1, Y));
    // last row
    /*
    i = rows - 1;
    tripletList.push_back(T(j * rows + i, j * rows + i - rows, X));
    tripletList.push_back(T(j * rows + i, j * rows + i - 1, Y));
    tripletList.push_back(T(j * rows + i, j * rows + i, -X - Y));
    */
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A.makeCompressed();

    solver.compute(A.block(0, 0, rows * cols - 1, rows * cols - 1));
}

void Pressure::rhs(Eigen::VectorXd &sv, const Eigen::ArrayXXd &velo_x, const Eigen::ArrayXXd &velo_y,
        const double &dt, const double &dx, const double &dy) const
{
    Index rows = velo_x.rows();
    Index cols = velo_y.cols();
    ArrayXXd as = ((velo_x.rightCols(cols) - velo_x.leftCols(cols)) / dx + (velo_y.topRows(rows) - velo_y.bottomRows(rows)) / dy) / dt;
    Map<VectorXd> ms(as.data(), as.size() - 1);
    sv = ms;
}