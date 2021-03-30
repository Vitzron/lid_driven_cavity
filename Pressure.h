/*
 * File:    Pressure.h
 * Desc:    header file of class Pressure
 */

#ifndef PRESSURE_
#define PRESSURE_

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Pressure
{
    public:
        /*
        Pressure                (const Eigen::Index &rows, const Eigen::Index &cols, 
                                 const double &dx, const double &dy);
        */
        Pressure                () {};
        ~Pressure               () = default;
        void set_solver         (const Eigen::Index &rows, const Eigen::Index &cols, 
                                 const double &dx, const double &dy);
        void projection         (Eigen::ArrayXXd &p, const Eigen::ArrayXXd &velo_x, const Eigen::ArrayXXd &velo_y,
                                 const double &dt, const double &dx, const double &dy) const;
    private:
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        void rhs                (Eigen::VectorXd &sv, const Eigen::ArrayXXd &velo_x, const Eigen::ArrayXXd &velo_y,
                                 const double &dt, const double &dx, const double &dy) const;
};

#endif