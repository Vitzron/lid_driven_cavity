/*
 * File: WENO.h
 * Desc: header file of class WENO
 */

#ifndef WENO_
#define WENO_

#include <Eigen/Dense>

class WENO
{
    public:
        WENO            () = default;
        ~WENO           () = default;
        void update     (Eigen::ArrayXXd &phi,
                         const Eigen::ArrayXXd &velo_x, const Eigen::ArrayXXd &velo_y,
                         const double &dt, const double &dx, const double &dy,
                         const size_t &num_gc = 3) const;
    private:
        void set_phi_s_d(Eigen::ArrayXXd &phi_s_d, const Eigen::ArrayXXd &v1,
                         const Eigen::ArrayXXd &v2, const Eigen::ArrayXXd &v3,
                         const Eigen::ArrayXXd &v4, const Eigen::ArrayXXd &v5) const;  
};

#endif