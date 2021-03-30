/*
 * File:    U_BC.h
 * Desc:    header file of class U_BC
 */

#ifndef U_BC_
#define U_BC_

#include <Eigen/Dense>

class U_BC
{
    public:
        U_BC            () = default;
        ~U_BC           () = default;
        void non_slip_bc(Eigen::ArrayXXd &phi, const size_t &num_gc = 3) const;
};

#endif