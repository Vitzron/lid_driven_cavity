/*
 * File:    V_BC.h
 * Desc:    header file of class V_BC
 */

#ifndef V_BC_
#define V_BC_

#include <Eigen/Dense>

class V_BC
{
    public:
        V_BC            () = default;
        ~V_BC           () = default;
        void non_slip_bc(Eigen::ArrayXXd &phi, const size_t &num_gc = 3) const;
};

#endif