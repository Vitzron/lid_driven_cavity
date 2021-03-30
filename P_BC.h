/*
 * File:    P_BC.h
 * Desc:    header file of class P_BC
 */

#ifndef P_BC_
#define P_BC_

#include <Eigen/Dense>

class P_BC
{
    public:
        P_BC        () = default;
        ~P_BC       () = default;
        void bc     (Eigen::ArrayXXd &phi) const;
        void bc     (Eigen::ArrayXXd &phi, Eigen::ArrayXXd &phi_x, Eigen::ArrayXXd &phi_y) const;
};

#endif