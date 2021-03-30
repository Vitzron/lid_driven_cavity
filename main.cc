#include <iostream>
#include <fstream>
#include "Grid.h"
#include "Grid_Utilities.h"
#include "WENO.h"
#include "Pressure.h"
#include "U_BC.h"
#include "V_BC.h"
#include "P_BC.h"
using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
    // set simulation condition
    size_t rows = atoi(argv[1]);
    size_t cols = atoi(argv[1]);
    double Re = atof(argv[2]);

    const size_t num_gc = 3;
    Array2i dims(cols, rows);
    Array2d lo(0., 0.);
    Array2d hi(1.0, 1.0);
    Grid grid(dims, lo, hi, num_gc);
    double dx = grid.get_dx()(0);
    double dy = grid.get_dx()(1);
    cout << "**************************grid**************************"
         << endl << grid 
         << "**************************grid**************************"
         << endl;
    
    Grid_Utilities gu(grid);
    ArrayXXd X = gu.get_X_ct();
    ArrayXXd Y = gu.get_Y_ct();
    ArrayXXd X_fx = gu.get_X_fx();
    ArrayXXd Y_fx = gu.get_Y_fx();
    ArrayXXd X_fy = gu.get_X_fy();
    ArrayXXd Y_fy = gu.get_Y_fy();

    ofstream output;
    const char *X_file = "X.dat";
    const char *Y_file = "Y.dat";

    output.open(X_file, std::ios::out);
    if (!output)
    {
        cerr << "CANNOT open file " << X_file << "!" << endl;
        exit(EXIT_FAILURE);
    }
    output << X_fy;
    output.flush();
    output.close();
    cout << "SUCCESSFULLY write X into file " << X_file << "!" << endl;
    output.open(Y_file, std::ios::out);
    if (!output)
    {
        cerr << "CANNOT open file " << Y_file << "!" << endl;
        exit(EXIT_FAILURE);
    }
    output << Y_fx;
    output.flush();
    output.close();
    cout << "SUCCESSFULLY write Y into file " << Y_file << "!" << endl;

    ArrayXXd u = ArrayXXd::Zero(rows + 2 * num_gc, cols + 2 * num_gc - 1);
    ArrayXXd v = ArrayXXd::Zero(rows + 2 * num_gc - 1, cols + 2 * num_gc);
    ArrayXXd p = ArrayXXd::Zero(rows + 2, cols + 2);

    WENO weno = WENO();
    Pressure pressure;
    pressure.set_solver(rows, cols, dx, dy);
    U_BC u_bc = U_BC();
    V_BC v_bc = V_BC();
    P_BC p_bc = P_BC();
    u_bc.non_slip_bc(u);

    const double CFL = 0.1;
    const double dt_0 = CFL * min(dx, dy);
    double dt = dt_0;
    cout << "dx: " << dx << "; dy: " << dy << "; dt: " << dt << endl;

    double T = 0.0;
    size_t i = 0;

    ArrayXXd velo_x, velo_y;
    ArrayXXd u_t, v_t;
    ArrayXXd u_0, v_0;
    ArrayXXd div;
    const double tolerance = 1.0e-8;
    ArrayXXd u_p, v_p;
    u_p = u + 1.0;
    v_p = v + 1.0;

    while ((u_p - u).abs().maxCoeff() > tolerance || (v_p - v).abs().maxCoeff() > tolerance)
    {
        cout << "step: " << i << "; dt: " << dt << "; T: " << T << "; variation of velo_x: " << (u - u_p).block(1, 1, rows, cols - 1).abs().maxCoeff() 
             << "; variation of  velo_y: " << (v - v_p).block(1, 1, rows - 1, cols).abs().maxCoeff() << endl;
        
        u_p = u;
        v_p = v;

        u_0 = u;
        v_0 = v;

        for (size_t substep = 0; substep != 1; ++substep)
        {
            u_t = u;
            v_t = v;

            //*********************advection*********************
            // u
            velo_x = u_t.block(3, 3, rows, cols - 1);
            velo_y = (v_t.block(2, 3, rows, cols - 1) + v_t.block(2, 4, rows, cols - 1) + v_t.block(3, 3, rows, cols - 1) + v_t.block(3, 4, rows, cols - 1)) / 4.0;
            weno.update(u, velo_x, velo_y, dt, dx, dy);
            u_bc.non_slip_bc(u);

            // v
            velo_x = (u_t.block(3, 2, rows - 1, cols) + u_t.block(3, 3, rows - 1, cols) + u_t.block(4, 2, rows - 1, cols) + u_t.block(4, 3, rows - 1, cols)) / 4.0;
            velo_y = v_t.block(3, 3, rows - 1, cols);
            weno.update(v, velo_x, velo_y, dt, dx, dy);
            v_bc.non_slip_bc(v);

            //*********************nonadvection (i)*********************
            velo_x = u.block(2, 2, rows + 2, cols + 1);
            velo_y = v.block(2, 2, rows + 1, cols + 2);
            u.block(3, 3, rows, cols - 1) += dt / Re * ((velo_x.block(1, 2, rows, cols - 1) + velo_x.block(1, 0, rows, cols - 1) - 2.0 * velo_x.block(1, 1, rows, cols - 1)) / dx / dx 
                + (velo_x.block(0, 1, rows, cols - 1) + velo_x.block(2, 1, rows, cols - 1) - 2.0 * velo_x.block(1, 1, rows, cols - 1)) / dy / dy);
            u_bc.non_slip_bc(u);
            v.block(3, 3, rows - 1, cols) += dt / Re * ((velo_y.block(1, 2, rows - 1, cols) + velo_y.block(1, 0, rows - 1, cols) - 2.0 * velo_y.block(1, 1, rows - 1, cols)) / dx / dx 
                + (velo_y.block(0, 1, rows - 1, cols) + velo_y.block(2, 1, rows - 1, cols) - 2.0 * velo_y.block(1, 1, rows - 1, cols)) / dy / dy);
            v_bc.non_slip_bc(v);

            /*********************nonadvection (ii)*********************/
            velo_x = u.block(2, 2, rows + 2, cols + 1);
            velo_y = v.block(2, 2, rows + 1, cols + 2);
            pressure.projection(p, velo_x, velo_y, dt, dx, dy);
            p_bc.bc(p);

            // update velocity
            u.block(3, 3, rows, cols - 1) -= dt * (p.block(1, 2, rows, cols - 1) - p.block(1, 1, rows, cols - 1)) / dx;
            u_bc.non_slip_bc(u);
            v.block(3, 3, rows - 1, cols) -= dt * (p.block(1, 1, rows - 1, cols) - p.block(2, 1, rows - 1, cols)) / dy;
            v_bc.non_slip_bc(v);

            if (substep == 1)
            {
                u = 0.75 * u_0 + 0.25 * u;
                v = 0.75 * v_0 + 0.25 * v;
            }
            if (substep == 2)
            {
                u = 1.0 / 3.0 * u_0 + 2.0 / 3.0 * u;
                v = 1.0 / 3.0 * v_0 + 2.0 / 3.0 * v;
            }
        }
        velo_x = u.block(2, 2, rows + 2, cols + 1);
        velo_y = v.block(2, 2, rows + 1, cols + 2);
        div = (velo_x.block(1, 1, rows, cols) - velo_x.block(1, 0, rows, cols)) / dx + (velo_y.block(0, 1, rows, cols) - velo_y.block(1, 1, rows, cols)) / dy;
        cout << "max of div: " << div.abs().maxCoeff() << endl;

        dt = min(dt_0, CFL / (velo_x.abs().maxCoeff() / dx + velo_y.abs().maxCoeff() / dy + 2.0 / Re / max(dx, dy) / max(dx, dy)));
        T += dt;
        ++i;
        if (T > 100.0) break;
    }

    const char *u_file = "u.dat";
    const char *v_file = "v.dat";
    output.open(u_file, std::ios::out);
    if (!output)
    {
        cerr << "CANNOT open file " << u_file << "!" << endl;
        exit(EXIT_FAILURE);
    }
    output << u;
    output.flush();
    output.close();
    cout << "SUCCESSFULLY write u into file " << u_file << "!" << endl;
    output.open(v_file, std::ios::out);
    if (!output)
    {
        cerr << "CANNOT open file " << v_file << "!" << endl;
        exit(EXIT_FAILURE);
    }
    output << v;
    output.flush();
    output.close();
    cout << "SUCCESSFULLY write v into file " << v_file << "!" << endl;

    return 0;
}