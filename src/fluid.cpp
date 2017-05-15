#include "fluid.h"

using namespace std;

namespace latbolt {

fluid::fluid(std::string filename, int Nx, int Ny, double tau): Nx(Nx), Ny(Ny), tau(tau) {
    dx = 1;
    dt = 1;

    outFile = h5cpp::h5file(filename, h5cpp::io::w);

    fi = tensor3(Nx, Ny, Nv);
    fi_prev = tensor3(Nx, Ny, Nv);
    rho = tensor(Nx,Ny);
    Ux = tensor(Nx,Ny);
    Uy = tensor(Nx,Ny);

    fi.setZero();
    fi_prev.setZero();
    rho.setZero();
    Ux.setZero();
    Uy.setZero();
}

void fluid::collide_and_stream() {
    static double omega = dt/tau;
    static double omega_p = 1 - omega;

    for (int i = 0; i != Nx; i++) {
        for (int j = 0; i != Ny; j++) {
            for (int k = 0; k != Nv; k++) {
                double fi_eq = wi[k]*rho(i,j) * (1 + cs2*(Ux(i,j)*cix[k] + Uy(i,j)*ciy[k])
                                                   + 0.5*cs4*pow((Ux(i,j)*cix[k] + Uy(i,j)*ciy[k]), 2)
                                                   - 0.5*cs2*pow((Ux(i,j)*Ux(i,j) + Uy(i,j)*Uy(i,j)), 2) );
                double fi_star = omega_p * fi(i,j,k) + omega * fi_eq;

                fi(i+cix[k], j+ciy[k], k) = fi_star;

            }
        }
    }
}

void fluid::update_macroscopic() {
    for (int i = 0; i != Nx; i++) {
        for (int j = 0; i != Ny; j++) {
            rho(i,j) = fi(i,j,0) + fi(i,j,1) + fi(i,j,2) + fi(i,j,3) + fi(i,j,4)
                       + fi(i,j,5) + fi(i,j,6) + fi(i,j,7) + fi(i,j,8) + fi(i,j,9);

            Ux(i,j) = 1/rho(i,j) * (fi(i,j,1) + fi(i,j,5) + fi(i,j,8)
                                    - fi(i,j,4) - fi(i,j,6) - fi(i,j,7));

            Uy(i,j) = 1/rho(i,j) * (fi(i,j,2) + fi(i,j,5) + fi(i,j,6)
                                    - fi(i,j,4) - fi(i,j,7) - fi(i,j,8));
        }
    }
}

void fluid::update() {
    update_macroscopic();
    collide_and_stream();
}

}
