#include "fluid.h"
#include <vector>

using namespace std;

namespace latbolt {

fluid::fluid(std::string filename, int Nx, int Ny, double tau): Nx(Nx), Ny(Ny), tau(tau) {
    dx = 1;
    dt = 1;
    t = 0;
    viscosity = cs2*(tau - dt/2);

    outFile = h5cpp::h5file(filename, h5cpp::io::w);
    init_file();

    fi = tensor3(Nx, Ny, Nv);
    fi_new = tensor3(Nx, Ny, Nv);
    rho = tensor(Nx,Ny);
    Ux = tensor(Nx,Ny);
    Uy = tensor(Nx,Ny);


    rho.setConstant(1);
    Ux.setZero();
    Uy.setZero();

    for (int k = 0; k < Nv; k++)
        fi.chip(k,2) = wi[k] * rho;
    fi_new = fi;
}

void fluid::collide_and_stream() {
    static double omega = dt/tau;
    static double omega_p = 1 - omega;

    for (int i = 1; i < Nx-1; i++) {
        for (int j = 1; i < Ny-1; j++) {
            for (int k = 0; k < Nv; k++) {
                // equilibrium
                double fi_eq = wi[k]*rho(i,j) * (1 + cs2*(Ux(i,j)*cix[k] + Uy(i,j)*ciy[k])
                                                   + 0.5*cs4*pow((Ux(i,j)*cix[k] + Uy(i,j)*ciy[k]), 2)
                                                   - 0.5*cs2*(Ux(i,j)*Ux(i,j) + Uy(i,j)*Uy(i,j)) );

                // compute stress tensor here if needed
                
                // collision
                double fi_star = omega_p * fi(i,j,k) + omega * fi_eq;

                // streaming
                fi_new(i+cix[k], j+ciy[k], k) = fi_star;

            }
        }
    }
    fi = fi_new;
}

void fluid::update_macroscopic() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; i < Ny; j++) {
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

    // apply boundary conditions here
    
    t += t;

    // force computation here
}

void fluid::write_rho() {
    auto dset = outFile.open_dataset("rho");
    dset.append(rho.data());
}

void fluid::write_U() {
    auto dset = outFile.open_dataset("Ux");
    dset.append(Ux.data());

    dset = outFile.open_dataset("Uy");
    dset.append(Uy.data());
}

void fluid::init_file() {
    vector<hsize_t> dims = {hsize_t(Nx),hsize_t(Ny),1};
    vector<hsize_t> max_dims = {hsize_t(Nx),hsize_t(Ny),h5cpp::inf};
    vector<hsize_t> chunk_dims = dims;
    h5cpp::dspace ds(dims, max_dims, chunk_dims, false);

    auto dset = outFile.create_dataset("rho", h5cpp::dtype::Double, ds);
    dset = outFile.create_dataset("Ux", h5cpp::dtype::Double, ds);
    dset = outFile.create_dataset("Uy", h5cpp::dtype::Double, ds);
}

}
