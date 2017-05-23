#include "fluid.h"
#include <vector>
#include <iostream>

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

    //rho(50,50) = 3;
    for (int k = 0; k < Nv; k++)
        fi.chip(k,2) = wi[k] * rho;
    fi_new = fi;

    //for (int k = 0; k != Nv; k++)
        //fi(50,50,k) += 0.1;
}

void fluid::collide_and_stream() {
    static double omega = dt/tau;
    static double omega_p = 1 - omega;

    for (int i = 0; i < Nx; i++) {
        for (int j = 1; j < Ny-1; j++) {
            for (int k = 0; k < Nv; k++) {
                // equilibrium
                double fi_eq = wi[k]*rho(i,j) * (1 + cs2*(Ux(i,j)*cix[k] + Uy(i,j)*ciy[k])
                                                   + 0.5*cs4*pow((Ux(i,j)*cix[k] + Uy(i,j)*ciy[k]), 2)
                                                   - 0.5*cs2*(Ux(i,j)*Ux(i,j) + Uy(i,j)*Uy(i,j)) );

                // compute stress tensor here if needed
                
                // collision
                double fi_star = omega_p * fi(i,j,k) + omega * fi_eq;

                unsigned int xi = (Nx + i + cix[k]) % Nx;
                unsigned int yi = j + ciy[k];
                // streaming
                //fi_new(i+cix[k], j+ciy[k], k) = fi_star;
                fi_new(xi,yi, k) = fi_star;

            }
        }
    }
}

void fluid::update_macroscopic() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            //rho(i,j) = fi(i,j,0) + fi(i,j,1) + fi(i,j,2) + fi(i,j,3) + fi(i,j,4)
                       //+ fi(i,j,5) + fi(i,j,6) + fi(i,j,7) + fi(i,j,8);

            //Ux(i,j) = 1/rho(i,j) * (fi(i,j,1) + fi(i,j,5) + fi(i,j,8)
                                    //- fi(i,j,4) - fi(i,j,6) - fi(i,j,7));

            //Uy(i,j) = 1/rho(i,j) * (fi(i,j,2) + fi(i,j,5) + fi(i,j,6)
                                    //- fi(i,j,4) - fi(i,j,7) - fi(i,j,8));

            rho(i,j) = 0;
            Ux(i,j) = 0;
            Uy(i,j) = 0;
            for (int k = 0; k < Nv; k++) {
                rho(i,j) += fi(i,j,k);
                Ux(i,j) += cix[k] * fi(i,j,k);
                Uy(i,j) += ciy[k] * fi(i,j,k);
            }
        }
    }
}

void fluid::update() {

    update_macroscopic();
    collide_and_stream();

    // apply boundary conditions here
    double u_in = 0.003*3;
    double rho_in = 1.01;
    double rho_out = 1.0;
    for (int j = 1; j < Ny-1; j++) {
        // velocity BC
        double rho_wall = 1/(1-u_in) * (fi_new(0,j,0) + fi_new(0,j,2) + fi_new(0,j,4)
                 + 2*(fi_new(0,j,3) + fi_new(0,j,6) + fi_new(0,j,7)) );
        fi_new(0,j,1) = fi_new(0,j,3);
        fi_new(0,j,5) = fi_new(0,j,7) - 0.5*(fi_new(0,j,2) - fi_new(0,j,4)) + 1.0/6*rho_wall*u_in;
        fi_new(0,j,8) = fi_new(0,j,6) + 0.5*(fi_new(0,j,2) - fi_new(0,j,4)) + 1.0/6*rho_wall*u_in;

        // pressure BC
        // left
        double ux = 1 - 1/(rho_in) * (fi_new(0,j,0) + fi_new(0,j,2) + fi_new(0,j,4)
                                 + 2*(fi_new(0,j,3) + fi_new(0,j,6) + fi_new(0,j,7)) );
        //fi_new(0,j,1) = fi_new(0,j,3) + 2/3*rho_in*ux;
        //fi_new(0,j,5) = fi_new(0,j,7) - 0.5*(fi_new(0,j,2) - fi_new(0,j,4)) + 1.0/6*rho_in*ux;
        //fi_new(0,j,8) = fi_new(0,j,6) + 0.5*(fi_new(0,j,2) - fi_new(0,j,4)) + 1.0/6*rho_in*ux;

        // right
        ux = 1 - 1/(rho_out) * (fi_new(Nx-1,j,0) + fi_new(Nx-1,j,2) + fi_new(Nx-1,j,4)
                                 + 2*(fi_new(Nx-1,j,1) + fi_new(Nx-1,j,5) + fi_new(Nx-1,j,8)) );
        fi_new(Nx-1,j,3) = fi_new(Nx-1,j,1) + 2/3*rho_out*ux;
        fi_new(Nx-1,j,6) = fi_new(Nx-1,j,5) - 0.5*(fi_new(Nx-1,j,2) - fi_new(Nx-1,j,4)) + 1.0/6*rho_out*ux;
        fi_new(Nx-1,j,7) = fi_new(Nx-1,j,8) + 0.5*(fi_new(Nx-1,j,2) - fi_new(Nx-1,j,4)) + 1.0/6*rho_out*ux;

        //rho(0,j) = rho_wall;
    }
    //fi = fi_new;
    //t += dt;
    //return;
    
    // top and bottom walls
    for (int i = 1; i < Nx-1; i++) {
        fi_new(i,1,2) = fi_new(i,0,4);
        fi_new(i,1,5) = fi_new(i-1,0,7);
        fi_new(i,1,6) = fi_new(i+1,0,8);
        fi_new(i,0,4) = wi[4]; 
        fi_new(i-1,0,7) = wi[7]; 
        fi_new(i+1,0,8) = wi[8]; 

        fi_new(i,Ny-2,4) = fi_new(i,Ny-1,2);
        fi_new(i,Ny-2,7) = fi_new(i+1,Ny-1,5);
        fi_new(i,Ny-2,8) = fi_new(i-1,Ny-1,6);
        fi_new(i,Ny-1,2) = wi[2]; 
        fi_new(i+1,Ny-1,5) = wi[5]; 
        fi_new(i-1,Ny-1,6) = wi[6]; 
    }

    // bottom left corner
    fi_new(0,1,2) = fi_new(0,0,4);
    fi_new(0,1,6) = fi_new(1,0,8);
    fi_new(0,0,4) = wi[4];
    fi_new(1,0,8) = wi[8];

    // bottom right corner
    fi_new(Nx-1,1,2) = fi_new(Nx-1,0,4);
    fi_new(Nx-1,1,5) = fi_new(Nx-2,0,7);
    fi_new(Nx-1,0,4) = wi[4];
    fi_new(Nx-2,0,7) = wi[7];

    // top left corner
    fi_new(0,Ny-2,4) = fi_new(0,Ny-1,2);
    fi_new(0,Ny-2,7) = fi_new(1,Ny-1,5);
    fi_new(0,Ny-1,2) = wi[2];
    fi_new(1,Ny-1,5) = wi[5];

    // top right corner
    fi_new(Nx-1,Ny-2,4) = fi_new(Nx-1,Ny-1,2);
    fi_new(Nx-1,Ny-2,8) = fi_new(Nx-2,Ny-1,6);
    fi_new(Nx-1,Ny-1,2) = wi[2];
    fi_new(Nx-2,Ny-1,6) = wi[6];


    // box
    // left and right
    int ia = 40, ib = 60;
    int ja = 65, jb = 85;
    for (int j = ja; j < jb; j++) {
        fi_new(ia,j,3) = fi_new(ia+1,j,1);
        fi_new(ia,j,6) = fi_new(ia+1,j-1,8);
        fi_new(ia,j,7) = fi_new(ia+1,j+1,5);

        fi_new(ia+1,j,1) = wi[1]; 
        fi_new(ia+1,j-1,8) = wi[8]; 
        fi_new(ia+1,j+1,5) = wi[5]; 

        fi_new(ib,j,1) = fi_new(ib-1,j,3);
        fi_new(ib,j,8) = fi_new(ib-1,j+1,6);
        fi_new(ib,j,5) = fi_new(ib-1,j-1,7);

        fi_new(ib-1,j,3) = wi[3]; 
        fi_new(ib-1,j+1,6) = wi[6]; 
        fi_new(ib-1,j-1,7) = wi[7]; 
    }
    // top and bottom
    for (int i = ia; i < ib; i++) {
        fi_new(i,ja,4) = fi_new(i,ja+1,2);
        fi_new(i,ja,7) = fi_new(i+1,ja+1,5);
        fi_new(i,ja,8) = fi_new(i-1,ja+1,6);

        fi_new(i,ja+1,2) = wi[2]; 
        fi_new(i+1,ja+1,5) = wi[5]; 
        fi_new(i-1,ja+1,6) = wi[6]; 

        fi_new(i,jb,2) = fi_new(i,jb-1,4);
        fi_new(i,jb,6) = fi_new(i+1,jb-1,8);
        fi_new(i,jb,5) = fi_new(i-1,jb-1,7);

        fi_new(i,jb-1,4) = wi[4]; 
        fi_new(i+1,jb-1,8) = wi[8]; 
        fi_new(i-1,jb-1,7) = wi[7]; 
    }


    // copy xmax to xmin, xmin to xmax
    
    fi = fi_new;
    t += dt;

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
    vector<hsize_t> dims = {hsize_t(Nx),hsize_t(Ny),0};
    vector<hsize_t> max_dims = {hsize_t(Nx),hsize_t(Ny),h5cpp::inf};
    vector<hsize_t> chunk_dims = {hsize_t(Nx),hsize_t(Ny),1};
    h5cpp::dspace ds(dims, max_dims, chunk_dims, false);

    auto dset = outFile.create_dataset("rho", h5cpp::dtype::Double, ds);
    dset = outFile.create_dataset("Ux", h5cpp::dtype::Double, ds);
    dset = outFile.create_dataset("Uy", h5cpp::dtype::Double, ds);
}

}
