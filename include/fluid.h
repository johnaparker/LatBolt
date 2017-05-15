#ifndef GUARD_fluid_h
#define GUARD_fluid_h

#include "vec.h"
#include <string>
#include <array>

namespace latbolt {

class fluid {
public:
    fluid(std::string filename, int Nx, int Ny, double tau);

    void collide_and_stream();
    void update_macroscopic();
    void update();

    void write_rho();
    void write_U();

private:
    void init_file();

private:
    int Nx,Ny;
    double tau;             // relaxation
    double viscosity;       // kinematic shear viscosity
    double dx, dt;
    double t;

    constexpr static int Nv = 9;     // 9 velocities
    constexpr static int cs2 = 3;     // speed^-2
    constexpr static int cs4 = 9;     // speed^-4

    tensor3 fi;            // populations
    tensor3 fi_new;       // new populations (temporary storage)

    tensor rho;            // density
    tensor Ux;             // velocity field
    tensor Uy;

    std::array<const double, Nv> wi = {4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36};     //weights
    std::array<const int, Nv> cix = {0, 1, 0, -1, 0, 1, -1, -1, 1};     //weights
    std::array<const int, Nv> ciy = {0, 0, 1, 0, -1, 1, 1, -1, -1};     //weights

    h5cpp::h5file outFile;

};

}

#endif

