#include <latbolt.h>
#include <iostream>

using namespace std;
using namespace latbolt;

int main() {
    fluid f("test.h5", 100, 100, 1.6);

    for (int i = 0; i < 1000; i++) {
        f.write_rho();
        f.update();
    }

}
