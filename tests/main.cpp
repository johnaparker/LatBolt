#include <latbolt.h>
#include <iostream>

using namespace std;
using namespace latbolt;

int main() {
    fluid f("test.h5", 300, 150, 1.0);

    for (int i = 0; i < 1000; i++) {
        f.write_U();
        f.update();
    }

}
