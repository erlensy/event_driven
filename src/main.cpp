#include <iostream>
#include "particles.h"

int main() {
    Particles ps{300, 1.0, 0.01, 1.0};
    //Particles ps{};
    ps.assert_no_overlap();
    ps.write_to_file("../data/skrt.txt");
    ps.simulate();
}

