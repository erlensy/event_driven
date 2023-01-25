#include <iostream>
#include "particles.h"

int main() {
    Particles ps{20, 11.0, 0.05, 0.1};
    //Particles ps{};
    ps.assert_no_overlap();
    ps.write_to_file("../data/skrt.txt");
    ps.simulate();
}

