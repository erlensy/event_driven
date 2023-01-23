#include <iostream>
#include "particles.h"

int main() {
    Particles ps{10, 10.0, 0.1, 1.0};
    ps.assert_no_overlap();
    ps.write_to_file("../data/skrt.txt");
}

