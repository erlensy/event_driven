#include "includes.h"
#include "gas.h"

int main() {
    // {int N, double m, double r, double v0}
    Gas ps{300, 1.0, 0.01, 1.0};
    
    // {int frames, double timestep}
    ps.simulate(1000, 0.1);
}
