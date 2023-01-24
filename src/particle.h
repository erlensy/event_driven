#pragma once

#include "utils.h"

class Particle {
    public:
        // position
        vec pos;

        // velocity
        vec vel;

        // radius
        double r;

        // mass
        double m;

        // constructor
        Particle(double x, double y, double vx, double vy, double r, double m);

        // helper functions
        double dist_squared_to(Particle& particle);
};
