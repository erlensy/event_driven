#pragma once

#include <cmath>

class Particle {
    public:
        // position
        double x;
        double y;

            // velocity
        double vx;
        double vy;

        // radius
        double r;

        // mass
        double m;

        // constructor
        Particle(double x, double y, double vx, double vy, double r, double m);

        // helper functions
        double dist_squared_to(Particle& particle);
};
