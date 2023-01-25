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

        // collision count
        int count;

        // constructor
        Particle(double x, double y, double vx, double vy, double r, double m);

        double get_collision_time_with(Particle* p);

        double get_collision_time_with_horizontal_wall();

        double get_collision_time_with_vertical_wall();

        void move_forward(double delta_t);

        void resolve_collision_with(Particle* p);

        void resolve_collision_with_horizontal_wall();

        void resolve_collision_with_vertical_wall();
        // helper functions
        double dist_squared_to(Particle* p);
};
