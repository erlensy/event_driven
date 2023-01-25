#pragma once
#include "includes.h"
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

        // move particle forward in time dt
        void move_forward(double dt);
    
        // collision time and resolve collision functions
        // particle with particle collision
        double get_collision_time_with(Particle* p);
        void resolve_collision_with(Particle* p);
    
        // particle with horizontal wall collision
        double get_collision_time_with_hw();
        void resolve_collision_with_hw();

        // particle with vertical wall collision
        double get_collision_time_with_vw();
        void resolve_collision_with_vw();

        // helper functions
        double dist_squared_to(Particle* p);
        friend std::ostream& operator<<(std::ostream& os, Particle* p);
};
