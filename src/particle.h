#pragma once
#include "includes.h"
#include "utils.h"

class Particle {
    private:
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

    public:
        // constructor
        Particle(double x, double y, double vx, double vy, double r, double m);

        // move particle forward in time dt
        void move_forward(double dt);
    
        // collision time and resolve collision functions
        // particle with particle collision
        double get_collision_time_with(Particle* p);
        void resolve_collision_with(Particle* p, double epsilon);
    
        // particle with horizontal wall collision
        double get_collision_time_with_hw();
        void resolve_collision_with_hw(double epsilon);

        // particle with vertical wall collision
        double get_collision_time_with_vw();
        void resolve_collision_with_vw(double epsilon);

        // helper functions
        double dist_squared_to(Particle* p);
    
        // let << and class Gas access private variables
        friend std::ostream& operator<<(std::ostream& os, Particle* p);
        friend class Gas;
};
