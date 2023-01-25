#pragma once
#include "includes.h"
#include "particle.h"

const std::map<std::string, int> collision_type_map {
    // maps type of collision to an integer
    // used in struct Collision
    {"particle-particle", 0},
    {"particle-vertical_wall", 1},
    {"particle-horizontal_wall", 2}
};

struct Collision {
    // time of collision
    double t;
    
    // particles involved in collision
    // p2 is a NULL ptr if particle-wall collision
    Particle* p1;
    Particle* p2;
    
    // refers to collision_type_map
    int collision_type;
    
    // collision counts of each particle at time t
    int p1_count;
    int p2_count;
};

// overload compare operator because of collision struct
// is used in the priority queue.
bool operator<(const Collision& c1, const Collision& c2);

class Gas {
    public:
        // number of particles
        int N;
        
        // vector which holds pointers to each particle
        std::vector<Particle*> particles;
        
        // priority queue with collisions
        std::priority_queue<Collision> pq;

        // constructors
        Gas(int N, double m, double r, double v0);
        
        // main function which simulates gas
        void simulate(int frames, double timestep);
    
        // calculate all possible collisions and
        // add them to the priority queue
        void initialize();
    
        // move all particle dt forward in time
        void move_forward(double dt);

        void get_collisions(Particle* p, double t);
    
        void get_collisions_walls(Particle* p, double t);
        void get_collisions_pp(Particle* p, double t);

        void manage_p_vw_collision(Particle* p, double t);
        void manage_p_hw_collision(Particle* p, double t);
        void manage_p_p_collision(Particle* p1, Particle* p2, double t);
        

        // helper functions
        void assert_no_overlap();
    
        ~Gas();
};

