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
    private:
        // number of particles
        int N;
        
        // vector which holds pointers to each particle
        std::vector<Particle*> particles;
        
        // priority queue with collisions
        std::priority_queue<Collision> pq;
    
        // collision parameter
        double epsilon;
    
    public:
        // constructors
        Gas(int N, double epsilon, std::vector<double>& m, std::vector<double>& r, std::vector<double>& v0,
            double x_min, double x_max, double y_min, double y_max);
    
       
        // helper functions for constructor
        void make_particle(double x, double y, double vx, double vy, double r, double m);
    
        void make_particles_with_no_overlap(int n, std::vector<double> v0,
                                            std::vector<double> r, std::vector<double> m,
                                            double x_min, double x_max,
                                            double y_min, double y_max);
        
        // main function which simulates gas
        void simulate(int frames, double timestep, double end);
    
        // move all particle dt forward in time
        void move_forward(double dt);
        
        // wrapper for calling correct collision case
        void manage_collision(int collision_type, double t);
        
        // check if a collision is valid
        bool valid_collision(const Collision* c);

        // wrapper for calling both get_collisions_walls
        // and get_collisions_pp
        void get_collisions(Particle* p, double t);
    
        void get_all_collisions(double t);
        
        // add collisions with particle and walls to pq
        void get_collisions_walls(Particle* p, double t);
    
        // add collisions between particle and all other
        // particles to pq
        void get_collisions_pp(Particle* p, double t);
        
        // resolve and add collisions to the pq,
        // for each collision case
        void manage_p_vw_collision(Particle* p, double t);
        void manage_p_hw_collision(Particle* p, double t);
        void manage_p_p_collision(Particle* p1, Particle* p2, double t);
        
        // used for testing if particles have overlap
        void assert_no_overlap();
    
        // save to file functions
        void save_particles(std::ofstream& out);
        std::ofstream initialize_ofstream(std::string filename, double timestep, int frames);
    
        
        // destructor which deallocates memory
        ~Gas();
};

