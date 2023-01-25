#pragma once

#include <queue>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include "particle.h"
#include "utils.h"

const std::map<std::string, int> collision_type_map {
    {"particle-particle", 0},
    {"particle-vertical_wall", 1},
    {"particle-horizontal_wall", 2}
};

struct Collision {
    double t;
    Particle* p1;
    Particle* p2;
    int collision_type;
    int p1_count;
    int p2_count;
};

bool operator<(const Collision& c1, const Collision& c2);

class Particles {
    private:
        // number of particles
        int N;

        std::vector<Particle*> particles;

        std::priority_queue<Collision> pq;

    public:
        // initialize random non-overlapping particles
        Particles(int N, double m, double r, double v0);
        Particles();

        void move_forward(double delta_t);
        void initialize();

        void get_collisions(Particle* p);
        void simulate();

        void assert_no_overlap();
        void write_to_file(std::string filename);

        ~Particles(); 
};

