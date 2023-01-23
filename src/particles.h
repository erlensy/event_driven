#pragma once

#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include "particle.h"
#include "utils.h"

class Particles {
    private:
        // number of particles
        int N;

        std::vector<Particle> particles;

    public:
        // initialize random non-overlapping particles
        Particles(int N, double m, double r, double v0);

        void assert_no_overlap();
        void write_to_file(std::string filename);
};

