#include "particles.h"

#include <random>

Particles::Particles(int N, double m, double r, double v0) {
    for (int i = 0; i < N; i++) {
        
        // create particle
        double theta = get_rand(0.0, 2.0 * M_PI);
        double vx = v0 * cos(theta);
        double vy = v0 * sin(theta);
        
        // check if overlap
        bool overlap = true;
        do {
            double x = get_rand(0.0 + r, 1.0 - r);
            double y = get_rand(0.0 + r, 1.0 - r);
            Particle p{x, y, vx, vy, r, m};

            int count = 0;
            for (Particle p2 : particles) {
                if (p2.dist_squared_to(p) > pow(2.0 * r, 2)) {
                    count += 1;
                }
                else {
                    break;
                }
            }
            if (count == particles.size()) {
                overlap = false;
                particles.push_back(p);
            }
        }
        while (overlap);
    }
}

void Particles::get_collisions() {
    // particle-particle collision
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            Particle p1 = particles[i];
            Particle p2 = particles[j];    
            vec delta_vel = p2.vel - p1.vel;
            vec delta_pos = p2.pos - p1.pos;
            double c1 = delta_vel * delta_pos;
            double c2 = delta_vel * delta_vel;
            double c3 = delta_pos * delta_pos;
            if (c1 < 0) {
                double d = c1 * c1 - (c2 * c2) * (c3 * c3 - p1.dist_squared_to(p2));
                if (d > 0) {
                    double delta_t = - (c1 + sqrt(d)) / c2;
                }
            }
        }
    }

    // particle-wall collision
    for (Particle p : particles) {
        // vertical wall
        if (p.vel.x != 0.0) {
            if (p.vel.x > 0) {
                double delta_t = (1.0 - p.r - p.pos.x) / p.vel.x;
            }
            else {
                double delta_t = (p.r - p.pos.x) / p.vel.x;
            }
        }

        // horizontal wall
        if (p.vel.y != 0.0) {
            if (p.vel.y > 0) {
                double delta_t = (1.0 - p.r - p.pos.y) / p.vel.y;
            }
            else {
                double delta_t = (p.r - p.pos.y) / p.vel.y;
            }
        }
}

void Particles::write_to_file(std::string filename) {
    std::ofstream out; out.open(filename);
    for (Particle p : particles) {
        out << p.pos.x << " " << p.pos.y << " " << p.vel.x << " " << p.vel.y << " " << p.r << " " << p.m << "\n";
    }
    out.close();
}

void Particles::assert_no_overlap() {
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double d = particles[i].dist_squared_to(particles[j]);
            assert(d > pow(particles[i].r + particles[j].r, 2));
        }
    }
}
