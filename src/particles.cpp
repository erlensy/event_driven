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

void Particles::write_to_file(std::string filename) {
    std::ofstream out; out.open(filename);
    for (Particle p : particles) {
        out << p.x << " " << p.y << " " << p.vx << " " << p.vy << " " << p.r << " " << p.m << "\n";
    }
    out.close();
}

void Particles::assert_no_overlap() {
    for (int i = 0; i < particles.size(); i++) {
        for (int j = i + 1; j < particles.size(); j++) {
            double d = particles[i].dist_squared_to(particles[j]);
            assert(d > pow(particles[i].r + particles[j].r, 2));
        }
    }
}
