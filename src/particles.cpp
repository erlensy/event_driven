#include "particles.h"

#include <random>

bool operator<(const Collision& c1, const Collision& c2) {
    if (c1.t > c2.t) {
        return true;
    }
    return false;
}

Particles::Particles(int N, double m, double r, double v0) : N{N} {
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
            
            Particle* p = new Particle{x, y, vx, vy, r, m};

            int count = 0;
            for (int j = 0; j < particles.size(); j++) {
                if (particles[j]->dist_squared_to(p) > pow(particles[j]->r + r, 2)) {
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

Particles::Particles() {
    N = 2;
    Particle* p1 = new Particle({0.25, 0.50, 0.2, 0.0, 0.1, 1.0});
    Particle* p2 = new Particle({0.75, 0.50, -0.5, -0.0, 0.1, 1.0});
    particles.push_back(p1);
    particles.push_back(p2);
}

void Particles::initialize() {
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double t = particles[i]->get_collision_time_with(particles[j]);
            if (t != -1.0) {
                Collision c{t, particles[i], particles[j], collision_type_map.at("particle-particle"), 0, 0}; 
                pq.push(c);
            }
        }
    }

    for (int i = 0; i < N; i++) {
        double t_h = particles[i]->get_collision_time_with_horizontal_wall();
        if (t_h != -1.0) {
            Collision c{t_h, particles[i], NULL, collision_type_map.at("particle-horizontal_wall"), 0, 0}; 
            pq.push(c);
        }
            
        double t_v = particles[i]->get_collision_time_with_vertical_wall();
        if (t_v != -1.0) {
            Collision c{t_h, particles[i], NULL, collision_type_map.at("particle-vertical_wall"), 0, 0}; 
            pq.push(c);
        }
    }
}

void Particles::get_collisions(Particle* p) {
    for (int i = 0; i < N; i++) {
        if (particles[i] != p) {
            double t = p->get_collision_time_with(particles[i]);
            if (t != -1.0) {
                Collision c{t, p, particles[i], collision_type_map.at("particle-particle"), p->count, particles[i]->count}; 
                pq.push(c);
            }
        }
    }
    double t_h = p->get_collision_time_with_horizontal_wall();
    if (t_h != -1.0) {
        Collision c{t_h, p, NULL, collision_type_map.at("particle-horizontal_wall"), p->count, 0}; 
        pq.push(c);
    }
        
    double t_v = p->get_collision_time_with_vertical_wall();
    if (t_v != -1.0) {
        Collision c{t_h, p, NULL, collision_type_map.at("particle-vertical_wall"), p->count, 0}; 
        pq.push(c);
    }
}

void Particles::simulate() {
    // initialize step
    write_to_file("../data/0.txt");
    initialize();
    int i = 1;
    while (!pq.empty()) {
        move_forward(pq.top().t);
        if (pq.top().collision_type == collision_type_map.at("particle-particle")) {
            if (pq.top().p1_count == pq.top().p1->count && pq.top().p2_count == pq.top().p2->count) {
                pq.top().p1->resolve_collision_with(pq.top().p2);
                pq.top().p1->count += 1;
                pq.top().p2->count += 1;

                // get new collision
                get_collisions(pq.top().p1);
                get_collisions(pq.top().p2);
            }
        }
        else if (pq.top().collision_type == collision_type_map.at("particle-vertical_wall")) {
            if (pq.top().p1_count == pq.top().p1->count) {
                pq.top().p1->resolve_collision_with_vertical_wall();
                pq.top().p1->count += 1;

                // get new collision
                get_collisions(pq.top().p1);
            }
        }
        else {
            if (pq.top().p1_count == pq.top().p1->count) {
                pq.top().p1->resolve_collision_with_horizontal_wall();
                pq.top().p1->count += 1;

                // get new collision
                get_collisions(pq.top().p1);
            }
        }
        pq.pop();
        write_to_file("../data/" + std::to_string(i) + ".txt");
        i++;
        if (i > 10) {
            break;
        }
    }
}

void Particles::move_forward(double delta_t) {
    for (int i = 0; i < N; i++) {
        particles[i]->move_forward(delta_t);
    }
}

void Particles::write_to_file(std::string filename) {
    std::ofstream out; 
    out.open(filename);
    for (int i = 0; i < N; i++) {
        out << particles[i]->pos.x << " " << particles[i]->pos.y << " " << particles[i]->vel.x << " " << particles[i]->vel.y << " " << particles[i]->r << " " << particles[i]->m << "\n";
    }
    out.close();
}

void Particles::assert_no_overlap() {
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double d = particles[i]->dist_squared_to(particles[j]);
            assert(d > pow(particles[i]->r + particles[j]->r, 2));
        }
    }
}

Particles::~Particles() {
    for (int i = 0; i < N; i++) {
        delete particles[i];
    }
    particles.clear();
}
