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

void Particles::simulate() {
    initialize();

    double t = 0.0;
    double timestep = 0.001;
    double end = timestep;
    int frames = 1000;

    std::ofstream out; 
    out.open("../data/skrt.txt");
    out << "N timestep frames\n";
    out << "-----------------------------\n";
    out << N << " " << timestep << " " << frames << "\n\n";
    out << "x y vx vy r m\n";
    out << "-----------------------------\n";
    for (int i = 0; i < N; i++) {
        out << particles[i];
    }
    out << "\n";
    int frames_count = 1;
    while (frames_count <= frames - 1) {
        // move forward
        double delta_t = pq.top().t - t;
        move_forward(delta_t);
        t = pq.top().t;
        for (int i = 0; i < N; i++) {
            out << particles[i];
        }
        out << "\n";
        frames_count++;

        // resolve and calcualte new collisions
        if (pq.top().collision_type == collision_type_map.at("particle-particle")) {
            p_p_collision(pq.top().p1, pq.top().p2, t);
        }
        else if (pq.top().collision_type == collision_type_map.at("particle-vertical_wall")) {
            p_vw_collision(pq.top().p1, t);
        }
        else {
            p_hw_collision(pq.top().p1, t);
        }
        
        // check if next collision is valid 
        while (true) {
            pq.pop();
            if (pq.top().collision_type == collision_type_map.at("particle-particle")) {
                if (pq.top().p1_count == pq.top().p1->count && pq.top().p2_count == pq.top().p2->count) {
                    break;
                }
            }
            else if (pq.top().p1_count == pq.top().p1->count) {
                break;
            }
        }
    }
    out.close();
}

void Particles::p_p_collision(Particle* p1, Particle* p2, double t) {
    pq.top().p1->resolve_collision_with(pq.top().p2);
    pq.top().p1->count += 1;
    pq.top().p2->count += 1;

    get_collisions(pq.top().p1, t);
    get_collisions(pq.top().p2, t);
}

void Particles::p_hw_collision(Particle* p, double t) {
    p->resolve_collision_with_horizontal_wall();
    p->count += 1;
    get_collisions(p, t);
}

void Particles::p_vw_collision(Particle* p, double t) {
    p->resolve_collision_with_vertical_wall();
    p->count += 1;
    get_collisions(p, t);
}

void Particles::get_collisions(Particle* p, double t) {
    get_collisions_pp(p, t);
    get_collisions_walls(p, t);
}

void Particles::get_collisions_pp(Particle* p, double t) {
    for (int i = 0; i < N; i++) {
        if (particles[i] != p) {
            double delta_t = p->get_collision_time_with(particles[i]);
            if (delta_t != -1.0) {
                Collision c{t + delta_t, p, particles[i], collision_type_map.at("particle-particle"), p->count, particles[i]->count}; 
                pq.push(c);
            }
        }
    }
}

void Particles::get_collisions_walls(Particle* p, double t) {
    double delta_t_h = p->get_collision_time_with_horizontal_wall();
    if (delta_t_h != -1.0) {
        Collision c{t + delta_t_h, p, NULL, collision_type_map.at("particle-horizontal_wall"), p->count, 0}; 
        pq.push(c);
    }
    double delta_t_v = p->get_collision_time_with_vertical_wall();
    if (delta_t_v != -1.0) {
        Collision c{t + delta_t_v, p, NULL, collision_type_map.at("particle-vertical_wall"), p->count, 0}; 
        pq.push(c);
    }
}

void Particles::initialize() {
    for (int i = 0; i < N; i++) {
        get_collisions(particles[i], 0);
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
        out << particles[i];
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
