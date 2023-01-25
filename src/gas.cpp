#include "gas.h"

bool operator<(const Collision& c1, const Collision& c2) {
    // compare operator needed in priority queue,
    // made such that pq.top() refers to the earliest time
    if (c1.t > c2.t) {
        return true;
    }
    return false;
}

Gas::Gas(int N, double m, double r, double v0) : N{N} {
    std::vector<double> v0_(N, v0);
    std::vector<double> m_(N, m);
    std::vector<double> r_(N, r);
    make_particles_with_no_overlap(N, v0_, r_, m_, r, 1.0 - r, r, 0.5);
}

void Gas::make_particles_with_no_overlap(int n, std::vector<double> v0,
                                         std::vector<double> r, std::vector<double> m,
                                         double x_min, double x_max,
                                         double y_min, double y_max) {
    for (int i = 0; i < n; i++) {
        double theta = get_rand(0.0, 2.0 * M_PI);
        double vx = v0[i] * cos(theta);
        double vy = v0[i] * sin(theta);
        bool overlap = true;
        do {
            double x = get_rand(x_min, x_max);
            double y = get_rand(y_min, y_max);
            Particle* p = new Particle{x, y, vx, vy, r[i], m[i]};
            int count = 0;
            for (int j = 0; j < particles.size(); j++) {
                if (particles[j]->dist_squared_to(p) > pow(particles[j]->r + r[i], 2)) {
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

void Gas::simulate(int frames, double timestep) {
    double t = 0.0;
    for (int i = 0; i < N; i++) {
        get_collisions(particles[i], t);
    }

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
        double dt = pq.top().t - t;
        move_forward(dt);
        t = pq.top().t;
        for (int i = 0; i < N; i++) {
            out << particles[i];
        }
        out << "\n";
        frames_count++;

        // resolve and calcualte new collisions
        manage_collision(pq.top().collision_type, t);
        
        // check if next collision is valid
        do {
            pq.pop();
        }
        while (!(valid_collision(&pq.top())));
    }
    out.close();
}

void Gas::move_forward(double dt) {
    // move all particles dt forward in time
    for (int i = 0; i < N; i++) {
        particles[i]->move_forward(dt);
    }
}

void Gas::manage_collision(int collision_type, double t) {
    // call correct function dependening on collision type
    if (collision_type == collision_type_map.at("particle-particle")) {
        manage_p_p_collision(pq.top().p1, pq.top().p2, t);
    }
    else if (collision_type == collision_type_map.at("particle-vertical_wall")) {
        manage_p_vw_collision(pq.top().p1, t);
    }
    else {
        manage_p_hw_collision(pq.top().p1, t);
    }
}

bool Gas::valid_collision(const Collision* c) {
    // return true if c is a valid collision (collision count is equal
    // in collision c and for the particle)
    if (c->collision_type == collision_type_map.at("particle-particle")) {
        if (c->p1_count == c->p1->count && c->p2_count == c->p2->count) {
            return true;
        }
    }
    else if (c->p1_count == c->p1->count) {
        return true;
    }
    return false;
}

void Gas::manage_p_p_collision(Particle* p1, Particle* p2, double t) {
    // resolve and add new collisions to priority queue
    // only use for particle-particle collision
    pq.top().p1->resolve_collision_with(pq.top().p2);
    pq.top().p1->count += 1;
    pq.top().p2->count += 1;

    get_collisions(pq.top().p1, t);
    get_collisions(pq.top().p2, t);
}

void Gas::manage_p_hw_collision(Particle* p, double t) {
    // resolve and add new collisions to priority queue
    // only use for particle-horizontal wall collision
    p->resolve_collision_with_hw();
    p->count += 1;
    get_collisions(p, t);
}

void Gas::manage_p_vw_collision(Particle* p, double t) {
    // resolve and add new collisions to priority queue
    // only use for particle-vertical wall collision
    p->resolve_collision_with_vw();
    p->count += 1;
    get_collisions(p, t);
}

void Gas::get_collisions(Particle* p, double t) {
    get_collisions_pp(p, t);
    get_collisions_walls(p, t);
}

void Gas::get_collisions_pp(Particle* p, double t) {
    // add collisions between particle p and all other particles
    for (int i = 0; i < N; i++) {
        if (particles[i] != p) {
            double dt = p->get_collision_time_with(particles[i]);
            // if dt is -1.0 collision time is infinite
            if (dt != -1.0) {
                Collision c{t + dt, p, particles[i], collision_type_map.at("particle-particle"), p->count, particles[i]->count};
                pq.push(c);
            }
        }
    }
}

void Gas::get_collisions_walls(Particle* p, double t) {
    // add collisions between particle p and walls
    double dt_h = p->get_collision_time_with_hw();
    if (dt_h != -1.0) {
        Collision c{t + dt_h, p, NULL, collision_type_map.at("particle-horizontal_wall"), p->count, 0};
        pq.push(c);
    }
    double dt_v = p->get_collision_time_with_vw();
    if (dt_v != -1.0) {
        Collision c{t + dt_v, p, NULL, collision_type_map.at("particle-vertical_wall"), p->count, 0};
        pq.push(c);
    }
}

void Gas::assert_no_overlap() {
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double d = particles[i]->dist_squared_to(particles[j]);
            assert(d > pow(particles[i]->r + particles[j]->r, 2));
        }
    }
}

Gas::~Gas() {
    for (int i = 0; i < N; i++) {
        delete particles[i];
    }
    particles.clear();
}
