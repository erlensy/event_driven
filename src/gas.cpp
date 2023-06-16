#include "gas.h"

bool operator<(const Collision& c1, const Collision& c2) {
    // compare operator needed in priority queue,
    // made such that pq.top() refers to the earliest time
    if (c1.t > c2.t) {
        return true;
    }
    return false;
}

Gas::Gas(int N, double epsilon, std::vector<double>& m, std::vector<double>& r, std::vector<double>& v0,
         double x_min, double x_max, double y_min, double y_max) : N{N}, epsilon{epsilon} {
    // make n particles with random posistions
    for (int i = 0; i < N; i++) {
        double theta = get_rand(0.0, 2.0 * M_PI);
        double vx = v0[i] * cos(theta);
        double vy = v0[i] * sin(theta);
        bool overlap = true;
        
        // while overlap draw new positions
        int failed_count = 0;
        do {
            double x = get_rand(x_min + r[i], x_max - r[i]);
            double y = get_rand(y_min + r[i], y_max - r[i]);
    
            // create particle
            Particle* p = new Particle{x, y, vx, vy, r[i], m[i]};
            
            // check if overlap
            int count = 0;
            for (int j = 0; j < particles.size(); j++) {
                if (particles[j]->dist_squared_to(p) > pow(particles[j]->r + r[i], 2)) {
                    count += 1;
                }
                else {
                    break;
                }
            }
            
            // if no overlap, add particle to particles
            if (count == particles.size()) {
                overlap = false;
                particles.push_back(p);
            }
            else {
                delete p;
                p = nullptr;
            }
            
            if (failed_count > 1000000) {
                throw std::runtime_error("Could not make particles without overlap. Try fewer particles or less radius");
            }
            failed_count++;
        }
        while (overlap);
    }
}

void Gas::make_particle(double x, double y, double vx, double vy, double r, double m) {
    // make one particle and add it to particles
    Particle* p = new Particle{x, y, vx, vy, r, m};
    particles.push_back(p);
    N += 1;
}

void Gas::simulate(int frames, double timestep, double end) {
    assert(particles.size() != 0);
    std::ofstream out = initialize_ofstream("../data/data.txt", timestep, frames);
    
    // time variables
    double t = 0.0;
    double rest = 0.0;
     
    // add all collisions to priority queue
    get_all_collisions(t);
    double ct = 0.0;
    
    // run simulation until end time
    while (t < end) {
        // calculate availabe time until next collision
        double dt = pq.top().t - t;
        int collision_type = pq.top().collision_type;
        auto p1 = pq.top().p1;
        auto p2 = pq.top().p2;
        int p1_count = pq.top().p1_count;
        int p2_count = pq.top().p2_count;
        
        // move particles forward and save at a multiplum of timesteps
        rest = move_forward_and_save_at_correct_times(out, t, dt, timestep, rest);
        t += dt;
        
        // resolve and calcualte new collisions
        manage_collision(collision_type, p1, p2, t);

        // check if next collision is valid
        do {
            // if not valid
            pq.pop();
            collision_type = pq.top().collision_type;
            p1 = pq.top().p1;
            p2 = pq.top().p2;
            p1_count = pq.top().p1_count;
            p2_count = pq.top().p2_count;
        }
        while (!(valid_collision(collision_type, p1_count, p1, p2_count, p2)));
    }
    out.close();
}

double Gas::move_forward_and_save_at_correct_times(std::ofstream& out, double t, double dt, double timestep, double rest) {
    if (rest < dt) {
        // move forward what we missed with the last bit,
        // such that we get outputs which are multiplum of timesteps
        move_forward(rest);
        t += rest;
        dt -= rest;
        
        // move particles forward multiplum of timesteps
        int n = dt / timestep;
        for (int i = 0; i < n; i++) {
            save_particles(out, t);
            move_forward(timestep);
            t += timestep;
        }
        save_particles(out, t);
        
        // move particles forward exactly to the next collision
        move_forward(dt - n * timestep);
        return rest;
    }
    // if (rest > dt) we just move to next collision without saving
    move_forward(dt);
    return rest - dt;
}

void Gas::move_forward(double dt) {
    // move all particles dt forward in time
    for (int i = 0; i < N; i++) {
        particles[i]->move_forward(dt);
    }
}

void Gas::manage_collision(int collision_type, Particle* p1, Particle* p2, double t) {
    // call correct function dependening on collision type
    
    if (collision_type == collision_type_map.at("particle-particle")) {
        manage_p_p_collision(p1, p2, t);
    }
    else if (collision_type == collision_type_map.at("particle-vertical_wall")) {
        manage_p_vw_collision(p1, t);
    }
    else {
        manage_p_hw_collision(p1, t);
    }
}

bool Gas::valid_collision(int collision_type, int p1_count, Particle* p1, int p2_count, Particle* p2) {
    // return true if c is a valid collision (collision count is equal
    // in collision c and for the particle)
    if (collision_type == collision_type_map.at("particle-particle")) {
        if (p1_count == p1->count && p2_count == p2->count) {
            return true;
        }
    }
    else if (p1_count == p1->count) {
        return true;
    }
    return false;
}

void Gas::manage_p_p_collision(Particle* p1, Particle* p2, double t) {
    // resolve and add new collisions to priority queue
    // only use for particle-particle collision
    p1->resolve_collision_with(p2, epsilon);

    p1->count += 1;
    p2->count += 1;
    
    get_collisions(p1, t);
    get_collisions(p2, t);
}

void Gas::manage_p_hw_collision(Particle* p, double t) {
    // resolve and add new collisions to priority queue
    // only use for particle-horizontal wall collision
    p->resolve_collision_with_hw(epsilon);
    p->count += 1;
    get_collisions(p, t);
}

void Gas::manage_p_vw_collision(Particle* p, double t) {
    // resolve and add new collisions to priority queue
    // only use for particle-vertical wall collision
    p->resolve_collision_with_vw(epsilon);
    p->count += 1;
    get_collisions(p, t);
}

void Gas::get_all_collisions(double t) {
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double dt = particles[i]->get_collision_time_with(particles[j]);
            // if dt is -1.0 collision time is infinite
            if (dt != -1.0) {
                Collision c{t + dt, particles[i], particles[j], collision_type_map.at("particle-particle"), particles[i]->count, particles[j]->count};
                pq.push(c);
            }
        }
    }
    for (int i = 0; i < N; i++) {
        get_collisions_walls(particles[i], t);
    }
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

void Gas::save_particles(std::ofstream& out, double& t) {
    for (int j = 0; j < N; j++) {
        out << t << " " << particles[j];
    }
    out << "\n";
}

std::ofstream Gas::initialize_ofstream(std::string filename, double timestep, int frames) {
    std::ofstream out;
    out.open(filename);
    out << "N timestep frames\n";
    out << "-----------------------------\n";
    out << N << " " << timestep << " " << frames << "\n\n";
    out << "x y vx vy r m\n";
    out << "-----------------------------\n";
    return out;
}

Gas::~Gas() {
    for (auto& p : particles) {
        delete p;
        p = nullptr;
    }
    particles.erase(std::remove(particles.begin(), particles.end(), nullptr), particles.end());
}
