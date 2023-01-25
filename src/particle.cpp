#include "particle.h"

Particle::Particle(double x, double y, double vx, double vy, double r, double m) : pos{x, y}, vel{vx, vy}, r{r}, m{m}, count{0} {}

void Particle::move_forward(double dt) {
    // move *this particle dt forward in time
    vec v = vel * dt;
    pos = pos + v;
}

double Particle::get_collision_time_with(Particle* p) {
    // returns collision time between *this particle
    // and particle p. returns -1.0 if no collision is found
    
    vec delta_vel = p->vel - vel;
    vec delta_pos = p->pos - pos;
    double c1 = delta_vel * delta_pos;

    // infinite collision time = no collision
    if (c1 >= 0) {
        return - 1.0;
    }

    double c2 = delta_vel * delta_vel;
    double c3 = delta_pos * delta_pos;
    double R = r + p->r;
    double d = c1 * c1 - c2 * (c3 - R*R);

    // infinite collision time = no collision
    if (d <= 0) {
        return - 1.0;
    }
    
    return -(c1 + sqrt(d)) / c2;
}

void Particle::resolve_collision_with(Particle* p) {
    // updates velocities of *this particle and particle p
    // should only be called if these two particles
    // currently are colliding
    double RR = (r + p->r) * (r + p->r);
    double rho1 = p->m / (m + p->m);
    double rho2 = m / (m + p->m);
    vec delta_pos = p->pos - pos;
    vec delta_vel = p->vel - vel;
    double c1 = delta_pos * delta_vel;
    vec vel1 = 2.0 * rho1 * c1 / RR * delta_pos;
    vec vel2 = 2.0 * rho2 * c1 / RR * delta_pos;
    vel = vel + vel1;
    p->vel = p->vel - vel2;
}

double Particle::get_collision_time_with_hw() {
    // returns collision time with a horizontal wall
    // returns -1.0 if no collision is found
    if (vel.y == 0.0) {
        return - 1.0;
    }
    if (vel.y > 0) {
        return (1.0 - r - pos.y) / vel.y;
    }
    return (r - pos.y) / vel.y;
}

void Particle::resolve_collision_with_hw() {
    // updates velocity of *this particle under
    // a collision with a horizontal wall
    vel = vec{vel.x, -vel.y};
}

double Particle::get_collision_time_with_vw() {
    // returns collision time with a vertical wall
    // returns -1.0 if no collision is found
    if (vel.x == 0.0) {
        return - 1.0;
    }
    if (vel.x > 0) {
        return (1.0 - r - pos.x) / vel.x;
    }
    return (r - pos.x) / vel.x;
}

void Particle::resolve_collision_with_vw() {
    // updates velocity of *this particle under
    // a collision with a vertical wall
    vel = vec{-vel.x, vel.y};
}

double Particle::dist_squared_to(Particle* p) {
    // calculate distance squared from *this
    // particle to particle p
    vec delta_pos = pos - p->pos;
    return delta_pos * delta_pos;
}

std::ostream& operator<<(std::ostream& os, Particle* p) {
    // overload operator << to print particle information
    os << p->pos.x << " " << p->pos.y << " " <<
          p->vel.x << " " << p->vel.y << " " <<
          p->r << " " << p->m << "\n";
    return os;
}
