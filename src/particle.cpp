#include "particle.h"

Particle::Particle(double x, double y, double vx, double vy, double r, double m) : pos{x, y}, vel{vx, vy}, r{r}, m{m}, count{0} {}

double Particle::get_collision_time_with(Particle* p) {
    vec delta_vel = p->vel - vel;
    vec delta_pos = p->pos - pos;
    double c1 = delta_vel * delta_pos;

    // no collision
    if (c1 >= 0) {
        return - 1.0;
    }

    double c2 = delta_vel * delta_vel;
    double c3 = delta_pos * delta_pos;
    double R = r + p->r;
    double d = c1 * c1 - c2 * (c3 - R*R);

    // no collision
    if (d <= 0) {
        return - 1.0;
    }

    return -(c1 + sqrt(d)) / c2;
}

double Particle::get_collision_time_with_horizontal_wall() {
    if (vel.y == 0.0) {
        return - 1.0;
    }
    if (vel.y > 0) {
        return (1.0 - r - pos.y) / vel.y;
    }
    return (r - pos.y) / vel.y;
}

double Particle::get_collision_time_with_vertical_wall() {
    if (vel.x == 0.0) {
        return - 1.0;
    }
    if (vel.x > 0) {
        return (1.0 - r - pos.x) / vel.x;
    }
    return (r - pos.x) / vel.x;
}

void Particle::move_forward(double delta_t) {
    vec v = vel * delta_t;
    pos = pos + v;
}

void Particle::resolve_collision_with(Particle* p) {
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

void Particle::resolve_collision_with_horizontal_wall() {
    vel = vec{vel.x, -vel.y};
}

void Particle::resolve_collision_with_vertical_wall() {
    vel = vec{-vel.x, vel.y};
}

double Particle::dist_squared_to(Particle* p) {
    vec delta_pos = pos - p->pos;
    return delta_pos * delta_pos;
}

std::ostream& operator<<(std::ostream& os, Particle* p) {
    os << p->pos.x << " " << p->pos.y << " " << p->vel.x << " " << p->vel.y << " " << p->r << " " << p->m << "\n";
    return os;
}
