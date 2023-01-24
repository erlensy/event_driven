#include "particle.h"

Particle::Particle(double x, double y, double vx, double vy, double r, double m) : pos{x, y}, vel{vx, vy}, r{r}, m{m} {}

double Particle::dist_squared_to(Particle& particle) {
    vec delta_pos = pos - particle.pos;
    return delta_pos * delta_pos;
}
