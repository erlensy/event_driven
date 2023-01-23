#include "particle.h"

Particle::Particle(double x, double y, double vx, double vy, double r, double m) : x{x}, y{y}, vx{vx}, vy{vy}, r{r}, m{m} {}

double Particle::dist_squared_to(Particle& particle) {
    return pow(x - particle.x, 2) + pow(y - particle.y, 2);
}
