#include "utils.h"

// "true" randomnes for seed generation 
std::random_device rd{};

// seed random number generator
std::mt19937 gen{rd()};

// uniform distribution between 0.0 and 1.0
std::uniform_real_distribution<double> uniform_distribution{0.0, 1.0};

double get_rand(double min, double max) {
    return min + (max - min) * uniform_distribution(gen);
}

vec operator+(vec& v1, vec& v2) {
    return vec{v1.x + v2.x, v1.y + v2.y};
}

vec operator-(vec& v1, vec& v2) {
    return vec{v1.x - v2.x, v1.y - v2.y};
}

double operator*(vec& v1, vec& v2) {
    return v1.x * v2.x + v1.y * v2.y;
} 
