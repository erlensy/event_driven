#pragma once
#include <random>
#include <iostream>

double get_rand(double min, double max);

struct vec {
    double x;
    double y;
};

vec operator+(vec& v1, vec& v2);

vec operator-(vec& v1, vec& v2);

double operator*(vec& v1, vec& v2);

vec operator*(vec& v, double c);

vec operator*(double c, vec& v);
