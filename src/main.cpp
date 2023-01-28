#include "includes.h"
#include "gas.h"

void problem_1() {
    int N = 5000;
    std::vector<double> m(N, 1.0);
    std::vector<double> r(N, 0.001);
    std::vector<double> v0(N, 1.0);
    Gas g(N, 1.0, m, r, v0, 0.0, 1.0, 0.0, 1.0);
    std::cout << "Created particles\n";
    g.simulate(1000, 0.1, 1.0);
}

void problem_2() {
    int N = 5000;
    int m_0 = 1.0;
    std::vector<double> m(m_0, N);
    for (int i = 0; i < N / 2; i++) {
        m[i] = 4 * m_0;
    }
    std::vector<double> r(N, 0.001);
    std::vector<double> v0(N, 1.0);
    Gas g(N, 1.0, m, r, v0, 0.0, 1.0, 0.0, 1.0);
    g.simulate(1000, 0.1, 1.0);
}

void problem_3() {
    int N = 5000;
    int m_0 = 1.0;
    std::vector<double> m(N, m_0);
    for (int i = 0; i < N / 2; i++) {
        m[i] = 4 * m_0;
    }
    std::vector<double> r(N, 0.001);
    std::vector<double> v0(N, 1.0);
    Gas g(N, 1.0, m, r, v0, 0.0, 1.0, 0.0, 1.0);
    g.simulate(1000, 0.1, 1.0);
}

void problem_4() {
    int N = 999;
    int m_0 = 1.0;
    std::vector<double> m(N, m_0);
    std::vector<double> r(N, 0.005);
    std::vector<double> v0(N, 1.0);
    Gas g(N, 0.5, m, r, v0, 0.0, 1.0, 0.0, 0.5);
    g.make_particle(0.5, 0.75, 0, -1.0, 0.1, 10.0);
    std::cout << "Created particles\n";
    g.simulate(5000, 0.01, 1.0);
    std::cout << "Finished simulation\n";
}

void test_1() {
    int N = 1000;
    std::vector<double> m(N, 1.0);
    std::vector<double> r(N, 0.01);
    std::vector<double> v0(N, 1.0);
    Gas g(N, 1.0, m, r, v0, 0.0, 1.0, 0.0, 1.0);
    std::cout << "Created particles\n";
    double end = 1.0;
    g.simulate(1000, 0.001, 10.0);
}

int main() {
    test_1();
}
