#ifndef __MAIN_HPP__
#define __MAIN_HPP__

#include <array>

#define DIFFUSION

constexpr int NX = 83;
constexpr int NZ = 42;
constexpr double DX = 400.;
constexpr double DZ = 400.;
constexpr double DT = 2.;
constexpr double TIMEEND = 4800.;
// diffusion coefficient
constexpr double KX = 100;
constexpr double KZ = 10;
// Robert-Asselin filter
constexpr double EPSILON = 0.1;
constexpr double ALPHA = 0.53;

// symmetric thermal distribution
// theta perturbation amplititude
constexpr double TP_AMP = 3.; // K
constexpr double ZCNT = 3000; // m
constexpr double RADZ = 4000; // m
constexpr double RADX = 4000; // m
constexpr double IMID = (NX - 1) / 2.;

constexpr double P_0 = 100000.0;
// pressure at model surface
constexpr double P_SFC = 96500.0;

// ========== Physical Constants ==========
// heat capacity at constant pressure of dry air
constexpr double C_PD = 1004.0; // J kg^-1 K^-1
// universal gas constant of dry air
constexpr double R_D = 287.0;
constexpr double G = 9.81;

using Matrix = std::array<std::array<double, NZ>, NX>;
using Vector = std::array<double, NZ>;

// 2d grid
struct Grid {
  // _p means prime/perturbation (deviate from base state)
  Matrix u, w, theta_p, pi_p, p_p, qv_p, qc, qr;
};

// base state
struct Base {
  Vector u, theta, theta_v, qv, pi, p, rhou, rhow;
};

#endif