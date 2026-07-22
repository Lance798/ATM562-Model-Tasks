#ifndef __MAIN_HPP__
#define __MAIN_HPP__

#include <array>

#define DIFFUSION

constexpr int NX = 1620;
constexpr int NZ = 420;
constexpr double DX = 40.;
constexpr double DZ = 40.;
constexpr double DT = 0.1;
constexpr double TIMEEND = 3600.;
// diffusion coefficient
constexpr double KX = 100.;
constexpr double KZ = 100.;
// Robert-Asselin filter
constexpr double EPSILON = 0.1;
constexpr double ALPHA = 0.53;

// symmetric thermal distribution
// theta perturbation amplititude
constexpr double TP_AMP = 3.; // K
constexpr double ZCNT = 3000; // m
constexpr double RADZ = 1000; // m
constexpr double RADX = 8000; // m
constexpr double IMID = (NX - 1) / 2. + 20;

constexpr double P_0 = 100000.0; // Pa
// pressure at model surface
constexpr double P_SFC = 96500.0; // Pa
// height of tropopause
constexpr double Z_TR = 12000; // m
// temperature at tropopause
constexpr double T_TR = 213; // K
// potential temperature at tropopause
constexpr double THETA_TR = 343; // K
// cloud water threshold
constexpr double QC_0 = 1e-3; // kg/kg
// conversion rate
constexpr double K_1 = 1e-3; // s^-1
// accretion rate
constexpr double K_2 = 2.2; // m^{7/8} kg^{7/8} s^-1

// ========== Physical Constants ==========
// heat capacity at constant pressure of dry air
constexpr double C_PD = 1004.0; // J kg^-1 K^-1
// universal gas constant of dry air
constexpr double R_D = 287.0;
// gravitational constant
constexpr double G = 9.81;
// latent heat of vaporization
constexpr double L_V = 2.5e6;

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