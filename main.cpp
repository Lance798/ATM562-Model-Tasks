#include "grads_io.hpp"
#include <cmath>
#include <iostream>
#include <numbers>
#include <ostream>

constexpr int NX = 83;
constexpr int NZ = 42;
constexpr double DX = 400.;
constexpr double DZ = 400.;

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

// 2d grid
struct {
  // _p means prime/perturbation (deviate from base state)
  double u[NX][NZ], w[NX][NZ], theta_p[NX][NZ], pi_p[NX][NZ], p_p[NX][NZ];
} past, now, future;

// base state
struct {
  double theta[NZ], theta_v[NZ], qv[NZ], pi[NZ], p[NZ], rhou[NZ], rhow[NZ];
} base;

// prognostic

void setup_base_state() {
  const double PI_SFC = std::pow(P_SFC / P_0, R_D / C_PD);
  for (int k = 0; k < NZ; k++) {
    base.theta[k] = 300.; // TODO: ??

    base.qv[k] = 0.;
    // if (k >= 1) {
    //   double z = (k - 1) * DZ + DZ * 0.5;
    //   if (z <= 4000.)
    //     base.qv[k] = 0.0161 - 0.000003375 * z;
    //   else if (z <= 8000.)
    //     base.qv[k] = 0.0026 - 0.00000065 * (base.qv[k] - 4000.);
    //   else
    //     base.qv[k] = 0;
    // }

    base.theta_v[k] = base.theta[k] * (1 + .61 * base.qv[k]);
    if (k == 1) {
      base.pi[k] = PI_SFC - G / (C_PD * base.theta_v[k]) * 0.5 * DZ;
    } else if (k >= 2) {
      base.pi[k] =
          base.pi[k - 1] -
          G / (C_PD * ((base.theta_v[k] + base.theta_v[k - 1]) / 2.0)) * DZ;
    }
    base.p[k] = std::pow(base.pi[k], C_PD / R_D) * P_0;

    base.rhou[k] = base.p[k] / (R_D * base.theta_v[k] * base.pi[k]);
  }
  // boundary conditios
  base.theta[0] = base.theta[1];
  base.theta[NZ - 1] = base.theta[NZ - 2];
  base.theta_v[0] = base.theta_v[1];
  base.theta_v[NZ - 1] = base.theta_v[NZ - 2];
  base.qv[0] = base.qv[1];
  base.qv[NZ - 1] = base.qv[NZ - 2];
  base.rhou[0] = base.rhou[1];
  base.rhou[NZ - 1] = base.rhou[NZ - 2];
  base.pi[0] = base.pi[1];
  base.pi[NZ - 1] = base.pi[NZ - 2];
}

int main() {
  // setup base state
  std::cout << "START" << std::endl;

  setup_base_state();

  for (int i = 1; i < NX - 1; i++) {
    for (int k = 1; k < NZ - 1; k++) {
      double z_t = DZ * (k - 0.5);
      double rad = std::sqrt(std::pow((z_t - ZCNT) / RADZ, 2.) +
                             std::pow(DX * (i - IMID) / RADX, 2.));
      if (rad <= 1.)
        now.theta_p[i][k] =
            0.5 * TP_AMP * (std::cos(rad * std::numbers::pi) + 1);
      else
        now.theta_p[i][k] = 0.;
    }
  }

  for (int i = 1; i < NX - 1; i++) {
    now.pi_p[i][NZ - 2] = 0.;
    now.pi_p[i][NZ - 1] = now.pi_p[i][NZ - 2];
    // integrate downward, assume pi'=0 at top
    for (int k = NZ - 3; k > 0; k--) {
      double theta_p_top = now.theta_p[i][k + 1];
      double theta_p_btm = now.theta_p[i][k];
      double theta_p_avg = (theta_p_top + theta_p_btm) / 2;
      double theta_top = base.theta[k + 1];
      double theta_btm = base.theta[k];
      double theta_avg = (theta_top + theta_btm) / 2;

      now.pi_p[i][k] = now.pi_p[i][k + 1] -
                       G / C_PD * theta_p_avg / std::pow(theta_avg, 2.) * DZ;
      // why minus??
      // now.p_p
    }
    now.pi_p[i][0] = now.pi_p[i][1];
  }

  for (int i = 0; i < NX; i++) {
    for (int k = 0; k < NZ; k++) {
      future.theta_p[i][k] = now.theta_p[i][k];
      // TODO: p_p感覺有問題 在grid top?
      now.p_p[i][k] = now.pi_p[i][k] * C_PD * base.rhou[k] * base.theta_v[k];
    }
  }
  GradsWriter out("model", NX, NZ, 0., DX, 0., DZ / 1000, {"th", "w"});
  out.put_field(now.theta_p, 1, 1);
  out.put_field(now.p_p, 1, 1);
  out.put_field(now.pi_p, 1, 1);
  out.end_snapshot();
  return 0;
}