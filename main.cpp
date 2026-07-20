#include "grads_io.hpp"
#include <cmath>
#include <iostream>
#include <numbers>
#include <ostream>

constexpr int NX = 83;
constexpr int NZ = 42;
constexpr double DX = 400.;
constexpr double DZ = 400.;
constexpr double DT = 2.;
constexpr double TIMEEND = 1200.;

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

inline double avg(double a, double b) { return 0.5 * (a + b); }

void setup_base_state() {
  const double PI_SFC = std::pow(P_SFC / P_0, R_D / C_PD);
  for (int k = 1; k < NZ - 1; k++) {
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

    if (k > 1)
      base.rhow[k] = avg(base.rhou[k], base.rhou[k - 1]);
    else if (k == 1)
      base.rhow[k] = P_SFC / (R_D * base.theta[k] * PI_SFC);
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
  base.rhow[0] = base.rhow[1];
  base.rhow[NZ - 1] = base.rhow[NZ - 2];
  base.pi[0] = base.pi[1];
  base.pi[NZ - 1] = base.pi[NZ - 2];
  base.p[0] = base.p[1];
  base.p[NZ - 1] = base.p[NZ - 2];
}

int main() {
  // setup base state
  base = {};
  past = {};
  now = {};
  future = {};
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

  for (int i = 0; i < NX; i++) {
    now.theta_p[i][0] = now.theta_p[i][1];
    now.theta_p[i][NZ - 1] = now.theta_p[i][NZ - 2];
  }
  // horizontal periodic boundary
  for (int k = 1; k < NZ - 1; k++) {
    now.theta_p[0][k] = now.theta_p[NX - 2][k];
    now.theta_p[NX - 1][k] = now.theta_p[1][k];
  }

  for (int i = 1; i < NX - 1; i++) {
    now.pi_p[i][NZ - 2] = 0.;
    now.pi_p[i][NZ - 1] = now.pi_p[i][NZ - 2];
    // integrate downward, assume pi'=0 at top
    for (int k = NZ - 3; k > 0; k--) {
      now.pi_p[i][k] =
          now.pi_p[i][k + 1] -
          G * now.theta_p[i][k] / (C_PD * std::pow(base.theta_v[k], 2.)) * DZ;
      // theta or theta_v??
    }
    now.pi_p[i][0] = now.pi_p[i][1];
  }

  for (int i = 1; i < NX - 1; i++) {
    for (int k = 1; k < NZ - 1; k++) {
      now.p_p[i][k] = now.pi_p[i][k] * C_PD * base.rhou[k] * base.theta_v[k];
    }
  }

  for (int i = 0; i < NX; i++) {
    now.w[i][NZ - 2] = 0;
    now.w[i][1] = 0;

    now.u[i][0] = now.u[i][1];
    now.u[i][NZ - 1] = now.u[i][NZ - 2];
    now.w[i][0] = now.w[i][1];
    now.w[i][NZ - 1] = now.w[i][NZ - 2];
    now.theta_p[i][0] = now.theta_p[i][1];
    now.theta_p[i][NZ - 1] = now.theta_p[i][NZ - 2];
    now.pi_p[i][0] = now.pi_p[i][1];
    now.pi_p[i][NZ - 1] = now.pi_p[i][NZ - 2];
    now.p_p[i][0] = now.p_p[i][1];
    now.p_p[i][NZ - 1] = now.p_p[i][NZ - 2];
  }

  // horizontal periodic boundary
  for (int k = 1; k < NZ - 1; k++) {
    now.u[0][k] = now.u[NX - 2][k];
    now.u[NX - 1][k] = now.u[1][k];
    now.w[0][k] = now.w[NX - 2][k];
    now.w[NX - 1][k] = now.w[1][k];
    now.theta_p[0][k] = now.theta_p[NX - 2][k];
    now.theta_p[NX - 1][k] = now.theta_p[1][k];
    now.pi_p[0][k] = now.pi_p[NX - 2][k];
    now.pi_p[NX - 1][k] = now.pi_p[1][k];
    now.p_p[0][k] = now.p_p[NX - 2][k];
    now.p_p[NX - 1][k] = now.p_p[1][k];
  }

  past = now;
  future = now;

  GradsWriter out("model", NX - 2, NZ - 2, 0., DX / 1000., 0., DZ / 1000.,
                  {"u", "w", "theta_p", "p_p"});

  double t = 0;
  while (t <= TIMEEND) {
    std::cout << "t=" << t << std::endl;
    out.put_field(now.u, 1, 1);
    out.put_field(now.w, 1, 1);
    out.put_field(now.theta_p, 1, 1);
    out.put_field(now.p_p, 1, 1);
    out.end_snapshot();
    double dt = t == 0 ? DT : 2. * DT;
    for (int i = 1; i < NX - 1; i++) {
      for (int k = 1; k < NZ - 1; k++) {
        future.u[i][k] =
            past.u[i][k] +
            dt * (-(std::pow(now.u[i + 1][k], 2.) -
                    std::pow(now.u[i - 1][k], 2.)) /
                      (2. * DX) -
                  (1 / base.rhou[k]) *
                      (base.rhow[k + 1] * avg(now.u[i][k + 1], now.u[i][k]) *
                           avg(now.w[i][k + 1], now.w[i - 1][k + 1]) -
                       base.rhow[k] * avg(now.u[i][k], now.u[i][k - 1]) *
                           avg(now.w[i][k], now.w[i - 1][k])) /
                      DZ -
                  C_PD * base.theta_v[k] *
                      (now.pi_p[i][k] - now.pi_p[i - 1][k]) / DX);
        future.w[i][k] =
            past.w[i][k] +
            dt * (-(avg(now.u[i + 1][k], now.u[i + 1][k - 1]) *
                        avg(now.w[i][k], now.w[i + 1][k]) -
                    avg(now.u[i][k], now.u[i][k - 1]) *
                        avg(now.w[i][k], now.w[i - 1][k])) /
                      DX -
                  (1 / base.rhow[k]) *
                      (base.rhow[k + 1] * std::pow(now.w[i][k + 1], 2.) -
                       base.rhow[k - 1] * std::pow(now.w[i][k - 1], 2.)) /
                      (2. * DZ) -
                  C_PD * avg(base.theta_v[k], base.theta_v[k - 1]) *
                      (now.pi_p[i][k] - now.pi_p[i][k - 1]) / DZ +
                  G * now.theta_p[i][k] / base.theta[k]);
        future.theta_p[i][k] =
            past.theta_p[i][k] +
            dt * (-(now.u[i + 1][k] *
                        avg(now.theta_p[i + 1][k], now.theta_p[i][k]) -
                    now.u[i][k] *
                        avg(now.theta_p[i][k], now.theta_p[i - 1][k])) /
                      DX -
                  (1 / base.rhou[k]) *
                      (base.rhow[k + 1] * now.w[i][k + 1] *
                           avg(now.theta_p[i][k + 1], now.theta_p[i][k]) -
                       base.rhow[k] * now.w[i][k] *
                           avg(now.theta_p[i][k], now.theta_p[i][k - 1])) /
                      DZ -
                  avg(now.w[i][k + 1], now.w[i][k]) *
                      (base.theta[k + 1] - base.theta[k - 1]) / (2 * DZ));
        double cs = 50.;
        future.pi_p[i][k] =
            past.pi_p[i][k] +
            dt * (-std::pow(cs, 2.) /
                  (base.rhou[k] * C_PD * std::pow(base.theta_v[k], 2.)) *
                  (base.rhou[k] * base.theta_v[k] *
                       (now.u[i + 1][k] - now.u[i][k]) / DX +
                   (base.rhow[k + 1] *
                        avg(base.theta_v[k + 1], base.theta_v[k]) *
                        now.w[i][k + 1] -
                    base.rhow[k] * avg(base.theta_v[k], base.theta_v[k - 1]) *
                        now.w[i][k]) /
                       DZ));

        future.p_p[i][k] =
            future.pi_p[i][k] * C_PD * base.rhou[k] * base.theta_v[k];
      }
    }

    // boundary conditon
    for (int i = 0; i < NX; i++) {
      future.w[i][NZ - 2] = 0;
      future.w[i][1] = 0;
      now.w[i][NZ - 2] = 0;
      now.w[i][1] = 0;

      future.u[i][0] = future.u[i][1];
      future.u[i][NZ - 1] = future.u[i][NZ - 2];
      future.w[i][0] = 0.;
      future.w[i][NZ - 1] = 0;
      future.theta_p[i][0] = future.theta_p[i][1];
      future.theta_p[i][NZ - 1] = future.theta_p[i][NZ - 2];
      future.pi_p[i][0] = future.pi_p[i][1];
      future.pi_p[i][NZ - 1] = future.pi_p[i][NZ - 2];
      future.p_p[i][0] = future.p_p[i][1];
      future.p_p[i][NZ - 1] = future.p_p[i][NZ - 2];
    }

    // horizontal periodic boundary
    for (int k = 1; k < NZ - 1; k++) {
      future.u[0][k] = future.u[NX - 2][k];
      future.u[NX - 1][k] = future.u[1][k];
      future.w[0][k] = future.w[NX - 2][k];
      future.w[NX - 1][k] = future.w[1][k];
      future.theta_p[0][k] = future.theta_p[NX - 2][k];
      future.theta_p[NX - 1][k] = future.theta_p[1][k];
      future.pi_p[0][k] = future.pi_p[NX - 2][k];
      future.pi_p[NX - 1][k] = future.pi_p[1][k];
      future.p_p[0][k] = future.p_p[NX - 2][k];
      future.p_p[NX - 1][k] = future.p_p[1][k];
    }

    past = now;
    now = future;
    t += DT;
  }

  return 0;
}