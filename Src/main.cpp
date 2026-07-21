#include "main.hpp"
#include "grads_io.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <ostream>

Grid past, now, future;
Base base;

inline double avg(double a, double b) { return 0.5 * (a + b); }
inline double qvs(double p, double t) {
  return 380.0 / p * std::exp(17.27 * (t - 273.) / (t - 36.));
}

Matrix total(const Matrix &_p, const Vector &base) {
  Matrix total;
  for (int i = 1; i < NX - 1; i++) {
    for (int k = 1; k < NZ - 1; k++) {
      total[i][k] = _p[i][k] + base[k];
    }
  }
  return total;
}

void setup_base_state() {
  const double PI_SFC = std::pow(P_SFC / P_0, R_D / C_PD);
  for (int k = 1; k < NZ - 1; k++) {
    base.theta[k] = 300.;

    // base.qv[k] = 0.;
    if (k >= 1) {
      double z = (k - 1) * DZ + DZ * 0.5;
      if (z <= 4000.)
        base.qv[k] = 0.0161 - 0.000003375 * z;
      else if (z <= 8000.)
        base.qv[k] = 0.0026 - 0.00000065 * (z - 4000.);
      else
        base.qv[k] = 0;
    }

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

    base.u[k] = 0;
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
  base.u[0] = base.u[1];
  base.u[NZ - 1] = base.u[NZ - 2];
}

void apply_kessler() {}

void apply_bc(Grid &g) {

  // boundary conditon
  for (int i = 0; i < NX; i++) {
    g.w[i][NZ - 2] = 0;
    g.w[i][1] = 0;

    g.u[i][0] = g.u[i][1];
    g.u[i][NZ - 1] = g.u[i][NZ - 2];
    g.w[i][0] = 0.;
    g.w[i][NZ - 1] = 0;
    g.theta_p[i][0] = g.theta_p[i][1];
    g.theta_p[i][NZ - 1] = g.theta_p[i][NZ - 2];
    g.pi_p[i][0] = g.pi_p[i][1];
    g.pi_p[i][NZ - 1] = g.pi_p[i][NZ - 2];
    g.p_p[i][0] = g.p_p[i][1];
    g.p_p[i][NZ - 1] = g.p_p[i][NZ - 2];
    g.qv_p[i][0] = g.qv_p[i][1];
    g.qv_p[i][NZ - 1] = g.qv_p[i][NZ - 2];
  }

  // horizontal periodic boundary
  for (int k = 1; k < NZ - 1; k++) {
    g.u[0][k] = g.u[NX - 2][k];
    g.u[NX - 1][k] = g.u[1][k];
    g.w[0][k] = g.w[NX - 2][k];
    g.w[NX - 1][k] = g.w[1][k];
    g.theta_p[0][k] = g.theta_p[NX - 2][k];
    g.theta_p[NX - 1][k] = g.theta_p[1][k];
    g.pi_p[0][k] = g.pi_p[NX - 2][k];
    g.pi_p[NX - 1][k] = g.pi_p[1][k];
    g.p_p[0][k] = g.p_p[NX - 2][k];
    g.p_p[NX - 1][k] = g.p_p[1][k];
    g.qv_p[0][k] = g.qv_p[NX - 2][k];
    g.qv_p[NX - 1][k] = g.qv_p[1][k];
  }
}

void apply_temporal_filter() {
  for (int i = 1; i < NX - 1; i++) {
    for (int k = 1; k < NZ - 1; k++) {
      now.u[i][k] +=
          0.5 * EPSILON * (future.u[i][k] - 2 * now.u[i][k] + past.u[i][k]);
      now.w[i][k] +=
          0.5 * EPSILON * (future.w[i][k] - 2 * now.w[i][k] + past.w[i][k]);
      now.theta_p[i][k] +=
          0.5 * EPSILON *
          (future.theta_p[i][k] - 2 * now.theta_p[i][k] + past.theta_p[i][k]);
      now.pi_p[i][k] +=
          0.5 * EPSILON *
          (future.pi_p[i][k] - 2 * now.pi_p[i][k] + past.pi_p[i][k]);
      now.qv_p[i][k] +=
          0.5 * EPSILON *
          (future.qv_p[i][k] - 2 * now.qv_p[i][k] + past.qv_p[i][k]);
    }
  }
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
      now.u[i][k] = base.u[k];
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
      now.qv_p[i][k] = (now.theta_p[i][k] != 0)
                           ? qvs(now.p_p[i][k] + base.p[k],
                                 now.theta_p[i][k] + base.theta[k]) -
                                 base.qv[k]
                           : 0.;
    }
  }

  apply_bc(now);

  past = now;
  future = now;

  GradsWriter out("model", NX - 2, NZ - 2, 0., DX / 1000., 0., DZ / 1000.,
                  {"u", "w", "theta_p", "theta", "p_p", "qv_p", "qv"});

  double t = 0;
  while (t <= TIMEEND) {
    std::cout << "t=" << t << std::endl;
    out.put_field(now.u, 1, 1);
    out.put_field(now.w, 1, 1);
    out.put_field(now.theta_p, 1, 1);
    out.put_field(total(now.theta_p, base.theta), 1, 1);
    out.put_field(now.p_p, 1, 1);
    out.put_field(now.qv_p, 1, 1);
    out.put_field(total(now.qv_p, base.qv), 1, 1);
    out.end_snapshot();
    double dt = t == 0 ? DT : 2. * DT;
    for (int i = 1; i < NX - 1; i++) {
      for (int k = 1; k < NZ - 1; k++) {
        future.u[i][k] =
            past.u[i][k] +
            dt *
                (-(std::pow(now.u[i + 1][k], 2.) -
                   std::pow(now.u[i - 1][k], 2.)) /
                     (2. * DX) -
                 (1 / base.rhou[k]) *
                     (base.rhow[k + 1] * avg(now.u[i][k + 1], now.u[i][k]) *
                          avg(now.w[i][k + 1], now.w[i - 1][k + 1]) -
                      base.rhow[k] * avg(now.u[i][k], now.u[i][k - 1]) *
                          avg(now.w[i][k], now.w[i - 1][k])) /
                     DZ -
                 C_PD * base.theta_v[k] *
                     (now.pi_p[i][k] - now.pi_p[i - 1][k]) / DX
#ifdef DIFFUSION
                 +
                 KX * (past.u[i + 1][k] - 2 * past.u[i][k] + past.u[i - 1][k]) /
                     (std::pow(DX, 2.)) +
                 KZ *
                     ((past.u[i][k + 1] - base.u[k + 1]) -
                      2 * (past.u[i][k] - base.u[k]) +
                      (past.u[i][k - 1] - base.u[k - 1])) /
                     (std::pow(DZ, 2.))
#endif
                );
        future.w[i][k] =
            past.w[i][k] +
            dt *
                (-(avg(now.u[i + 1][k], now.u[i + 1][k - 1]) *
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
                 G * (now.theta_p[i][k] / base.theta[k] + .61 * now.qv_p[i][k] -
                      (now.qr[i][k] + now.qc[i][k]))
#ifdef DIFFUSION
                 +
                 KX * (past.w[i + 1][k] - 2 * past.w[i][k] + past.w[i - 1][k]) /
                     (std::pow(DX, 2.)) +
                 KZ * (past.w[i][k + 1] - 2 * past.w[i][k] + past.w[i][k - 1]) /
                     (std::pow(DZ, 2.))
#endif
                );
        future.theta_p[i][k] =
            past.theta_p[i][k] +
            dt *
                (-(now.u[i + 1][k] *
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
                     (base.theta[k + 1] - base.theta[k - 1]) / (2 * DZ)
#ifdef DIFFUSION
                 +
                 KX *
                     (past.theta_p[i + 1][k] - 2 * past.theta_p[i][k] +
                      past.theta_p[i - 1][k]) /
                     (std::pow(DX, 2.)) +
                 KZ *
                     (past.theta_p[i][k + 1] - 2 * past.theta_p[i][k] +
                      past.theta_p[i][k - 1]) /
                     (std::pow(DZ, 2.))
#endif
                );
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
                        base.rhow[k] *
                            avg(base.theta_v[k], base.theta_v[k - 1]) *
                            now.w[i][k]) /
                           DZ)
#ifdef DIFFUSION
                  + KX *
                        (past.pi_p[i + 1][k] - 2 * past.pi_p[i][k] +
                         past.pi_p[i - 1][k]) /
                        (std::pow(DX, 2.)) +
                  KZ *
                      (past.pi_p[i][k + 1] - 2 * past.pi_p[i][k] +
                       past.pi_p[i][k - 1]) /
                      (std::pow(DZ, 2.))
#endif
                 );

        future.qv_p[i][k] =
            past.qv_p[i][k] +
            dt *
                (-(now.u[i + 1][k] * avg(now.qv_p[i][k], now.qv_p[i + 1][k]) -
                   now.u[i][k] * avg(now.qv_p[i][k], now.qv_p[i - 1][k])) /
                     DX -
                 (1 / base.rhou[k]) *
                     (base.rhow[k + 1] * now.w[i][k + 1] *
                          avg(now.qv_p[i][k + 1], now.qv_p[i][k]) -
                      base.rhow[k] * now.w[i][k] *
                          avg(now.qv_p[i][k], now.qv_p[i][k - 1])) /
                     DZ -
                 avg(now.w[i][k], now.w[i][k + 1]) *
                     (base.qv[k + 1] - base.qv[k - 1]) / (2. * DZ)
#ifdef DIFFUSION
                 +
                 KX *
                     (past.qv_p[i + 1][k] - 2 * past.qv_p[i][k] +
                      past.qv_p[i - 1][k]) /
                     (std::pow(DX, 2.)) +
                 KZ *
                     (past.qv_p[i][k + 1] - 2 * past.qv_p[i][k] +
                      past.qv_p[i][k - 1]) /
                     (std::pow(DZ, 2.))
#endif
                );
        if (future.qv_p[i][k] + base.qv[k] < 0)
          future.qv_p[i][k] = -base.qv[k];
      }
    }

    apply_kessler();
    apply_bc(future);
    apply_temporal_filter();
    apply_bc(now);

    for (int i = 1; i < NX - 1; i++)
      for (int k = 1; k < NZ - 1; k++)
        future.p_p[i][k] =
            future.pi_p[i][k] * C_PD * base.rhou[k] * base.theta_v[k];

    past = now;
    now = future;
    t += DT;
  }

  return 0;
}