#include "main.hpp"
#include "grads_io.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <ostream>

Grid past, now, future;
Base base;

inline double sqr(double a) { return a * a; }
inline double avg(double a, double b) { return 0.5 * (a + b); }
inline double calc_temp(const Grid &g, int i, int k) {
  return (base.theta[k] + g.theta_p[i][k]) * (base.pi[k] + g.pi_p[i][k]);
}
inline double calc_qvs(const Grid &g, int i, int k) {
  double t = calc_temp(g, i, k);
  double p = base.p[k] + g.p_p[i][k];
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
    double z = (k - 1) * DZ + DZ * 0.5;

    if (z <= Z_TR)
      base.theta[k] = 300. + 43. * std::pow(z / Z_TR, 1.25);
    else
      base.theta[k] = THETA_TR * std::exp(G * (z - Z_TR) / (C_PD * T_TR));

    // base.qv[k] = 0.;
    if (z <= 4000.)
      base.qv[k] = 0.0161 - 0.000003375 * z;
    else if (z <= 8000.)
      base.qv[k] = 0.0026 - 0.00000065 * (z - 4000.);
    else
      base.qv[k] = 0;

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

    if (z > 6000.)
      base.u[k] = 0.;
    else
      base.u[k] = -10. + z * (10. / 6000.);
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

void handle_kessler(Grid &g, double dt) {
  for (int i = 1; i < NX - 1; i++) {
    for (int k = 1; k < NZ - 1; k++) {
      // autoconversion of qc -> qr
      if (g.qc[i][k] > QC_0) {
        double a = std::min(g.qc[i][k] - QC_0, (g.qc[i][k] - QC_0) * K_1 * dt);
        g.qr[i][k] += a;
        g.qc[i][k] -= a;
      }

      // accretion of rain water
      double b = K_2 * g.qc[i][k] * std::pow(base.rhou[k] * g.qr[i][k], 7. / 8);
      if (g.qc[i][k] - b * dt < 0) {
        g.qr[i][k] += g.qc[i][k];
        g.qc[i][k] = 0.;
      } else {
        g.qc[i][k] -= b * dt;
        g.qr[i][k] += b * dt;
      }

      // evaporiation of rain water
      double qvs = calc_qvs(g, i, k);
      double c_vent = 1.6 + 30.39 * std::pow(base.rhou[k] * g.qr[i][k], 0.2046);
      double e = (1 / base.rhou[k]) *
                 ((1 - (base.qv[k] + g.qv_p[i][k]) / qvs) * c_vent *
                  std::pow(base.rhou[k] * g.qr[i][k], 0.525)) /
                 (2.03e4 + (9.58e6) / (base.p[k] * qvs));
      if (g.qr[i][k] - e * dt < 0) {
        g.qv_p[i][k] += g.qr[i][k];
        g.theta_p[i][k] -= L_V * g.qr[i][k] / (C_PD * base.pi[k]);
        g.qr[i][k] = 0.;
      } else {
        g.qr[i][k] -= e * dt;
        g.qv_p[i][k] += e * dt;
        g.theta_p[i][k] -= L_V * (e * dt) / (C_PD * base.pi[k]);
      }
    }
  }
}

void handler_water(Grid &g) {
  // check RH, if RH > 100, qv -> qc, elif RH<100 qc -> qv
  for (int i = 1; i < NX - 1; i++) {
    for (int k = 1; k < NZ - 1; k++) {
      double qv = g.qv_p[i][k] + base.qv[k];
      double t = calc_temp(g, i, k);
      double qvs = calc_qvs(g, i, k);
      double phi = qvs * (17.27 * 237. * L_V / (C_PD * sqr(t - 36.)));
      double c = (qv - qvs) / (1 + phi);
      // RH > 100%
      if (c > 0) {
        g.qv_p[i][k] -= c;
        g.qc[i][k] += c;
        g.theta_p[i][k] += L_V * c / (C_PD * base.pi[k]);
      } else if (c < 0 && g.qc[i][k] > 0) {
        double evap = std::min(-c, g.qc[i][k]);
        g.qv_p[i][k] += evap;
        g.qc[i][k] -= evap;
        g.theta_p[i][k] -= L_V * evap / (C_PD * base.pi[k]);
      }
    }
  }
}

void apply_open_bc(double dt) {
  const double cstar = 30.0;
  for (int k = 1; k < NZ - 1; k++) {
    future.u[1][k] = past.u[1][k] - std::min(now.u[1][k] - cstar, 0.0) *
                                        (past.u[2][k] - past.u[1][k]) *
                                        (dt / DX);
    future.u[NX - 1][k] =
        past.u[NX - 1][k] - std::max(now.u[NX - 1][k] + cstar, 0.0) *
                                (past.u[NX - 1][k] - past.u[NX - 2][k]) *
                                (dt / DX);
  }
}

void apply_bc(Grid &g) {

  // boundary conditon
  for (int i = 0; i < NX; i++) {
    // g.w[i][NZ - 2] = 0;
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
    g.qc[i][0] = g.qc[i][1];
    g.qc[i][NZ - 1] = g.qc[i][NZ - 2];
    g.qr[i][0] = g.qr[i][1];
    g.qr[i][NZ - 1] = g.qr[i][NZ - 2];
  }

  // lateral zero-gradient boundary
  for (int k = 1; k < NZ - 1; k++) {
    g.w[0][k] = g.w[1][k];
    g.w[NX - 1][k] = g.w[NX - 2][k];
    g.theta_p[0][k] = g.theta_p[1][k];
    g.theta_p[NX - 1][k] = g.theta_p[NX - 2][k];
    g.pi_p[0][k] = g.pi_p[1][k];
    g.pi_p[NX - 1][k] = g.pi_p[NX - 2][k];
    g.p_p[0][k] = g.p_p[1][k];
    g.p_p[NX - 1][k] = g.p_p[NX - 2][k];
    g.qv_p[0][k] = g.qv_p[1][k];
    g.qv_p[NX - 1][k] = g.qv_p[NX - 2][k];
    g.qc[0][k] = g.qc[1][k];
    g.qc[NX - 1][k] = g.qc[NX - 2][k];
    g.qr[0][k] = g.qr[1][k];
    g.qr[NX - 1][k] = g.qr[NX - 2][k];
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
      double rad =
          std::sqrt(sqr((z_t - ZCNT) / RADZ) + sqr(DX * (i - IMID) / RADX));
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
      now.pi_p[i][k] = now.pi_p[i][k + 1] - G * now.theta_p[i][k] /
                                                (C_PD * sqr(base.theta_v[k])) *
                                                DZ;
      // theta or theta_v??
    }
    now.pi_p[i][0] = now.pi_p[i][1];
  }

  for (int i = 1; i < NX - 1; i++) {
    for (int k = 1; k < NZ - 1; k++) {
      now.p_p[i][k] = now.pi_p[i][k] * C_PD * base.rhou[k] * base.theta_v[k];
      now.qv_p[i][k] =
          (now.theta_p[i][k] != 0) ? calc_qvs(now, i, k) - base.qv[k] : 0.;
    }
  }

  apply_bc(now);

  past = now;
  future = now;

  GradsWriter out(
      "model", NX - 2, NZ - 2, 0., DX / 1000., 0., DZ / 1000.,
      {"u", "w", "theta_p", "theta", "p_p", "qv_p", "qv", "qc", "qr", "pi"});

  double t = 0;
  while (t <= TIMEEND) {
    std::cout << "t=" << t << std::endl;
    if (static_cast<int>(t) % 50 == 0) {
      std::cout << "OUT" << std::endl;
      out.put_field(now.u, 1, 1);
      out.put_field(now.w, 1, 1);
      out.put_field(now.theta_p, 1, 1);
      out.put_field(total(now.theta_p, base.theta), 1, 1);
      out.put_field(now.p_p, 1, 1);
      out.put_field(now.qv_p, 1, 1);
      out.put_field(total(now.qv_p, base.qv), 1, 1);
      out.put_field(now.qc, 1, 1);
      out.put_field(now.qr, 1, 1);
      out.put_field(total(now.pi_p, base.pi), 1, 1);
      out.end_snapshot();
    }
    double dt = t == 0 ? DT : 2. * DT;
    for (int i = 1; i < NX - 1; i++) {
      for (int k = 1; k < NZ - 1; k++) {
        future.u[i][k] =
            past.u[i][k] +
            dt *
                (-(sqr(now.u[i + 1][k]) - sqr(now.u[i - 1][k])) / (2. * DX) -
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
                     (sqr(DX)) +
                 KZ *
                     ((past.u[i][k + 1] - base.u[k + 1]) -
                      2 * (past.u[i][k] - base.u[k]) +
                      (past.u[i][k - 1] - base.u[k - 1])) /
                     (sqr(DZ))
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
                     (base.rhow[k + 1] * sqr(now.w[i][k + 1]) -
                      base.rhow[k - 1] * sqr(now.w[i][k - 1])) /
                     (2. * DZ) -
                 C_PD * avg(base.theta_v[k], base.theta_v[k - 1]) *
                     (now.pi_p[i][k] - now.pi_p[i][k - 1]) / DZ +
                 G * (avg(now.theta_p[i][k], now.theta_p[i][k - 1]) /
                          avg(base.theta[k], base.theta[k - 1]) +
                      .61 * avg(now.qv_p[i][k], now.qv_p[i][k - 1]) -
                      (avg(now.qr[i][k], now.qr[i][k - 1]) +
                       avg(now.qc[i][k], now.qc[i][k - 1])))
#ifdef DIFFUSION
                 +
                 KX * (past.w[i + 1][k] - 2 * past.w[i][k] + past.w[i - 1][k]) /
                     (sqr(DX)) +
                 KZ * (past.w[i][k + 1] - 2 * past.w[i][k] + past.w[i][k - 1]) /
                     (sqr(DZ))
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
                     (sqr(DX)) +
                 KZ *
                     (past.theta_p[i][k + 1] - 2 * past.theta_p[i][k] +
                      past.theta_p[i][k - 1]) /
                     (sqr(DZ))
#endif
                );
        const double cs = 50.;
        future.pi_p[i][k] =
            past.pi_p[i][k] +
            dt * (-sqr(cs) / (base.rhou[k] * C_PD * sqr(base.theta_v[k])) *
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
                        (sqr(DX)) +
                  KZ *
                      (past.pi_p[i][k + 1] - 2 * past.pi_p[i][k] +
                       past.pi_p[i][k - 1]) /
                      (sqr(DZ))
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
                     (sqr(DX)) +
                 KZ *
                     (past.qv_p[i][k + 1] - 2 * past.qv_p[i][k] +
                      past.qv_p[i][k - 1]) /
                     (sqr(DZ))
#endif
                );

        future.qc[i][k] =
            past.qc[i][k] +
            dt * (-(now.u[i + 1][k] * avg(now.qc[i][k], now.qc[i + 1][k]) -
                    now.u[i][k] * avg(now.qc[i][k], now.qc[i - 1][k])) /
                      DX -
                  (1 / base.rhou[k]) *
                      (base.rhow[k + 1] * now.w[i][k + 1] *
                           avg(now.qc[i][k + 1], now.qc[i][k]) -
                       base.rhow[k] * now.w[i][k] *
                           avg(now.qc[i][k], now.qc[i][k - 1])) /
                      DZ
#ifdef DIFFUSION
                  + KX *
                        (past.qc[i + 1][k] - 2 * past.qc[i][k] +
                         past.qc[i - 1][k]) /
                        (sqr(DX)) +
                  KZ *
                      (past.qc[i][k + 1] - 2 * past.qc[i][k] +
                       past.qc[i][k - 1]) /
                      (sqr(DZ))
#endif
                 );
        // terminal velocity of raindrop
        double v_t = 6; // m/s
        future.qr[i][k] =
            past.qr[i][k] +
            dt * (-(now.u[i + 1][k] * avg(now.qr[i][k], now.qr[i + 1][k]) -
                    now.u[i][k] * avg(now.qr[i][k], now.qr[i - 1][k])) /
                      DX -
                  (1 / base.rhou[k]) *
                      (base.rhow[k + 1] * (now.w[i][k + 1] - v_t) *
                           avg(now.qr[i][k + 1], now.qr[i][k]) -
                       base.rhow[k] * (now.w[i][k] - v_t) *
                           avg(now.qr[i][k], now.qr[i][k - 1])) /
                      DZ
#ifdef DIFFUSION
                  + KX *
                        (past.qr[i + 1][k] - 2 * past.qr[i][k] +
                         past.qr[i - 1][k]) /
                        (sqr(DX)) +
                  KZ *
                      (past.qr[i][k + 1] - 2 * past.qr[i][k] +
                       past.qr[i][k - 1]) /
                      (sqr(DZ))
#endif
                 );

        // TODO: positive-definite
        if (future.qv_p[i][k] + base.qv[k] < 0)
          future.qv_p[i][k] = -base.qv[k];
        future.qc[i][k] = std::max(0., future.qc[i][k]);
        future.qr[i][k] = std::max(0., future.qr[i][k]);
      }
    }
    apply_open_bc(dt);
    apply_bc(future);
    handle_kessler(future, dt);
    handler_water(future);

    // for (int i = 1; i < NX - 1; i++) {
    //   for (int k = 1; k < NZ - 1; k++) {
    //     future.u[i][k] += base.u[k];
    //   }
    // }
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

  // TODO: analysis drift

  return 0;
}