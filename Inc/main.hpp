#ifndef __MAIN_HPP__
#define __MAIN_HPP__

#include <array>
#include <vector>

#define DIFFUSION

// TODO: I think there's not aligned between MODEL_NX and MODEL_NZ
//  may affect grads_wrtier
constexpr int MODEL_NX = 1620;
constexpr int MODEL_NZ = 420;
constexpr double DX = 40.;
constexpr double DZ = 40.;
constexpr double DT = 0.1;
constexpr double TIMEEND = 3600;
// write output data every OUT_FREQ sec
constexpr int OUT_PERIOD = 100;

// diffusion coefficient
constexpr double KX = 5.;
constexpr double KZ = 5.;
// Robert-Asselin filter
constexpr double EPSILON = 0.1;
constexpr double ALPHA = 0.53;

// symmetric thermal distribution
// theta perturbation amplititude
constexpr double TP_AMP = 3.; // K
constexpr double ZCNT = 3000; // m
constexpr double RADZ = 1000; // m
constexpr double RADX = 5000; // m
constexpr int IMID = (MODEL_NX - 1) / 2;

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

using Column = std::array<double, MODEL_NZ>;
using Matrix = std::vector<Column>;
namespace VarId {
enum VarName { U, W, THETA_P, PI_P, P_P, QV_P, QC, QR, NUM_VARS_ };
} // namespace VarId

// 2d grid
struct Grid {
  // _p means prime/perturbation (deviate from base state)
  // Matrix u, w, theta_p, pi_p, p_p, qv_p, qc, qr;
  std::array<Matrix, VarId::NUM_VARS_> vars;

  Matrix &u = vars[VarId::U];
  Matrix &w = vars[VarId::W];
  Matrix &theta_p = vars[VarId::THETA_P];
  Matrix &pi_p = vars[VarId::PI_P];
  Matrix &p_p = vars[VarId::P_P];
  Matrix &qv_p = vars[VarId::QV_P];
  Matrix &qc = vars[VarId::QC];
  Matrix &qr = vars[VarId::QR];

  void resize_nx(int n) {
    for (auto &var : vars) {
      var.resize(n);
    }
  }

  // 1. 預設建構子
  Grid() = default;

  // 2. 🚨 自訂「拷貝建構子」 (解決陷阱的關鍵！)
  // 這樣寫會只拷貝 vars 陣列，而下面的 reference 會自動綁定到自己全新的 vars
  // 身上！
  Grid(const Grid &other) : vars(other.vars) {}

  // 3. 🚨 自訂「賦值運算子」 (把你之前刪掉的加回來)
  Grid &operator=(const Grid &other) {
    vars = other.vars;
    return *this;
  }
};

// base state
struct Base {
  Column u, theta, theta_v, qv, pi, p, rhou, rhow;
};

extern Grid *past, *now, *future;
extern Base base;
extern unsigned grid_nx;

#endif