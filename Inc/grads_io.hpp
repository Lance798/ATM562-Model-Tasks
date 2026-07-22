// grads_io.hpp — zero-dependency model output: flat float32 .dat +
// auto-generated GrADS .ctl The SAME .dat feeds GrADS (via the .ctl) and Python
// (np.fromfile — see read_grads.py). C-ARRAY EDITION: put_field takes global
// arrays declared as  double f[NX][NZ]  indexed f[i][k].
//
// Record order (GrADS convention): x fastest, then z, then variable, then time.
// Usage once per output time:
//     out.put_field(th, 1, 1);     // same order as the names given to the
//     constructor out.put_field(w,  1, 1); out.end_snapshot();          //
//     appends one time level, rewrites .ctl (openable mid-run)

#ifndef __GRADS_IO_HPP__
#define __GRADS_IO_HPP__

#include "main.hpp"
#include <fstream>
#include <string>
#include <vector>

class GradsWriter {
public:
  // nx, nz      : number of REAL points to output (you exclude ghosts via i0,
  // k0) x0,dx,z0,dz : axis coordinates written into the .ctl (km gives tidy
  // GrADS axes)
  GradsWriter(const std::string &base, int nx, int nz, double x0, double dx,
              double z0, double dz, std::vector<std::string> var_names)
      : base_(base), nx_(nx), nz_(nz), x0_(x0), dx_(dx), z0_(z0), dz_(dz),
        vars_(std::move(var_names)),
        dat_(base + ".dat", std::ios::binary | std::ios::trunc) {}

  // f is a C array declared double f[FNX][FNZ], indexed f[i][k]. (i0,k0) =
  // first REAL point. Dimensions are deduced automatically, so arrays of any
  // size work.
  void put_field(const Matrix &f, int i0, int k0) {
    buf_.resize(static_cast<size_t>(nx_) * nz_);
    size_t n = 0;
    for (int k = 0; k < nz_; ++k)   // z outer
      for (int i = 0; i < nx_; ++i) // x fastest  (GrADS order)
        buf_[n++] = static_cast<float>(f[i0 + i][k0 + k]);
    dat_.write(reinterpret_cast<const char *>(buf_.data()),
               static_cast<std::streamsize>(buf_.size() * sizeof(float)));
  }

  void end_snapshot() {
    ++nt_;
    dat_.flush();
    write_ctl();
  }

private:
  void write_ctl() const {
    // DSET ^name.dat : ^ means "same directory as the .ctl"
    const std::string fname = base_.substr(base_.find_last_of('/') + 1);
    std::ofstream ctl(base_ + ".ctl", std::ios::trunc);
    ctl << "DSET ^" << fname << ".dat\n"
        << "TITLE 2D cloud model\n"
        << "UNDEF 1.e20\n"
        << "OPTIONS little_endian\n" // C++ stream write = no Fortran record
                                     // markers,
        << "XDEF " << nx_ << " LINEAR " << x0_ << ' ' << dx_
        << "\n" // so NO 'sequential'
        << "YDEF 1 LINEAR 0 1\n"
        << "ZDEF " << nz_ << " LINEAR " << z0_ << ' ' << dz_ << "\n"
        << "TDEF " << nt_ << " LINEAR 00Z01JAN2000 1mn\n" // TODO: not 1mn
        << "VARS " << vars_.size() << "\n";
    for (const auto &v : vars_)
      ctl << v << ' ' << nz_ << " 99 " << v << "\n";
    ctl << "ENDVARS\n";
  }

  std::string base_;
  int nx_, nz_, nt_ = 0;
  double x0_, dx_, z0_, dz_;
  std::vector<std::string> vars_;
  std::ofstream dat_;
  std::vector<float> buf_;
};

#endif