#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include <streambuf>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>

#define N         32
// Thermal Diffusivity
#define V         1.0
#define END_STEP  30
// CFL Number
#define MU        1.0

using namespace std;
using namespace Eigen;

void displayField (VectorXd, unsigned int, double);
double smootherstep (unsigned int);
void fourier2DSquare (Ref<VectorXd>);

int main() {
  unsigned int end_step = END_STEP;
  double v        = V,
         h        = 1.0 / (N + 1),
         dt       = MU * h * h / v,
         c        = 0.5;

  // Initialize the temperature vector and fill with the square wave
  VectorXd u = VectorXd::Zero(N * N);  
  
  fourier2DSquare (u);
  
  double mu = v * dt / (h * h);
  if (mu > 0.5) 
    cerr << "Warning! Mu = " << mu << "; Solution may contain unphysical oscillations!" << endl;
  MatrixXd BBlock = MatrixXd::Zero (N, N);
  BBlock.diagonal (0) = VectorXd::Constant (N, -4 * mu);
  BBlock.diagonal (1) = BBlock.diagonal (-1) = VectorXd::Constant (N - 1, mu);
  MatrixXd A = MatrixXd::Zero(N * N, N * N);
  A.block<N, N>(0, 0) = BBlock;
  for (unsigned int i = 0; i < N - 1; ++i) {
    A.block<N, N>(N * (i + 1), N * (i + 1)) = BBlock;
    A.block<N, N>(N * (i + 1), N * i) = A.block<N, N>(N * i, N * (i + 1)) =
      MatrixXd::Identity(N, N) * mu;
  }

  displayField(u, 0, c);

  BiCGSTAB<MatrixXd> solver;

  for (unsigned int timestep = 1; timestep < end_step; ++timestep) {
    c = 0.5;
    VectorXd rhs = (MatrixXd::Identity(N * N, N * N) + c * A) * u;
    MatrixXd B   = (MatrixXd::Identity(N * N, N * N) - (1 - c) * A);
    solver.compute(B);
    u = solver.solve(rhs);
    cerr << "residual: " << (B * u - rhs).norm() << endl;

    displayField(u, timestep, c);
  }

  return 0;
}

void displayField (VectorXd u, unsigned int timestep, double c) {
  stringstream filename;
  filename << "diffusion-" << setw (3) << setfill ('0') << timestep << "-" << c << ".dat";
  ofstream out_file (filename.str().c_str());

  for (int i = 0; i < N; ++i)
    out_file << u.segment<N>(i * N).transpose() << endl;
  
  out_file.close ();
}

double smootherstep (const unsigned int timestep) {
  unsigned int lower = 2,
               upper = 6;
  if (timestep < lower) 
    return 0;
  if (timestep > upper)
    return 0.5;
  double t = (double(timestep) - lower) / (upper - lower) ;
  return t * t * t * (t * (t * 6 - 15) + 10) * 0.5;
}

void fourier2DSquare (Ref<VectorXd> u) {
  unsigned int n = N;
  double h = 1.0 / (N + 1),
         t = 1.25e-4;
  Vector2d x = Vector2d::Constant (0.5 * h);
  VectorXd bk (n);
  for (unsigned int k = 0; k < n; ++k)
    bk [k] = 2.0 * (cos ((k + 1) * M_PI / 4) - cos (3 * (k + 1) * M_PI / 4)) / ((k + 1) * M_PI);
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      u(i * n + j) = 0;
        for (int k1 = 0; k1 < int(n); ++k1) {
          for (int k2 = 0; k2 < int(n); ++k2) {
            u(i * n + j) += bk[k1] * bk[k2] * std::exp (-(k1 + 1) * (k1 + 1) * V * t * M_PI * M_PI) * std::exp (-(k2 + 1) * (k2 + 1) * V * t * M_PI * M_PI) * sin ((k1 + 1) * M_PI * x[0]) * sin ((k2 + 1) * M_PI * x[1]);
          }
        }
      x(0) += h;
    }
    x(0) = (0.5 * h); x(1) += h;
  }
}
