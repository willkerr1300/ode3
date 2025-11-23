///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg]
  double d;   // diameter of ball [m]
  double b;   // b coeff for air resistance
  double c;   // c coeff for air resistance
};

// State vector y:
// y[0]: x
// y[1]: vx
// y[2]: y
// y[3]: vy
// y[4]: z
// y[5]: vz

// Derivatives:
// f0 = vx
// f1 = ax
// f2 = vy
// f3 = ay
// f4 = vz
// f5 = az

double f_x(double t, const vector<double> &y, void *params) {
    return y[1];
}

double f_vx(double t, const vector<double> &y, void *params) {
    Params *p = (Params*)params;
    double v = sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]);
    if (v == 0) return 0;
    double Fd = p->b * v + p->c * v * v;
    return -Fd * (y[1]/v) / p->m;
}

double f_y(double t, const vector<double> &y, void *params) {
    return y[3];
}

double f_vy(double t, const vector<double> &y, void *params) {
    Params *p = (Params*)params;
    double v = sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]);
    if (v == 0) return 0;
    double Fd = p->b * v + p->c * v * v;
    return -Fd * (y[3]/v) / p->m;
}

double f_z(double t, const vector<double> &y, void *params) {
    return y[5];
}

double f_vz(double t, const vector<double> &y, void *params) {
    Params *p = (Params*)params;
    double v = sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]);
    if (v == 0) return -p->g;
    double Fd = p->b * v + p->c * v * v;
    return -Fd * (y[5]/v) / p->m - p->g;
}

double f_stop(double t, const vector<double> &y, void *params) {
    if (y[0] >= 18.5) return 1;
    return 0;
}

int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d=0.075;
  
  pars.b = 1.6e-4 * pars.d;
  pars.c = 0.25 * pars.d * pars.d;
  
  void *p_par = (void*) &pars;

  double xend=18.5;       // meters to plate
  double z0=1.4;             // height of release [m]
  double theta0=1;         // angle of velocity at release (degrees)
  double target_z = 0.9;
  
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // Setup ODE system
  vector<pfunc_t> v_fun(6);
  v_fun[0] = f_x;
  v_fun[1] = f_vx;
  v_fun[2] = f_y;
  v_fun[3] = f_vy;
  v_fun[4] = f_z;
  v_fun[5] = f_vz;

  double vPitch = 0;   // m/s of pitch needed to land in strike zone at 0.9 meters
  
  // Root finding loop (Secant Method)
  double v1 = 30.0; // Guess 1
  double v2 = 60.0; // Guess 2
  double tol = 1e-4;
  int max_iter = 20;
  
  double z_err1 = 0;
  double z_err2 = 0;
  
  // Helper lambda to run simulation and get z error
  auto get_z_error = [&](double v_guess) -> double {
      vector<double> y0(6);
      double rad = theta0 * M_PI / 180.0;
      y0[0] = 0;
      y0[1] = v_guess * cos(rad);
      y0[2] = 0;
      y0[3] = 0;
      y0[4] = z0;
      y0[5] = v_guess * sin(rad);
      
      // Run solver - FIXED: use t0 variable instead of literal 0
      double t0 = 0.0;
      auto tg = RK4SolveNA(v_fun, y0, 100, t0, 2.0, p_par, f_stop, 1e-6);
      
      // Get final state
      int n = tg[0].GetN();
      double t_final, x_final, z_final;
      tg[0].GetPoint(n-1, t_final, x_final);
      tg[4].GetPoint(n-1, t_final, z_final);
      
      
      return z_final - target_z;
  };
  
  z_err1 = get_z_error(v1);
  z_err2 = get_z_error(v2);
  
  for(int i=0; i<max_iter; ++i) {
      if(abs(z_err2) < tol) {
          vPitch = v2;
          break;
      }
      
      // Secant step
      double v_new = v2 - z_err2 * (v2 - v1) / (z_err2 - z_err1);
      
      v1 = v2;
      z_err1 = z_err2;
      v2 = v_new;
      z_err2 = get_z_error(v2);
      
      vPitch = v2; // Update current best
  }

  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    // Re-run with final vPitch to show plot
    vector<double> y0(6);
    double rad = theta0 * M_PI / 180.0;
    y0[0] = 0;
    y0[1] = vPitch * cos(rad);
    y0[2] = 0;
    y0[3] = 0;
    y0[4] = z0;
    y0[5] = vPitch * sin(rad);
    
    // FIXED: use t0 variable instead of literal 0
    double t0 = 0.0;
    auto tg = RK4SolveNA(v_fun, y0, 100, t0, 2.0, p_par, f_stop, 1e-6);
    
    TCanvas *c1 = new TCanvas("c1","Baseball Trajectory", 800, 600);
    c1->Divide(1,2);
    c1->cd(1);
    tg[4].SetTitle("z vs t");
    tg[4].Draw("AL");
    
    c1->cd(2);
    // Plot z vs x
    TGraph *gr_zx = new TGraph(tg[0].GetN(), tg[0].GetY(), tg[4].GetY());
    gr_zx->SetTitle("z vs x;x [m];z [m]");
    gr_zx->Draw("AL");
    
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}