// used by kepcart

#ifndef _KEPCART_
#define _KEPCART_

// phasespace and orbital element structures

typedef struct
{
  double x, y, z, xd, yd, zd;
} PhaseState;

typedef struct
{
  double a, e, i, longnode, argperi, meananom;
} OrbitalElements;

#define GG 6.6732e-8      // grav constant cgs
#define CC 2.997924562e10 // speed o light, cm/s
#define Msol 1.989e33     // mass o sun,  g
#define AU 1.49597892e13  // astron unit, cm

/* in kepcart.cpp */

// handing probability distributions randoms
double rayleigh(double sigma);
double powerlaw(double xmin, double xmax, double gamma);

// kepcart conversions
double ecc_ano(double e, double l);    // solves kepler's eqn
double ecc_anohyp(double e, double l); // solves kepler's eqn
void keplerian(double GM, PhaseState state, OrbitalElements *orbel);
void cartesian(double GM, OrbitalElements orbel, PhaseState *state);
double kepler(double ecc, double mean_anom);

// for integrating with f,g functions
void kepstep(double t, double M1,
             double x, double y, double z,
             double vx, double vy, double vz,
             double *xnew, double *ynew, double *znew,
             double *vxnew, double *vynew, double *vznew);
double solvex(double r0dotv0, double alpha,
              double M1, double r0, double dt);
double C_prussing(double y);
double S_prussing(double y);

#endif
