
/*
   keplerian conversion routines 
   some old code modified from Holman's original
   added f+g functions
   generalized to hyperbolic
   also has some distribution function random generators
*/

/* routines:
   keplerian: returns orbital elements from cartesian coords
   cartesian: returns cartesian from orbital elements 
   ecc_ano:   solve Kepler's equation eccentric orbits
   ecc_anohyp:   solve Kepler's equation hyperbolic orbits
   kepler:    solve Kepler's equation general orbits
   kepstep:   do a kepstep using f+g functions
   solvex     needed by kepstep, universal coordinate diff kepler's eqn solver
   C_prussing, S_prussing needed by solvex
   rayleigh   Rayleigh probability distribution
   powerlaw   powerlaw probability distribution
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "kepcart.h"

#define SIGN(a) ((a) < 0 ? -1 : 1)

// parabolic case not implemented
// double kepler(double ecc, double mean_anom);

/* given PhaseState returns orbital elements */
void keplerian(double GM, PhaseState state, OrbitalElements *orbel)
{
  double rxv_x, rxv_y, rxv_z, hs, h;
  double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
  double rcosu, rsinu, u, trueanom, eccanom;

  /* find direction of angular momentum vector */
  rxv_x = state.y * state.zd - state.z * state.yd;
  rxv_y = state.z * state.xd - state.x * state.zd;
  rxv_z = state.x * state.yd - state.y * state.xd;
  hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z;
  h = sqrt(hs);

  r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
  vs = state.xd * state.xd + state.yd * state.yd + state.zd * state.zd;

  rdotv = state.x * state.xd + state.y * state.yd + state.z * state.zd;
  rdot = rdotv / r;

  orbel->i = acos(rxv_z / h);

  if (rxv_x != 0.0 || rxv_y != 0.0)
  {
    orbel->longnode = atan2(rxv_x, -rxv_y);
  }
  else
    orbel->longnode = 0.0;

  orbel->a = 1.0 / (2.0 / r - vs / GM); // could be negative

  ecostrueanom = hs / (GM * r) - 1.0;
  esintrueanom = rdot * h / GM;
  orbel->e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom);

  if (esintrueanom != 0.0 || ecostrueanom != 0.0)
  {
    trueanom = atan2(esintrueanom, ecostrueanom);
  }
  else
    trueanom = 0.0;

  cosnode = cos(orbel->longnode);
  sinnode = sin(orbel->longnode);

  /* u is the argument of latitude */
  rcosu = state.x * cosnode + state.y * sinnode;
  rsinu = (state.y * cosnode - state.x * sinnode) / cos(orbel->i);

  if (rsinu != 0.0 || rcosu != 0.0)
  {
    u = atan2(rsinu, rcosu);
  }
  else
    u = 0.0;

  orbel->argperi = u - trueanom;

  double foo = sqrt(fabs(1.0 - orbel->e) / (1.0 + orbel->e));
  if (orbel->e < 1.0)
  {
    eccanom = 2.0 * atan(foo * tan(trueanom / 2.0));
    orbel->meananom = eccanom - orbel->e * sin(eccanom);
    if (orbel->meananom > M_PI)
      orbel->meananom -= 2.0 * M_PI;
    if (orbel->meananom < -M_PI)
      orbel->meananom += 2.0 * M_PI;
    // only shift M if elliptic orbit
  }
  else
  {
    eccanom = 2.0 * atanh(foo * tan(trueanom / 2.0));
    orbel->meananom = orbel->e * sinh(eccanom) - eccanom;
  }
  if (orbel->argperi > M_PI)
    orbel->argperi -= 2.0 * M_PI;
  if (orbel->argperi < -M_PI)
    orbel->argperi += 2.0 * M_PI;
}

/* given orbital elements return PhaseState */
void cartesian(double GM, OrbitalElements orbel, PhaseState *state)
{
  double meanmotion, cosE, sinE, foo;
  double x, y, z, xd, yd, zd;
  double xp, yp, zp, xdp, ydp, zdp;
  double cosw, sinw, cosi, sini, cosnode, sinnode;
  double E0, rovera;
  double a = orbel.a;
  double e = orbel.e;
  double i = orbel.i;
  double longnode = orbel.longnode;
  double argperi = orbel.argperi;
  double meananom = orbel.meananom;
  /* double E1, E2, den; */

  /* compute eccentric anomaly */
  if (e < 1)
    E0 = ecc_ano(e, meananom);
  else
    E0 = ecc_anohyp(e, meananom);
  // E0 = kepler(e,meananom); // also works

  if (e < 1.0)
  {
    cosE = cos(E0);
    sinE = sin(E0);
  }
  else
  {
    cosE = cosh(E0);
    sinE = sinh(E0);
  }
  a = fabs(a);
  meanmotion = sqrt(GM / (a * a * a));
  foo = sqrt(fabs(1.0 - e * e));
  /* compute unrotated positions and velocities */
  rovera = (1.0 - e * cosE);
  if (e > 1.0)
    rovera *= -1.0;
  x = a * (cosE - e);
  y = foo * a * sinE;
  z = 0.0;
  xd = -a * meanmotion * sinE / rovera;
  yd = foo * a * meanmotion * cosE / rovera;
  zd = 0.0;
  if (e > 1.0)
    x *= -1.0;

  /* rotate by argument of perihelion in orbit plane*/
  cosw = cos(argperi);
  sinw = sin(argperi);
  xp = x * cosw - y * sinw;
  yp = x * sinw + y * cosw;
  zp = z;
  xdp = xd * cosw - yd * sinw;
  ydp = xd * sinw + yd * cosw;
  zdp = zd;

  /* rotate by inclination about x axis */
  cosi = cos(i);
  sini = sin(i);
  x = xp;
  y = yp * cosi - zp * sini;
  z = yp * sini + zp * cosi;
  xd = xdp;
  yd = ydp * cosi - zdp * sini;
  zd = ydp * sini + zdp * cosi;

  /* rotate by longitude of node about z axis */
  cosnode = cos(longnode);
  sinnode = sin(longnode);
  state->x = x * cosnode - y * sinnode;
  state->y = x * sinnode + y * cosnode;
  state->z = z;
  state->xd = xd * cosnode - yd * sinnode;
  state->yd = xd * sinnode + yd * cosnode;
  state->zd = zd;
}

/* ----------------Solve Kepler's equation
 iterate to get an estimate for eccentric anomaly (u) given the mean anomaly (l).
Appears to be accurate to level specified, I checked this 
and it works for u,lambda in all quadrants and at high eccentricity
*/
#define PREC_ecc_ano 1e-14 /* no reason that this must be very accurate in code at present */
double ecc_ano(double e, double l)
{
  double du, u0, l0;
  du = 1.0;
  u0 = l + e * sin(l) + 0.5 * e * e * sin(2.0 * l);
  // also see M+D equation 2.55
  /* supposed to be good to second order in e, from Brouwer+Clemence
                    u0 is first guess */
  while (fabs(du) > PREC_ecc_ano)
  {
    l0 = u0 - e * sin(u0);
    du = (l - l0) / (1.0 - e * cos(u0));
    u0 += du; /* this gives a better guess */
              // equation 2.58 from M+D
  }
  return u0;
}
// hyperbolic case
double ecc_anohyp(double e, double l)
{
  double du, u0, fh, dfh;
  du = 1.0;
  u0 = log(2.0 * l / e + 1.8); //danby guess
  while (fabs(du) > PREC_ecc_ano)
  {
    fh = e * sinh(u0) - u0 - l;
    dfh = e * cosh(u0) - 1.0;
    du = -fh / dfh;
    u0 += du;
  }
  return u0;
}

// return a random number r that is consistent with the
// Rayleigh probability distribution
// p(r) = (r/sigma^2) exp[-r^2/(2 sigma^2))]
double rayleigh(double sigma)
{
  double y = 0.0;
  double r2;
  y = (double)rand() / (double)RAND_MAX;
  r2 = -2.0 * log(y);
  // printf("rayleigh %.2e %.2e %.2e\n",y,r2,sigma);
  return sigma * sqrt(r2);
}

// return a random number x for p(x)
// a powerlaw distribution that goes between xmin, xmax with exponent gamma
// p(x) propto x^gamma
double powerlaw(double xmin, double xmax, double gamma)
{
  double diff, y, x;
  double u = (double)rand() / (double)RAND_MAX;
  if (gamma == -1.0)
  {
    x = xmin * pow(xmax / xmin, u);
  }
  else
  {
    double g1 = gamma + 1.0;
    diff = pow(xmax, g1) - pow(xmin, g1);
    y = u * diff + pow(xmin, g1);
    x = pow(y, 1.0 / g1);
  }
  return x;
}

// solving kepler's equation including the
// hyperbolic case, code downloaded from project pluto

// #include <math.h>

// #define PI 3.14159265358979323 in -lm
#define THRESH 1.e-8
#define CUBE_ROOT(X) (exp(log(X) / 3.))

//static double asinh( const double z);
// double kepler( const double ecc, double mean_anom);

// in -lm
// static double asinh( const double z)
// { return( log( z + sqrt( z * z + 1.))); }

double kepler(double ecc, double mean_anom)
{
  double curr, err, thresh;
  int is_negative = 0, n_iter = 0;

  if (!mean_anom)
    return (0.);

  if (ecc < .3) /* low-eccentricity formula from Meeus,  p. 195 */
  {
    curr = atan2(sin(mean_anom), cos(mean_anom) - ecc);
    /* one correction step,  and we're done */
    err = curr - ecc * sin(curr) - mean_anom;
    curr -= err / (1. - ecc * cos(curr));
    return (curr);
  }

  if (mean_anom < 0.)
  {
    mean_anom = -mean_anom;
    is_negative = 1;
  }

  curr = mean_anom;
  thresh = THRESH * fabs(1. - ecc);
  if (ecc > .8 && mean_anom < M_PI / 3. || ecc > 1.) /* up to 60 degrees */
  {
    double trial = mean_anom / fabs(1. - ecc);

    if (trial * trial > 6. * fabs(1. - ecc)) /* cubic term is dominant */
    {
      if (mean_anom < M_PI)
        trial = CUBE_ROOT(6. * mean_anom);
      else /* hyperbolic w/ 5th & higher-order terms predominant */
        trial = asinh(mean_anom / ecc);
    }
    curr = trial;
  }

  if (ecc < 1.)
  {
    err = curr - ecc * sin(curr) - mean_anom;
    while (fabs(err) > thresh)
    {
      n_iter++;
      curr -= err / (1. - ecc * cos(curr));
      err = curr - ecc * sin(curr) - mean_anom;
    }
  }
  else
  {
    err = ecc * sinh(curr) - curr - mean_anom;
    while (fabs(err) > thresh)
    {
      n_iter++;
      curr -= err / (ecc * cosh(curr) - 1.);
      err = ecc * sinh(curr) - curr - mean_anom;
    }
  }
  return (is_negative ? -curr : curr);
}

// double kepler( const double ecc, double mean_anom);
void kepstep(double t, double M1,
             double x, double y, double z,
             double vx, double vy, double vz,
             double *xnew, double *ynew, double *znew,
             double *vxnew, double *vynew, double *vznew);
double solvex(double r0dotv0, double alpha,
              double M1, double r0, double dt);
double C_prussing(double y);
double S_prussing(double y);

// given M1,x,y,z, vx,vy,vz, calculate new position and velocity at
// time t later
// using f and g functions and formulae from Prussing  + conway
// here M1 is actually G(M_1+M_2)

// do a keplerian time step using f,g functions.
// is robust, covering hyperbolic and parabolic as well as elliptical orbits
// M1 = GM
#define eps 1e-16
void kepstep(double dt, double M1,
             double x, double y, double z,
             double vx, double vy, double vz,
             double *xnew, double *ynew, double *znew,
             double *vxnew, double *vynew, double *vznew)
{
  double r0 = sqrt(x * x + y * y + z * z);   // current radius
  double v2 = (vx * vx + vy * vy + vz * vz); // current velocity
  double r0dotv0 = (x * vx + y * vy + z * vz);
  double alpha = (2.0 / r0 - v2 / M1);             // inverse of semi-major eqn 2.134 MD
                                                   // here alpha=1/a and can be negative
  double x_p = solvex(r0dotv0, alpha, M1, r0, dt); // solve universal kepler eqn
  double smu = sqrt(M1);
  double foo = 1.0 - r0 * alpha;
  double sig0 = r0dotv0 / smu;

  double x2, x3, alx2, Cp, Sp, r;
  x2 = x_p * x_p;
  x3 = x2 * x_p;
  alx2 = alpha * x2;
  Cp = C_prussing(alx2);
  Sp = S_prussing(alx2);
  r = sig0 * x_p * (1.0 - alx2 * Sp) + foo * x2 * Cp + r0; // eqn 2.42  PC
                                                           // f,g functions equation 2.38a  PC
  double f_p = 1.0 - (x2 / r0) * Cp;
  double g_p = dt - (x3 / smu) * Sp;
  // dfdt,dgdt function equation 2.38b PC
  double dfdt = x_p * smu / (r * r0) * (alx2 * Sp - 1.0);
  double dgdt = 1.0 - (x2 / r) * Cp;

  if (r0 > 0.0)
  {                             // error catch if a particle is at Sun
    *xnew = x * f_p + g_p * vx; // eqn 2.65 M+D
    *ynew = y * f_p + g_p * vy;
    *znew = z * f_p + g_p * vz;

    *vxnew = dfdt * x + dgdt * vx; //eqn 2.70 M+D
    *vynew = dfdt * y + dgdt * vy;
    *vznew = dfdt * z + dgdt * vz;
  }
  else
  {
    *xnew = x;
    *ynew = y;
    *znew = z;

    *vxnew = vx;
    *vynew = vy;
    *vznew = vz;
  }
}

  // use Laguerre method as outlined by Prusing+C eqn 2.43
  // return x universal variable
  // solving differential Kepler's equation
  // in universal variable

#define N_LAG 5.0 // integer for recommeded Laguerre method

double solvex(double r0dotv0, double alpha,
              double M1, double r0, double dt)
{
  double smu = sqrt(M1);
  double foo = 1.0 - r0 * alpha;
  double sig0 = r0dotv0 / smu;
  double x = M1 * dt * dt / r0; // initial guess could be improved

  double u = 1.0;
  for (int i = 0; i < 7; i++)
  { // while(fabs(u) > EPS){
    double x2, x3, alx2, Cp, Sp, F, dF, ddF, z;
    x2 = x * x;
    x3 = x2 * x;
    alx2 = alpha * x2;
    Cp = C_prussing(alx2);
    Sp = S_prussing(alx2);
    F = sig0 * x2 * Cp + foo * x3 * Sp + r0 * x - smu * dt; // eqn 2.41 PC
    dF = sig0 * x * (1.0 - alx2 * Sp) + foo * x2 * Cp + r0; // eqn 2.42 PC
    ddF = sig0 * (1.0 - alx2 * Cp) + foo * x * (1.0 - alx2 * Sp);
    z = fabs((N_LAG - 1.0) * ((N_LAG - 1.0) * dF * dF - N_LAG * F * ddF));
    z = sqrt(z);
    u = N_LAG * F / (dF + SIGN(dF) * z); // equation 2.43 PC
    x -= u;
    //     printf("s: %.3e\n",F); is good
  }
  return x;
}

// functions needed
double C_prussing(double y) // equation 2.40a Prussing + Conway
{
  if (fabs(y) < 1e-4)
    return 1.0 / 2.0 * (1.0 - y / 12.0 * (1.0 - y / 30.0 * (1.0 - y / 56.0)));
  double u = sqrt(fabs(y));
  if (y > 0.0)
    return (1.0 - cos(u)) / y;
  else
    return (cosh(u) - 1.0) / -y;
}

double S_prussing(double y) // equation 2.40b Prussing +Conway
{
  if (fabs(y) < 1e-4)
    return 1.0 / 6.0 * (1.0 - y / 20.0 * (1.0 - y / 42.0 * (1.0 - y / 72.0)));
  double u = sqrt(fabs(y));
  double u3 = u * u * u;
  if (y > 0.0)
    return (u - sin(u)) / u3;
  else
    return (sinh(u) - u) / u3;
}
