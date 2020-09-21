#ifndef consts_h
#define consts_h

const double theta = 0.6;

const double G = 6.67e-11;
const double tmax = 5e9 * (60 * 60 * 24 * 365.25);
const double dt = 0.5e6 * (60 * 60 * 24 * 365.25);

const double distConversion = 3.08567758e19;
const double speedConversion = 100000;
const double massConversion = 4.65e39;

const double snapshot = 100 * dt;

const double Rmin = 25 * 3.08567758e19;

const double rad = 0.0174532925;

const double Galaxy1Mass = 1e11 * 2e30;
const double Galaxy2Mass = 0.5 * 1e11 * 2e30;

const double damp = 50 * 3.08567758e16;

#endif

// Vect Gr1(Rmin, Rmin, Rmin);
// Vect Gr2(-Rmin, Rmin, -Rmin);
// Vect Gr4(-Rmin, -Rmin, 0);
// Vect Gr3(Rmin, -2 * Rmin, Rmin);

// bodies[0]->r = Gr1;
// bodies[1]->r = Gr2;
// bodies[2]->r = Gr3;
// bodies[3]->r = Gr4;

// bodies[0]->v = Galaxy1Vel * -5;
// bodies[1]->v = Galaxy1Vel * -3;
// bodies[2]->v = Galaxy2Vel * -7;
// bodies[3]->v = Galaxy2Vel * -1;
//debug conditions
