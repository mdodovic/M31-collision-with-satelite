#ifndef __GALAXY_H__
#define __GALAXY_H__

class Body;

void CreateGalaxy(Body **B, double GalaxyMass, Vect GalaxyPos, Vect GalaxyVel, double Scale, int Direction, int Start);
void CreatePointMass(Body **B, double GalaxyMass, Vect GalaxyPos, Vect GalaxyVel, int Start);
Vect GalaxySpeed(double Mass1, double Mass2, Vect Pos1, Vect Pos2, double Rmin);
#endif
