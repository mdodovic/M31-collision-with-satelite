#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>

#include "vect.h"
#include "consts.h"
#include "galaxy.h"
#include "bodies.h"
#include "node.h"

int elements = 0;

Vect Galaxy1Pos(-200 * distConversion, 0 * distConversion, 0);
Vect Galaxy1Vel(100000, 0, 0);
Vect Galaxy2Pos(200 * distConversion, 0 * distConversion, -15 * distConversion);
Vect Galaxy2Vel(-100000, 0, 0);
Vect InitialCenter(0, 0, 0);

double LastSnapTime = 0;

ifstream input;
ofstream out;
ofstream nodeOut;
ofstream setup;
ofstream energy;

string s;
string snp = "snapshot_";

double perc, currPerc, frame;
Vect VectDist;

inline double cube(double a)
{
	return a * a * a;
}

inline double sqr(double a)
{
	return a * a;
}

double PEn(Body **b)
{
	double P = 0;

	for (int i = 0; i < elements; ++i)
	{
		for (int j = 0; j < elements; ++j)
		{
			if (i != j)
			{
				Vect d = b[i]->r - b[j]->r;
				P -= G * b[i]->m * b[j]->m / d.Int();
			}
		}
	}
	return P;
}

double GetSize(Body **bodies)
{
	if (bodies == NULL)
	{
		return 0;
	}
	else
	{
		double sizemax = 0;
		for (int i = 0; i < elements; i++)
		{
			if (abs(bodies[i]->r.x) > sizemax)
				sizemax = abs(bodies[i]->r.x);
			if (abs(bodies[i]->r.y) > sizemax)
				sizemax = abs(bodies[i]->r.y);
			if (abs(bodies[i]->r.z) > sizemax)
				sizemax = abs(bodies[i]->r.z);
		}
		//added abs to return
		return sizemax;
	}
}

int snapshotCount = 0;

int main()
{
	double lastPercTime = 0;
	int currPerc = 0;

	input.open("galaxy");
	out.open("out.txt");
	nodeOut.open("node.txt");
	setup.open("setup.txt");
	energy.open("energy.txt");

	string s;
	string snp = "Snapshots/snapshot_";
	getline(input, s);
	int disk, halo, bulge;
	elements = 0;
	sscanf(s.c_str(), "%d", &disk);
	getline(input, s);
	sscanf(s.c_str(), "%d", &halo);
	getline(input, s);
	sscanf(s.c_str(), "%d", &bulge);

	elements = disk + halo + bulge;
	setup << elements << "\t";

	Body *bodies[elements];

	for (int i = 0; i < elements; i++)
	{
		bodies[i] = new Body();
	}

	for (int i = 0; i < elements; ++i)
	{
		double unM = 0;
		Vect unR, unV;
		getline(input, s);
		if (!s.empty())
		{
			sscanf(s.c_str(), "%lf %lf %lf %lf %lf %lf %lf", &unM, &unR.x, &unR.y, &unR.z, &unV.x, &unV.y, &unV.z);
		}
		if (i < elements / 2)
		{
			bodies[i]->m = unM * massConversion;
			bodies[i]->r = unR * distConversion;
			bodies[i]->v = unV * speedConversion;
		}
		else
		{
			bodies[i]->m = unM * massConversion;
			bodies[i]->r = unR * distConversion;
			bodies[i]->v = unV * speedConversion;
		}
		if (i < elements / 2)
		{
			bodies[i]->r = bodies[i]->r.RotateY(60 * rad);
			bodies[i]->v = bodies[i]->v.RotateY(60 * rad);
			bodies[i]->r += Galaxy1Pos;
			bodies[i]->v += Galaxy1Vel;
		}
		else
		{
			bodies[i]->r = bodies[i]->r.RotateY(120 * rad);
			bodies[i]->v = bodies[i]->v.RotateY(120 * rad);
			bodies[i]->r += Galaxy2Pos;
			bodies[i]->v += Galaxy2Vel;
		}
	}

	for (double t = 0; t <= tmax; t += dt)
	{
		//BARNES HUT!!!
		//HERE BE DRAGONS!

		double Nsize = GetSize(bodies);
		Node *root = new Node(InitialCenter, Nsize);
		for (int i = 0; i < elements; ++i)
		{
			root->InsertToNode(bodies[i]);
		}

		root->ComputeMassDistribution();

		if (t / tmax >= lastPercTime + 0.01)
		{
			cout << ++currPerc << "\r";
			lastPercTime = currPerc / 100.0;
			fflush(stdout);
		}

		if (t > LastSnapTime + snapshot || t == 0)
		{
			// cout << "pls" << endl;
			LastSnapTime = t;
			ofstream snap;
			s = std::to_string(snapshotCount);
			snap.open(snp + s);
			for (int i = 0; i < elements; i++)
			{
				snap << *bodies[i] << endl;
			}
			snapshotCount++;
			snap.close();

			double KinEn = 0;

			for (int i = 0; i < elements; i++)
			{
				KinEn += bodies[i]->m * bodies[i]->v.Int() * bodies[i]->v.Int() / 2;
			}
			double P = PEn(bodies);
			energy << t << "\t" << KinEn << "\t" << P << endl;
		}

		for (int i = 0; i < elements; ++i)
		{
			bodies[i]->F = root->CalcForce(bodies[i]);
			bodies[i]->a = bodies[i]->F / bodies[i]->m;
		}
		for (int i = 0; i < elements; i++)
		{
			bodies[i]->v += bodies[i]->a * dt;
			bodies[i]->r += bodies[i]->v * dt;
		}

		//flush
		for (int i = 0; i < elements; i++)
		{
			bodies[i]->a.x = 0;
			bodies[i]->a.y = 0;
			bodies[i]->a.z = 0;
		}
		delete root;
	}

	setup << snapshotCount;

	out.close();
	nodeOut.close();
	setup.close();
	energy.close();

	cout << "Done!" << endl;
	return 0;
}
