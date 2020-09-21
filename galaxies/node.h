#ifndef __NODE_H__
#define __NODE_H__

#include <iostream>
#include <fstream>
#include "vect.h"
#include "bodies.h"

extern int NodeCounter;
class Body;

Vect GetCenter(int octant);

class Node
{
public:
	Vect center;
	bool Contains(Vect r);
	double size;
	Node *Octants[8];
	Body *NodeBody;
	int NodeBodyCount;
	Vect centerOfMass;
	double mass;

	Node(Vect _center, double size);
	~Node();
	void InsertToNode(Body *NewBody);
	void ComputeMassDistribution();
	int BHcount();
	int GetOct(Body *body);
	Vect CalcForce(Body *calculatedBody);
	//double GetSize(Body **bodies);
};

ostream &operator<<(ostream &O, Node &B);
#endif