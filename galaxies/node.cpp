#include <iostream>
#include <stdlib.h>
#include "vect.h"
#include "consts.h"
#include "bodies.h"
#include <fstream>

const int TNE = 7;
const int TNW = 6;
const int TSE = 5;
const int TSW = 4;
const int BNE = 3;
const int BNW = 2;
const int BSE = 1;
const int BSW = 0;

int NodeCounter = 0;

//Contructor
Node::Node(Vect _center, double _size)
{
    center = _center;
    size = _size; // + (_size * 0.1);
    // cout << size << endl;
    for (int i = 0; i < 8; i++)
    {
        Octants[i] = NULL;
    }
    NodeCounter++;
}

//Destructor. Deletes each subnode.
Node::~Node()
{
    mass = 0;
    this->NodeBodyCount = 0;
    NodeBody = NULL;
    NodeCounter--;
    // NodeCounter--;

    // cout << center << endl;

    for (int i = 0; i < 8; i++)
    {
        delete Octants[i];
    }
}

//Inserts a body in to a Node, as described above each case a body may encounter.
void Node::InsertToNode(Body *NewBody)
{
    //If node has more bodies, check if desired node exists.
    //If it does, InsertToNode there, if not, create 4 subnodes,
    //then insert to desired node.
    if (this->NodeBodyCount > 1)
    {
        int oct = GetOct(NewBody);
        if (Octants[oct] == NULL)
        {
            Vect NewCenter = (GetCenter(oct) * (size / 2.0)) + center;
            Octants[oct] = new Node(NewCenter, (size / 2.0));
        }
        Octants[oct]->InsertToNode(NewBody);
    }
    else
    {
        //If node has one body, find it's new quadrant. If the new quadrant exists, place it there. If it doesn't,
        //create it and place it there.
        //Then take the new body and find it's own quadrant and check if it exists. If it does, place it there,
        //if not, create it and place it there.
        if (this->NodeBodyCount == 1)
        {
            int oct = GetOct(NodeBody);
            if (Octants[oct] == NULL)
            {
                Vect Smece = GetCenter(oct);
                // cout << Smece << "\t" << size << endl;
                Vect NewCenter = (GetCenter(oct) * (size / 2.0)) + center;
                Octants[oct] = new Node(NewCenter, (size / 2.0));
            }
            Octants[oct]->InsertToNode(NodeBody);

            oct = GetOct(NewBody);
            if (Octants[oct] == NULL)
            {
                Vect NewCenter = (GetCenter(oct) * (size / 2.0)) + center;
                Octants[oct] = new Node(NewCenter, size / 2.0);
            }
            Octants[oct]->InsertToNode(NewBody);

            NodeBody = NULL; //Make sure the previous body is no longer a leaf.
        }
        //If the node is empty, store the new body as the Node Body
        else
        {
            // cout << size;
            NodeBody = NewBody;
        }
    }
    //Increase the number of particles in the node;
    this->NodeBodyCount++;
}

//Not used, useful for debugging.
//Recursively counts bodies that are in a node.
// int Node::BHcount()
// {
//    int num = 0;
//    if (ne != NULL)
//    {
//       num += ne->BHcount();
//    }
//    if (nw != NULL)
//    {
//       num += nw->BHcount();
//    }
//    if (se != NULL)
//    {
//       num += se->BHcount();
//    }
//    if (sw != NULL)
//    {
//       num += sw->BHcount();
//    }
//    return num;
// }

//Returns the quadrant in which the forwarded body belongs.
//Returns number from 0 to 7.
int Node::GetOct(Body *body)
{
    // if (body->r.z > this->center.z)
    // {
    //    if (body->r.x > this->center.x)
    //    {
    //       if (body->r.y > this->center.y)
    //       {
    //          return 0; //TNE
    //       }
    //       else
    //       {
    //          return 1; //TSE
    //       }

    //    }
    //    else
    //    {
    //       if (body->r.y > this->center.y)
    //       {
    //          return 2; //TNW
    //       }
    //       else
    //       {
    //          return 3; //TSW
    //       }
    //    }
    // }
    // else
    // {
    //    if (body->r.x > this->center.x)
    //    {
    //       if (body->r.y > this->center.y)
    //       {
    //          return 4; //BNE
    //       }
    //       else
    //       {
    //          return 5; //BSE
    //       }

    //    }
    //    else
    //    {
    //       if (body->r.y > this->center.y)
    //       {
    //          return 6; //BNW
    //       }
    //       else
    //       {
    //          return 7; //BSW
    //       }
    //    }
    // }
    Vect posRelToNode = body->r - this->center;
    // cout << posRelToNode << '\t' << this << endl;
    return (posRelToNode.x > 0) |
           ((posRelToNode.y > 0) << 1) |
           ((posRelToNode.z > 0) << 2);
}

//If the node has 1 particle, set the node mass to that of the particle,
//and set its center of mass at the position of the particle.
//If the node has more than one particle, compute the mass distribution for all of the nodes existing and non-empty
//child nodes.
void Node::ComputeMassDistribution()
{
    if (this->NodeBodyCount == 1)
    {
        this->mass = NodeBody->m;
        this->centerOfMass = NodeBody->r;
    }
    else
    {
        double masses[8];
        for (int i = 0; i < 8; ++i)
        {
            masses[i] = 0;
        }

        for (int i = 0; i < 8; ++i)
        {
            if (Octants[i] != NULL)
            {
                Octants[i]->ComputeMassDistribution();
                masses[i] = Octants[i]->mass;
                centerOfMass += Octants[i]->centerOfMass * Octants[i]->mass;
            }
        }

        for (int i = 0; i < 8; ++i)
        {
            this->mass += masses[i];
        }

        this->centerOfMass /= this->mass;
    }
}

// Main point of the algorhythm.
// Calculates the force exerted on a body by traversing the BH tree and checking whether approximations can me used.
Vect Node::CalcForce(Body *calculatedBody)
{
    Vect zeros;
    Vect Force;
    //If the node has only one body, calculate the force acting between two bodies.
    if (this->NodeBodyCount == 1)
    {
        if (calculatedBody != NodeBody)
        {
            Vect dist = calculatedBody->r - NodeBody->r;
            if (dist.Int() < damp)
                return zeros;
            Force = dist * (-G * calculatedBody->m * this->mass / dist.CubeInt());
            return Force;
        }
        else
        {
            return zeros;
        }
    }
    //If the node has more than one body check whether approximations can be used.
    else
    {
        if (NodeBodyCount != 0)
        {
            Vect distFromCentOfMass = calculatedBody->r - this->centerOfMass;
            double distInt = distFromCentOfMass.Int();
            //If the body is sufficiently far away from the node, use the node as an approximation of a body
            //with the same mass as all of it's bodies, and the center of mass as its position.
            if (size / distInt < theta)
            {
                Vect dist;
                dist = calculatedBody->r - this->centerOfMass;
                if (dist.Int() < damp)
                    return zeros;
                Force = dist * (-G * calculatedBody->m * this->mass / dist.CubeInt());
                return Force;
            }
            //If the parameter isn't met, calculate the force acting within each subnode and sum it up.
            //Recursion is magic.
            else
            {
                // Vect Forces[8];
                for (int i = 0; i < 8; ++i)
                {
                    if (Octants[i] != NULL)
                    {
                        Force += Octants[i]->CalcForce(calculatedBody);
                    }
                }
                // for (int i = 0; i < 8; ++i)
                // {
                //    Force += Forces[i];
                // }

                return Force;
            }
        }
        else
        {
            Vect zeros;
            return zeros;
        }
    }
}

//Overloaded operator that allows the user to recursively output the center and size of each node and with subnodes.
ostream &operator<<(ostream &O, Node &B)
{
    O << B.center << "\t" << B.size << endl;
    // if (B.ne != NULL)
    //    O << *B.ne;
    // if (B.nw != NULL)
    //    O << *B.nw;
    // if (B.sw != NULL)
    //    O << *B.sw;
    // if (B.se != NULL)
    //    O << *B.se;

    for (int i = 0; i < 8; ++i)
    {
        if (B.Octants[i] != NULL)
        {
            O << *B.Octants[i];
        }
    }
    return O;
}

Vect GetCenter(int octant)
{
    Vect ort;
    int oct = octant;

    if (oct == TNE || oct == TSE || oct == BNE || oct == BSE)
        ort.x = 1;
    else
        ort.x = -1;

    if (oct == TNE || oct == TNW || oct == BNE || oct == BNW)
        ort.y = 1;
    else
        ort.y = -1;

    if (oct == TNE || oct == TSE || oct == TSW || oct == TNW)
        ort.z = 1;
    else
        ort.z = -1;

    return ort;
}
