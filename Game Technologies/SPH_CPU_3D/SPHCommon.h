#pragma once
#include <nclgl\Vector3.h>

#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sort.h>

#include <vector>
#include <list>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef unsigned int uint;

/*struct cpuint3
{
	int x;
	int y;
	int z;
};

struct uint3
{
	uint x;
	uint y;
	uint z;
};*/

class Particle
{
public:
	uint id;
	Vector3 pos;
	Vector3 vel;
	Vector3 acc;
	Vector3 ev;

	float dens;
	float pres;

	float surf_norm;

	Particle *next;
};