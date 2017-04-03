#pragma once


#include <nclgl\Matrix3.h>
#include <nclgl\Matrix4.h>
#include <nclgl\Vector2.h>
#include <vector>

#ifndef uint
typedef unsigned int uint;
#endif

struct FETriangle
{
	FETriangle() {}
	FETriangle(uint a, uint b, uint c, uint d, uint e, uint f)
	{
		v1 = a;
		v2 = b; 
		v3 = c;
		v4 = d;
		v5 = e;
		v6 = f;

		for (int i = 0; i < 9; ++i)
		{
			tangents[i] = 0;
			tan_multipliers[i] = 1.f;
		}
	}

	union
	{
		struct
		{
			uint v1, v2, v3, v4, v5, v6;
			uint t12, t13, t21, t23, t31, t32, tab, tbc, tca;
			float t12_mult, t13_mult, t21_mult, t23_mult, t31_mult, t32_mult, tab_mult, tbc_mult, tca_mult;
		};

		struct
		{
			uint phyxels[6]; //outer i,j,k -> inner ij, jk, ki
			uint tangents[9];
			float tan_multipliers[9];
		};
	};
};

struct FEVertDescriptor
{
	Vector3 pos;
	Vector3 ipos;
	Vector2 tCoord;
	bool isStatic;
};


#define GRAVITY				Vector3(0, -9.81f, 0)
#define STATIC_MASS_SCALING 1000.0f

#define WOOL_VISCOSE		0
#define WOOL				1
#define ACETATE				2
#define POLYESTER			3
#define SHEAR_RESISTANT		4
#define BEND_RESISTANT		5

#define CHOSEN_MATERIAL		WOOL



//const float V_SCALAR = 0.012f;		//Scalar of C to produce Voidt Viscosity tensor
const float V_SCALAR = 0.002f;		//Scalar of C to produce Voidt Viscosity tensor
const Vector3 WIND = Vector3(0, 0, 0);// Vector3(20.5f, 0.0f, 1.5f);

const float max_timestep = 1.0f / 920.0f;
const float min_timestep = max_timestep * 0.125 * 0.25f;



#if CHOSEN_MATERIAL == WOOL_VISCOSE
const float mass_density = 0.23;	// kg/m^2
const float C1111 = 245;			// Weft Stretch
const float C2222 = 366;			// Warp Stretch
const float C1212 = 2.38;			// Shear Modulus
const float C1122 = 61.1;			// Transverse Contraction
const float B1 = 0.013 * 10000;	// Weft Bend
const float B2 = 0.037 * 10000;	// Warp Bend
#endif




								///-------------------------------------
#if CHOSEN_MATERIAL == WOOL

								//const float mass_density = 0.26f;	// kg/m^2

								//YUSEF FABRIC 1

								/*const float C1111 = 2252.f;			// Weft Stretch		(0.0014)
								const float C2222 = 1249.f;			// Warp Stretch		(0.0025)
								const float C1212 = 0.8f;			// Shear Modulus
								const float C1122 = 225.7f;			// Transverse Contraction
								const float mass_density = 0.127f;	// kg/m^2
								const float B1 = 0.00999f ;	// Weft Bend
								const float B2 = 0.004041f;	// Warp Bend*/

								/*
								//Dropping1
								const float C1111 = 3782.f;			// Weft Stretch		(0.0014)
								const float C2222 = 3320.f;			// Warp Stretch		(0.0025)
								const float C1212 = 0.0f;			// Shear Modulus
								const float C1122 = 225.7f;			// Transverse Contraction
								const float mass_density = 0.127f;	// kg/m^2
								const float B1 = 0.0999f;	// Weft Bend
								const float B2 = 0.05041f;	// Warp Bend*/

								//Dropping2
const float C1111 = 1670.f;			// Weft Stretch		(0.0014)
const float C2222 = 2665.f;			// Warp Stretch		(0.0025)
const float C1212 = 0.2f;			// Shear Modulus
const float C1122 = 225.7f;			// Transverse Contraction
const float mass_density = 0.165f;	// kg/m^2
const float B1 = 0.106f;	// Weft Bend
const float B2 = 0.08f;	// Warp Bend

						/*
						//YUSEF FABRIC 2
						const float C1111 = 549.2f;			// Weft Stretch		(0.0170)
						const float C2222 = 1202.f;			// Warp Stretch		(0.0078)
						const float C1212 = 0.8f;			// Shear Modulus
						const float C1122 = 225.7f;			// Transverse Contraction
						const float mass_density = 0.164f;	// kg/m^2
						const float B1 = 0.106f;	// Weft Bend
						const float B2 = 0.08041f;	// Warp Bend
						*/

#endif
						//----------------------------------------




#if CHOSEN_MATERIAL == ACETATE
const float mass_density = 0.17;	// kg/m^2
const float C1111 = 3057;			// Weft Stretch
const float C2222 = 1534;			// Warp Stretch
const float C1212 = 1.22;			// Shear Modulus
const float C1122 = 459.1;			// Transverse Contraction
const float B1 = 0.055;// * 10000;		// Weft Bend
const float B2 = 0.092;// * 10000;		// Warp Bend
#endif
#if CHOSEN_MATERIAL == POLYESTER
const float mass_density = 0.26;	// kg/m^2
const float C1111 = 2400;			// Weft Stretch
const float C2222 = 3600;			// Warp Stretch
const float C1212 = 5.23;			// Shear Modulus
const float C1122 = 600;			// Transverse Contraction
const float B1 = 0.371 * 10000;		// Weft Bend
const float B2 = 0.480 * 10000;		// Warp Bend
#endif	
#if CHOSEN_MATERIAL == SHEAR_RESISTANT
const float mass_density = 0.2;	// kg/m^2
const float C1111 = 2000;			// Weft Stretch
const float C2222 = 2000;			// Warp Stretch
const float C1212 = 1000;			// Shear Modulus
const float C1122 = 400;			// Transverse Contraction
const float B1 = 100 * 10000;		// Weft Bend
const float B2 = 100 * 10000;		// Warp Bend
#endif

#if CHOSEN_MATERIAL == BEND_RESISTANT
const float mass_density = 0.2;	// kg/m^2
const float C1111 = 2000;			// Weft Stretch
const float C2222 = 2000;			// Warp Stretch
const float C1212 = 5.0;			// Shear Modulus
const float C1122 = 400;			// Transverse Contraction
const float B1 = 100.0 * 10000;		// Weft Bend
const float B2 = 100.0 * 10000;		// Warp Bend
#endif

const float C2121 = C1212;			// Shear Modulus
const float C2112 = C1212;
const float C1221 = C1212;

const float C2211 = C1122;			// Transverse Contraction

const float C2122 = 0.f;
const float C1211 = 0.f;
const float C1112 = 0.f;
const float C2111 = 0.f;
const float C2221 = 0.f;
const float C2212 = 0.f;
const float C1121 = 0.f;
const float C1222 = 0.f;




















