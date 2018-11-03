#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
double distance(const double * p1, const double * p2, int dim = 3)
{
	double tot = 0;
	for (int i = 0; i < dim; i++)
	{
		tot += (p1[i] - p2[i]) * (p1[i] - p2[i]);
	}
	return sqrt(tot);
}

double elev(const double * sat_pos, const double * user_pos)
{

	double dpos[3] = { 0 };
	double ori[3]{ 0,0,0 };
	dpos[0] = sat_pos[0] - user_pos[0];
	dpos[1] = sat_pos[1] - user_pos[1];
	dpos[2] = sat_pos[2] - user_pos[2];

	double user_distance_to_earth = distance(user_pos, ori);

	double mod = sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2]);
	if (fabs(user_distance_to_earth * mod < 1.0)) {
		return M_PI_2;
	}
	else {
		double m = dpos[0] * user_pos[0] + dpos[1] * user_pos[1] + dpos[2] * user_pos[2];
		double n = m / (mod * user_distance_to_earth);
		return M_PI_2 - acos(n);
	}
}
