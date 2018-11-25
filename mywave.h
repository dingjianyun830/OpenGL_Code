#pragma once
#include "myray.h"
#include "glm.h"

double thetas[2] = { 0.38,1.42 };
double amplitudes[3][2] = { 0.2,0.2,0.3,0.5,0.2,0.6 };
double omegas[3] = { 3.27,3.31,3.42 };
double ks[3] = { 1.091,1.118,1.1935 };

void wave(Vect v1, Vect v, double time)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			int k = omegas[i] * omegas[i] / 9.8;
			v1[X] -= cos(thetas[j])*amplitudes[i][j] * sin(k * (v[X]*cos(thetas[j]) + v[Z]*sin(thetas[j])) - omegas[i] * time);
			v1[Y] += amplitudes[i][j] * cos(k * (v[X]*cos(thetas[j]) + v[Z]*sin(thetas[j])) - omegas[i] * time);
			v1[Z] -= sin(thetas[j])*amplitudes[i][j] * sin(k * (v[X]*cos(thetas[j]) + v[Z]*sin(thetas[j])) - omegas[i] * time);
		}
	}
}
void myWave(GLMmodel *model, double time)
{
	static GLMgroup* group;
	static GLMtriangle* triangle;
	Vect V0, V1, V2;
	for (group = model->groups; group; group = group->next)
	{

		for (int i = 0; i < group->numtriangles; i++)
		{
			triangle = &model->triangles[group->triangles[i]];

			// get triangle vertices

			V0[X] = model->vertices[3 * triangle->vindices[0]];
			V0[Y] = model->vertices[3 * triangle->vindices[0] + 1];
			V0[Z] = model->vertices[3 * triangle->vindices[0] + 2];

			V1[X] = model->vertices[3 * triangle->vindices[1]];
			V1[Y] = model->vertices[3 * triangle->vindices[1] + 1];
			V1[Z] = model->vertices[3 * triangle->vindices[1] + 2];

			V2[X] = model->vertices[3 * triangle->vindices[2]];
			V2[Y] = model->vertices[3 * triangle->vindices[2] + 1];
			V2[Z] = model->vertices[3 * triangle->vindices[2] + 2];

			wave(V0, V0, time);
			wave(V1, V1, time);
			wave(V2, V2, time);

			model->vertices[3 * triangle->vindices[0]] = V0[X];
			model->vertices[3 * triangle->vindices[0] + 1] = V0[Y];
			model->vertices[3 * triangle->vindices[0] + 2] = V0[Z];

			model->vertices[3 * triangle->vindices[1]] = V1[X];
			model->vertices[3 * triangle->vindices[1] + 1] = V1[Y];
			model->vertices[3 * triangle->vindices[1] + 2] = V1[Z];

			model->vertices[3 * triangle->vindices[2]] = V2[X];
			model->vertices[3 * triangle->vindices[2] + 1] = V2[Y];
			model->vertices[3 * triangle->vindices[2] + 2] = V2[Z];
		}
	}
}