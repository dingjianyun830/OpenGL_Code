// raytracer.cpp : Defines the entry point for the console application.
//



#include <stdio.h>
#include <stdlib.h>
#include <GL\glut.h>
#include <iostream>

//#include "pch.h"


#include "myray.h"
#include "glm.h"
#include "mywave.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

extern Camera *ray_cam;       // camera info
extern int image_i, image_j;  // current pixel being shaded
extern bool wrote_image;      // has the last pixel been shaded?

							  // reflection/refraction recursion control

extern int maxlevel;          // maximum depth of ray recursion 
extern double minweight;      // minimum fractional contribution to color

							  // these describe the scene

extern vector < GLMmodel * > model_list;
extern vector < Sphere * > sphere_list;
extern vector < Light * > light_list;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// intersect a ray with the entire scene (.obj models + spheres)

// x, y are in pixel coordinates with (0, 0) the upper-left hand corner of the image.
// color variable is result of this function--it carries back info on how to draw the pixel

void trace_ray(int level, double weight, Ray *ray, Vect color)
{
	Intersection *nearest_inter = NULL;
	Intersection *inter = NULL;
	int i;

	// test for intersection with all .obj models

	for (i = 0; i < model_list.size(); i++) {
		inter = intersect_ray_glm_object(ray, model_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	if (nearest_inter) {
		//shade_ray_false_color_normal(nearest_inter, color);
		//shade_ray_intersection_mask(color);  
		//shade_ray_diffuse(ray, nearest_inter, color);
		//shade_ray_local(ray, nearest_inter, color);
		shade_ray_recursive(level, weight, ray, nearest_inter, color);
	}

	// color the ray using a default

	else
		shade_ray_background(ray, color);
}

//----------------------------------------------------------------------------
void shade_ray(Ray *ray, Intersection *inter, Vect color)
{
	Vect L, n;
	double diff_factor;

	// iterate over lights

	for (int i = 0; i < light_list.size(); i++) {

		// AMBIENT

		color[R] += inter->surf->amb[R] * light_list[i]->amb[R];
		color[G] += inter->surf->amb[G] * light_list[i]->amb[G];
		color[B] += inter->surf->amb[B] * light_list[i]->amb[B];

		// DIFFUSE 

		// INSERT YOUR CODE HERE
		VectSub(light_list[i]->P, inter->P, L);
		VectUnit(L);
		VectUnit(inter->N);
		double NdotL = VectDotProd(inter->N, L);
		color[R] += inter->surf->diff[R] * light_list[i]->diff[R] * NdotL;
		color[G] += inter->surf->diff[G] * light_list[i]->diff[G] * NdotL;
		color[B] += inter->surf->diff[B] * light_list[i]->diff[B] * NdotL;
	}

	// clamp color to [0, 1]

	VectClamp(color, 0, 1);
}
// test for ray-sphere intersection; return details of intersection if true

Intersection *intersect_ray_sphere(Ray *ray, Sphere *S)
{
	// INSERT YOUR CODE HERE (line below says "no" for all spheres, so replace it)

	Vect V0, V1, V2;
	Intersection *nearest_inter = NULL;
	Intersection *inter = NULL;
	
	VectSub(ray->orig, S->P, V0);
	double a = VectDotProd(ray->dir, ray->dir);
	double b = 2.0*VectDotProd(V0, ray->dir);
	double c = VectDotProd(V0, V0) - S->radius*S->radius;
	if ((b*b - 4 * a*c) > 0)
	{
		inter = make_intersection();
		inter->t = (-b - sqrt(b*b - 4 * a*c)) / (2.0*a);
		// set normal
		Vect Ingoing;
		VectAddS(inter->t, ray->dir, ray->orig, Ingoing);
		VectAddS(-1.0, ray->orig, Ingoing, inter->N);
		VectUnit(inter->N);

		inter->surf = S->surf;

		// this the first hit 

		if (!nearest_inter)
		{
			nearest_inter = inter;
		}

		// this is closer than any previous hit

		else if (inter->t < nearest_inter->t)
		{
			free(nearest_inter);
			nearest_inter = inter;
		}

		// something else is closer--move along

		else
		{
			free(inter);
		}
	}

	return nearest_inter;
}

//----------------------------------------------------------------------------

// only local, ambient + diffuse lighting (no specular, shadows, reflections, or refractions)

void shade_ray_diffuse(Ray *ray, Intersection *inter, Vect color)
{
	Vect L, n;
	double diff_factor;

	// iterate over lights

	for (int i = 0; i < light_list.size(); i++) {

		// AMBIENT

		color[R] += inter->surf->amb[R] * light_list[i]->amb[R];
		color[G] += inter->surf->amb[G] * light_list[i]->amb[G];
		color[B] += inter->surf->amb[B] * light_list[i]->amb[B];

		// DIFFUSE 

		// INSERT YOUR CODE HERE
		VectSub(light_list[i]->P, inter->P, L);
		VectUnit(L);
		VectUnit(inter->N);
		double NdotL = VectDotProd(inter->N, L);
		color[R] += inter->surf->diff[R] * light_list[i]->diff[R] * NdotL;
		color[G] += inter->surf->diff[G] * light_list[i]->diff[G] * NdotL;
		color[B] += inter->surf->diff[B] * light_list[i]->diff[B] * NdotL;
	}

	// clamp color to [0, 1]

	VectClamp(color, 0, 1);
}

//----------------------------------------------------------------------------

// same as shade_ray_diffuse(), but add specular lighting + shadow rays (i.e., full Phong illumination model)

void shade_ray_local(Ray *ray, Intersection *inter, Vect color)
{
	// INSERT YOUR CODE HERE (EXTRA CREDITS FOR SHADOW RAY) 
	Vect LightVec, RefectLight, ViewVec;
	double diff_factor;
	double spec_factor;

	// iterate over lights

	for (int i = 0; i < light_list.size(); i++) {

		// AMBIENT

		color[R] += inter->surf->amb[R] * light_list[i]->amb[R];
		color[G] += inter->surf->amb[G] * light_list[i]->amb[G];
		color[B] += inter->surf->amb[B] * light_list[i]->amb[B];

		// DIFFUSE 

		VectSub(light_list[i]->P, inter->P, LightVec);
		VectUnit(LightVec);
		VectUnit(inter->N);
		double NdotL = VectDotProd(inter->N, LightVec);
		color[R] += inter->surf->diff[R] * light_list[i]->diff[R] * NdotL;
		color[G] += inter->surf->diff[G] * light_list[i]->diff[G] * NdotL;
		color[B] += inter->surf->diff[B] * light_list[i]->diff[B] * NdotL;

		// SPECULAR

		Vect temp = { 0 };
		VectAddS(2 * VectDotProd(inter->N, LightVec), inter->N, temp, temp);
		VectSub(LightVec, temp, RefectLight);
		VectUnit(RefectLight);
		VectAddS(inter->t, ray->dir, ray->orig, ViewVec);
		VectUnit(ViewVec);
		double VdotR = VectDotProd(ViewVec, RefectLight);
		color[R] += inter->surf->spec[R] * light_list[i]->spec[R] * pow(VdotR, inter->surf->spec_exp);
		color[G] += inter->surf->spec[G] * light_list[i]->spec[G] * pow(VdotR, inter->surf->spec_exp);
		color[B] += inter->surf->spec[B] * light_list[i]->spec[B] * pow(VdotR, inter->surf->spec_exp);
	}

	// clamp color to [0, 1]

	VectClamp(color, 0, 1);
}

//----------------------------------------------------------------------------

// full shading model: ambient/diffuse/specular lighting, shadow rays, recursion for reflection, refraction

// level = recursion level (only used for reflection/refraction)

void shade_ray_recursive(int level, double weight, Ray *ray, Intersection *inter, Vect color)
{
	Surface *surf = (Surface *)malloc(sizeof(Surface));
	int i;

	// initialize color to Phong reflectance model

	shade_ray_local(ray, inter, color);

	// if not too deep, recurse

	if (level + 1 < maxlevel) {

		// add reflection component to color

		if (surf->reflectivity * weight > minweight) {

			// INSERT YOUR CODE HERE
			Vect InComing, OutGoing;
			VectAddS(inter->t, ray->dir, ray->orig, InComing);
			reflection_direction(InComing, inter->N, OutGoing);
			Ray *outRay = make_ray(inter->P, OutGoing);
			trace_ray(level, weight, outRay, color);
		}

		// add refraction component to color

		if (surf->transparency * weight > minweight) {

			// INSERT YOUR CODE HERE (EXTRA CREDITS FOR REFRACTION)
			double n1 = 1;
			double n2 = inter->surf->ior;
			Vect InComing, OutGoing;
			VectAddS(inter->t, ray->dir, ray->orig, InComing);
			double c1 = -VectDotProd(InComing, inter->N);
			double c2 = sqrt(1 - (n1 / n2)*(n1 / n2)*(1 - c1*c1));
			Vect temp = {0};
			VectAddS((n1 / n2)*c1 - c2, inter->N, temp, temp);
			VectAddS(n1 / n2, InComing, temp, OutGoing);
			Ray *outRay = make_ray(inter->P, OutGoing);
			trace_ray(level, weight, outRay, color);
		}
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// ray trace another pixel if the image isn't finished yet

void idle()
{
	int count = 0;
	while (count < 100)
	{
		if (image_j < ray_cam->im->h) {

			raytrace_one_pixel(image_i, image_j);

			image_i++;

			if (image_i == ray_cam->im->w) {
				image_i = 0;
				image_j++;
			}
		}
		else
		{
			myWave(model_list[0], double(count/10));
			count++;
		}
		glutPostRedisplay();
	}
	

	// write rendered image to file when done

}

//----------------------------------------------------------------------------

// show the image so far

void display(void)
{
	// draw it!

	glPixelZoom(1, -1);
	glRasterPos2i(0, ray_cam->im->h);

	glDrawPixels(ray_cam->im->w, ray_cam->im->h, GL_RGBA, GL_FLOAT, ray_cam->im->data);

	glFlush();
}

//----------------------------------------------------------------------------

void init()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, ray_cam->im->w, 0.0, ray_cam->im->h);
}

//----------------------------------------------------------------------------

void createObj(int h)
{
	FILE *fp1 = fopen("test2.obj", "wb");
	double ddd = 1.0 / (h-1.0);
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < h; j++)
		{			
			fprintf(fp1, "v %f %f %f\n", i*ddd, j*ddd, 0.0f);
		}
	}

	for (int i = 0; i < h - 1; i+=2)
	{
		for (int j = 0; j < h; j+=2)
		{
			fprintf(fp1, "f %d %d %d\n", i * h + j, i * h + j + 1, (i + 1) * h + j);
			fprintf(fp1, "f %d %d %d\n", i * h + j + 1, (i + 1) * h + j, (i + 1) * h + j + 1);
		}
	}
	fclose(fp1);
}

int main(int argc, char** argv)
{
	createObj(50);
	glutInit(&argc, argv);

	// initialize scene (must be done before scene file is parsed)
	Camera* new_cam = (Camera *)malloc(sizeof(Camera));
	ray_cam = (Camera *)malloc(sizeof(Camera));
	init_raytracing(ray_cam);

	if (argc == 2)
		parse_scene_file(argv[1], ray_cam);
	else {
		printf("missing .scene file\n");
		exit(1);
	}

	myWave(model_list[0], 1);
	for (int i = 0; i < 50; i++)
	{
		for (int j = 0; j < 50; j++)
		{
			model_list[0]->vertices[i * 50 + j];
		}
	}
	// opengl business

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ray_cam->im->w, ray_cam->im->h);
	glutInitWindowPosition(500, 300);
	glutCreateWindow("hw4_raytracing");
	init();

	glutDisplayFunc(display);
	glutIdleFunc(idle);

	glutMainLoop();

	return 0;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
