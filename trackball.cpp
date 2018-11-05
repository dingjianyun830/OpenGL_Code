#include "stdafx.h"
#include <iostream>
#include <fstream>
using namespace std;

#include <stdlib.h>
#include <math.h>

#include <sys/types.h>
#include </freeglut-3.0.0/include/GL/glut.h>
#include </freeglut-3.0.0/include/GL/freeglut.h>

#include "objLoader.h"


#define KEY_LEFT 100
#define KEY_UP 101
#define KEY_RIGHT 102
#define KEY_DOWN 103

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

int winWidth = 1024;
int winHeight = 1024;
bool firstTime = true;

WavefrontObj *obj_data;

// Trackball parameters initialization 
float angle = 0.0, axis[3], trans[3];

bool trackingMouse = false;
bool redrawContinue = false;
bool trackballMove = false;

float lastPos[3] = { 0.0, 0.0, 0.0 };
float newPos[3] = { 0.0, 0.0, 0.0 };
float curx, cury;
float startX, startY;

// Translation & Rotation
float x_trans = 0.0f; // translate object in x direction
float y_trans = 0.0f; // translate object in y direction
float zoom = 1.0f; // zoom for scaling

//
// GLUT keypress function
// 
void Specialkey(int key, int x, int y)
{
	// Write your own code below (Hint: write response for keypress for up/down/left/right arrow, which has been #define as KEY_UP/KEY_DOWN/KEY_LEFT/KEY_RIGHT)
	// This function is use up/left/down/right to implement the translation.
	switch (key)
	{
	case KEY_UP:
		y_trans += 1.0f;
		break;
	case KEY_DOWN:
		y_trans -= 1.0f;
		break;
	case KEY_LEFT:
		x_trans -= 1.0f;
		break;
	case KEY_RIGHT:
		x_trans += 1.0f;
		break;
	}
	// Write your own code above
   
    glutPostRedisplay();
}

void onKeyPress(unsigned char key, int x, int y)
{

	if (key == 'p')
	{
		obj_data->mode = GL_LINE_LOOP;
	}
	else if (key == 's')
	{
		glShadeModel(GL_SMOOTH);								// Set Smooth Shading 
		obj_data->mode = GL_POLYGON;
	}
	else if (key == 'f')
	{
		glShadeModel(GL_FLAT);								// Set Smooth Shading 
		obj_data->mode = GL_POLYGON;
	}
	else if (key == 'q')
	{
		delete obj_data;
		exit(0);
	}
	else if (key == 'c')
	{
		redrawContinue = true;
	}

	glutPostRedisplay();
}

void MouseWheel(int wheel, int direction, int x, int y)
{
	// Write your own code below (Hint: set zoom for mouse wheel)
	// This function uses mouse wheel to implement the Zoom.
	wheel = 0;
	if (direction >0)
	{
		zoom += 0.3f;
	}
	else
	{
		zoom -= 0.3f;
	}

	// Write your own code above	
	glutPostRedisplay();

}

float length(float *pos)
{
	return sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
}

float ComputeAngle(float* pos1, float* pos2)
{
	return acos(pos1[0] * pos2[0] + pos1[1] * pos2[1] + pos1[2] * pos2[2])*180.0 / M_PI;
}

void Normalize(float * pos)
{
	float d = length(pos);
	pos[0] = pos[0] / d;
	pos[1] = pos[1] / d;
	pos[2] = pos[2] / d;
}

void ComputeRotate(float* pos1, float* pos2)
{
	Normalize(pos1);
	Normalize(pos2);
	angle += ComputeAngle(pos1, pos2);

	axis[0] += pos1[1] * pos2[2] - pos1[2] * pos2[1];
	axis[1] += pos1[2] * pos2[0] - pos1[0] * pos2[2];
	axis[2] += pos1[0] * pos2[1] - pos1[1] * pos2[0];
	Normalize(axis);
}

void trackBallMapping(int x, int y, float *vec3)
{
	float d = 0;
	vec3[0] = (2.0 * x - winWidth) / winWidth;
	vec3[1] = (winHeight - 2.0 * y) / winHeight;
	vec3[2] = sqrt(1.0 - vec3[0]*vec3[0] - vec3[1]*vec3[1]);
}

void mouseMotion(int x, int y)
{
	// Write your own code below
	if (trackingMouse)
	{
		trackBallMapping(x, y, newPos);

		ComputeRotate(lastPos, newPos);
	}
	else if (trackballMove)
	{
		curx = (x - winWidth / 2.0f) / winWidth;
		cury = (y - winHeight / 2.0f) / winHeight;
		trans[0] += curx - startX;
		trans[1] += startY - cury;
		trans[2] += 0;
	}
	glutPostRedisplay();
	// Write your own code above
}

void mouseButton(int button, int state, int x, int y)
{
	// Write your own code below (Hint:holding down left button allows user to rotate cube)
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		trackBallMapping(x, y, lastPos);
		trackingMouse = true;
	}
	else if(button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
	{
		startX = (x - winWidth / 2.0f) / winWidth;
		startY = (y - winHeight / 2.0f) / winHeight;
		trackballMove = true;
	}
	else
	{
		trackingMouse = false;
		trackballMove = false;
	}
	// Write your own code above
}

void Init(int w, int h)
{
	glViewport(0, 0, w, h);
	glShadeModel(GL_SMOOTH);								// Set Smooth Shading 
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);			    	// Background Color 
	glClearDepth(1.0f);										// Depth buffer setup 
	glEnable(GL_DEPTH_TEST);								// Enables Depth Testing 
	glDepthFunc(GL_LEQUAL);									// The Type Of Depth Test To Do 
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);		// Use perspective correct interpolation if available

	glMatrixMode(GL_PROJECTION);							// Select The Projection Matrix
	glLoadIdentity();										// Reset The Projection Matrix
	double aspect = (double)h / w;
	glFrustum(-5, 5, -5 * aspect, 5 * aspect, 10, 500);          // Define perspective projection frustum
																 //gluPerspective(30, w/h, 10, 74);
	glTranslated(0.0, 0.0, -24);                          // Viewing transformation

	glMatrixMode(GL_MODELVIEW);								// Select The Modelview Matrix
	glLoadIdentity();										// Reset The Modelview Matrix

	if (firstTime)
	{
		glEnable(GL_LIGHTING);
		float frontColor[] = { 0.2f, 0.7f, 0.7f, 1.0f };

		glMaterialfv(GL_FRONT, GL_AMBIENT, frontColor);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, frontColor);
		glMaterialfv(GL_FRONT, GL_SPECULAR, frontColor);
		glMaterialf(GL_FRONT, GL_SHININESS, 100);

		float lightDirection[] = { 2.0f, 0.0f, 1.0f, 0 };
		float ambientIntensity[] = { 0.1f, 0.1f, 0.1f, 1.0f };
		float lightIntensity[] = { 0.9f, 0.9f, 0.9f, 1.0f };
		float lightSpec[] = { 1.0f, 1.0f, 1.0f, 1 };
		glLightfv(GL_LIGHT0, GL_AMBIENT, ambientIntensity);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, lightIntensity);
		glLightfv(GL_LIGHT0, GL_POSITION, lightDirection);
		glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec);
		glEnable(GL_LIGHT0);
		firstTime = false;
	}
}

void Draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		// Clear The Screen And The Depth Buffer

	// Write your own code below (Hint:Translation Rotation & Scaling)

	glPushMatrix();
	glTranslatef(x_trans, y_trans, 0);
	glRotatef(angle, axis[0], axis[1], axis[2]);
	glScalef(zoom, zoom, zoom);

	// Write your own code above

	if (obj_data != NULL)
		obj_data->Draw();
	else
		glutSolidTeapot(1.0);	//draw a teapot when no argument is provided
	glPopMatrix();

	// display a trackball
	glPushMatrix();
	// render the local min
	glColor3f(1.0f, 1.0f, 0.0f);
	glRotatef(angle, axis[0], axis[1], axis[2]);
	glutWireCube(10.0f);
	glPopMatrix();

	// Write your own code below (Hint:reset translate and scale value)
	if (redrawContinue)
	{
		zoom = 1;
		x_trans = 0;
		y_trans = 0;
		redrawContinue = false;
	}
	
	// Write your own code above
	glutSwapBuffers();
}


int main( int argc, char *argv[] ) 
{  

    // glut initialization functions:
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
    glutInitWindowSize(winWidth, winHeight);
    glutInitWindowPosition(100,100);
    glutCreateWindow("ImageViewer");

    Init(winWidth, winHeight);

    // display, onMouseButton, mouse_motion, onKeyPress, and resize are functions defined above
    glutDisplayFunc(Draw);       
    glutKeyboardFunc(onKeyPress);
	glutSpecialFunc(Specialkey);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMotion);
	glutMouseWheelFunc(MouseWheel);
    glutReshapeFunc(Init);

    if ( argc >= 2 ) 
		obj_data = new WavefrontObj( argv[1] );
    else
		obj_data = NULL;

    // start glutMainLoop -- infinite loop
    glutMainLoop();

    return 0;
}
