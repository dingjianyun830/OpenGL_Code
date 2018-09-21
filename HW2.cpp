// This code is the homework1 of course EE7755
// Author: Yuqi Ding(yding18@lsu.edu)

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <stdlib.h>
#include <vector>
#include <GL/glut.h>
#include "Mesh.h"
#include "Iterators.h"

using namespace std;

Mesh *cMesh = new Mesh();
vector<Point> cNormal;

float x_BCenter = 0.0f;
float y_BCenter = 0.0f;
float z_BCenter = 0.0f;
float BoundingBoxDiagonalAxisLength = 0.0f;

float ComputeDiag(float x, float y, float z)
{
	return sqrtf(x*x + y*y + z*z);
}

bool ComputeBoundingBox()
{
	float xMax = -9999.0f;
	float yMax = -9999.0f;
	float zMax = -9999.0f;
	float xMin = 9999.0f;
	float yMin = 9999.0f;
	float zMin = 9999.0f;
	float sumX = 0.0f;
	float sumY = 0.0f;
	float sumZ = 0.0f;
	for (MeshVertexIterator vit(cMesh); !vit.end(); ++vit)
	{
		Vertex *ver = *vit;
		
		if (ver->point()[0] >= xMax)
		{
			xMax = ver->point()[0];
		}

		if (ver->point()[0] <= xMin)
		{
			xMin = ver->point()[0];
		}

		if (ver->point()[1] >= yMax)
		{
			yMax = ver->point()[1];
		}

		if (ver->point()[1] <= yMin)
		{
			yMin = ver->point()[1];
		}

		if (ver->point()[2] >= zMax)
		{
			zMax = ver->point()[2];
		}

		if (ver->point()[2] < zMin)
		{
			zMin = ver->point()[2];
		}

		sumX += ver->point()[0];
		sumY += ver->point()[1];
		sumZ += ver->point()[2];
	}


	x_BCenter = sumX / cMesh->numVertices();
	y_BCenter = sumY / cMesh->numVertices();
	z_BCenter = sumZ / cMesh->numVertices();
	BoundingBoxDiagonalAxisLength = ComputeDiag(xMax - xMin, yMax - yMin, zMax - zMin);

	return 1;
}

double ComputeFaceArea(Face * f)
{
	Halfedge * he1 = f->he(); // one halfedge
	Halfedge * he2 = he1->next(); // its next halfedge inside this face
	Point & pt1 = he1->source()->point();
	Point & pt2 = he1->target()->point();
	Point & pt3 = he2->target()->point();
	Point cprod = (pt2 - pt1) ^ (pt3 - pt1);
	return cprod.norm() / 2.0;
}

void ComputeFaceCornerAngles(Face * f, double cAngles[3])
{
	Halfedge * hes[3];
	hes[0] = f->he();
	hes[1] = hes[0]->next();
	hes[2] = hes[0]->prev();
	Point & p0 = hes[0]->target()->point();
	Point & p1 = hes[1]->target()->point();
	Point & p2 = hes[2]->target()->point();
	double l20 = (p0 - p2).norm();
	double l01 = (p1 - p0).norm();
	double l12 = (p2 - p1).norm();
	cAngles[0] = acos((l20*l20 + l01*l01 - l12*l12) / (2 * l20*l01));
	cAngles[1] = acos((l01*l01 + l12*l12 - l20*l20) / (2 * l01*l12));
	cAngles[2] = acos((l12*l12 + l20*l20 - l01*l01) / (2 * l12*l20));
}

void ComputerFaceNormal()
{
	for (MeshFaceIterator fit(cMesh); !fit.end(); ++fit)
	{
		Face *f = *fit;
		Halfedge * he1 = f->he(); // one halfedge
		Halfedge * he2 = he1->next(); // its next halfedge inside this face
		Point & pt1 = he1->source()->point();
		Point & pt2 = he1->target()->point();
		Point & pt3 = he2->target()->point();
		cNormal.push_back((pt2 - pt1) ^ (pt3 - pt1));
	}
}
/////////////////////////////////////////////////////////////////////////////////
// rendering code
char keyControl;

int mouseButton;

float obj_angle_x = 0.0f;
float obj_angle_y = 0.0f;

float obj_trans[2] = { 0.0f, 0.0f };

float camera_zoom = 1.0f;

float x_cam = 0.0f;
float y_cam = 0.0f;
float z_cam = 0.0f;

float scale = 1.0f;

/* Some variables to measure mouse movement	*/
int mousePositionX0 = 0, mousePositionY0 = 0;

void initRendering()
{
	glClearColor(0, 0, 0, 1);
	keyControl = 0;
}

//Called when a key is pressed
void handleKeypress(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'r':
		keyControl = 'r';// rotation
		break;
	case 't':
		keyControl = 't';// translate
		break;
	case 'z':
		keyControl = 'z';//zoom
		break;
	case 27: //Escape key
		exit(0);
	}
}

//Called when the window is resized
void handleResize(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90, (double)w / (double)h, 1.0, 200.0);
}

void mouseOperation(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		mouseButton = GLUT_LEFT_BUTTON;
	}
	else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
	{
		mouseButton = GLUT_MIDDLE_BUTTON;
	}
	else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
	{
		mouseButton = GLUT_RIGHT_BUTTON;
	}

	mousePositionX0 = x;
	mousePositionY0 = y;
}

void mouseMove(int x, int y)
{
	float rFactor = 0.2f;
	float tFactor = 0.02f;
	float zFactor = 0.01f;

	if (mouseButton == GLUT_LEFT_BUTTON)
	{
		if (keyControl == 'r')
		{
			obj_angle_y += rFactor*(x - mousePositionX0);
			obj_angle_x += rFactor*(y - mousePositionY0);
		}
		else if (keyControl == 't')
		{
			obj_trans[0] += tFactor * (x - mousePositionX0);
			obj_trans[1] += tFactor * (mousePositionY0 - y);
		}
		else if (keyControl == 'z')
		{
			camera_zoom += zFactor * (y - mousePositionY0);
		}
	}

	if (mouseButton == GLUT_MIDDLE_BUTTON)
	{
		obj_trans[0] += tFactor * (x - mousePositionX0);
		obj_trans[1] += tFactor * (mousePositionY0 - y);
	}

	if (mouseButton == GLUT_RIGHT_BUTTON)
	{
		camera_zoom += zFactor * (y - mousePositionY0);
	}

	mousePositionX0 = x;
	mousePositionY0 = y;

	// post the current window
	glutPostRedisplay();
}

void Render_Mesh()
{
	// clean the window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0, 1.0, 1.0, 200.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// set the camera position
	gluLookAt(x_cam, y_cam, z_cam, x_BCenter, y_BCenter, z_BCenter, 0, 1, 0);
	glTranslatef(x_BCenter, y_BCenter, z_BCenter);

	// Set the T,R and S Matrix.
	glTranslatef(obj_trans[0], obj_trans[1], 0);
	glRotatef(obj_angle_x, 1, 0, 0);
	glRotatef(obj_angle_y, 0, 1, 0);
	glScalef(camera_zoom, camera_zoom, camera_zoom);

	// store the state of this object 
	glPushMatrix();

	// move the obj to the original point
	glTranslatef(-x_BCenter, -y_BCenter, -z_BCenter);
	glScalef(scale, scale, scale);

	//Set the color of the object

	//Begin Mesh
	glBegin(GL_TRIANGLES);
	for (MeshFaceIterator fit(cMesh); !fit.end(); ++fit)
	{
		Face * f = *fit;
		Halfedge *he1 = f->he();
		Halfedge *he2 = he1->next();
		Point &pt1 = he1->source()->point();
		Point &pt2 = he1->target()->point();
		Point &pt3 = he2->target()->point();
		glNormal3f(cNormal[f->index()].v[0], cNormal[f->index()].v[1], cNormal[f->index()].v[2]);
		glColor3f(1, 0, 0);
		glVertex3f(pt1.v[0], pt1.v[1], pt1[2]);
		glColor3f(0.9, 0, 0);
		glVertex3f(pt2.v[0], pt2.v[1], pt2[2]);
		glColor3f(0.8, 0, 0);
		glVertex3f(pt3.v[0], pt3.v[1], pt3[2]);
	}
	glEnd();

	glPopMatrix();

	/*
	// Original Coordinate System
	glPushMatrix();
	glBegin(GL_LINES);
	// x axis red
	glColor3f(1, 0, 0);
	glVertex3f(0, 0, 0);
	glVertex3f(x_cam, 0, 0);
	// y axis green
	glColor3f(0, 1, 0);
	glVertex3f(0, 0, 0);
	glVertex3f(0, y_cam, 0);
	// z axis blue
	glColor3f(0, 0, 1);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, z_cam);
	glEnd();
	glPopMatrix();
	*/
	glutSwapBuffers();
}

int main(int argc, char** argv)
{
	const char ObjFileName[128] = "C:\\Users\\yding18\\Desktop\\E7755\\OBJ Meshes\\bunny.obj";
	// read the obj file
	bool flag = cMesh->readOBJFile(ObjFileName);

	if (!flag)
	{
		cerr << "Fail to read the obj file " << ObjFileName << endl;
		system("pause");
		return -1;
	}

	cout << "The total number of the Faces: " << cMesh->numFaces() << endl;
	cout << "The total number of the Vertexs: " << cMesh->numVertices() << endl;

	cout << "Compute the Bounding Box..." << endl;
	ComputeBoundingBox();
	cout << "The Center(" << x_BCenter << "," << y_BCenter << "," << z_BCenter << ")" << endl;
	cout << "Set the camera position..." << endl;
	x_cam = x_BCenter;
	y_cam = y_BCenter;
	z_cam = z_BCenter + 1.5 * BoundingBoxDiagonalAxisLength;
	cout << "The Camera(" << x_cam << "," << y_cam << "," << y_cam << ")" << endl;

	cout << "Compute the Normal Vector ..." << endl;
	ComputerFaceNormal();
	// Set the OpenGL for 3D rendering
	//Initialize GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(400, 400);
	glutInitWindowPosition(100, 100);

	//Create the window
	//glutCreateWindow(ObjFileName);
	glutCreateWindow(argv[1]);

	initRendering();

	// call the display function. The object will be render.
	glutDisplayFunc(Render_Mesh);

	// enable the reshape window function.
	glutReshapeFunc(handleResize);

	// enable the mouse fuction
	glutMouseFunc(mouseOperation);

	// enable the mouse motion function
	glutMotionFunc(mouseMove);

	// enable the keyboard function
	glutKeyboardFunc(handleKeypress);

	glutMainLoop();

	return 0;
}