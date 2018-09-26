// This code is the homework2 of course EE7755
// Author: Yuqi Ding(yding18@lsu.edu)

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <stdlib.h>
#include <vector>
#include <GL\glut.h>
#include "Mesh.h"
#include "Iterators.h"

using namespace std;

Mesh *cMesh = new Mesh();
vector<Point> cNormalFace;
vector<Point> cNormalVertex;
vector<double> cGaussCurv;
vector<int> LocalMaxGC;
vector<int> LocalMinGC;
vector<int> ShapeEdge;
#define PI 3.1415926

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
		//cout << f->index() << endl;
		Halfedge * he1 = f->he(); // one halfedge
		Halfedge * he2 = he1->next(); // its next halfedge inside this face
		Point & pt1 = he1->source()->point();
		Point & pt2 = he1->target()->point();
		Point & pt3 = he2->target()->point();
		Point pt = (pt2 - pt1) ^ (pt3 - pt1);
		pt = pt / pt.norm();
		// add the face normal vector
		cNormalFace.push_back(pt);
	}
}

void ComputeVertexNormal()
{
	for (MeshVertexIterator vit(cMesh); !vit.end(); ++vit)
	{
		Vertex *ver = *vit;
		//cout<<ver->index()<<endl;
		Point cNormalV(0, 0, 0);
		double sumAngles = 0;
		for (VertexOutHalfedgeIterator vheit(ver); !vheit.end(); ++vheit)
		{
			// select the two halfedges
			Halfedge *he = *vheit;
			Halfedge *he1 = he->prev();

			// select 3 points, p1 is the current vertex
			Point & p0 = he1->source()->point();
			Point & p1 = he->source()->point();
			Point & p2 = he->target()->point();

			// compute the incident angle as weight
			double a = (p0 - p1).norm();
			double b = (p2 - p1).norm();
			double cAngles = acos(((p0 - p1)*(p2 - p1)) / (a * b));
			sumAngles += cAngles;

			// select the face normal vector
			int FaceIndex = he->face()->index();
			Point pt1 = cNormalFace[FaceIndex];
			pt1 *= cAngles;
			cNormalV += pt1;
		}
		cNormalV = cNormalV / sumAngles;
		cNormalV = cNormalV / cNormalV.norm();
		// add the vertex normal vector
		cNormalVertex.push_back(cNormalV);
	}
}

void ComputeGaussianCurv()
{
	for (MeshVertexIterator vit(cMesh); !vit.end(); ++vit)
	{
		Vertex *ver = *vit;
		//cout<<ver->index()<<endl;
		double GuassCurv = 0;
		double sumAngles = 0;
		for (VertexOutHalfedgeIterator vheit(ver); !vheit.end(); ++vheit)
		{
			// select the two halfedges
			Halfedge *he = *vheit;
			Halfedge *he1 = he->prev();

			// select 3 points, p1 is the current vertex
			Point & p0 = he1->source()->point();
			Point & p1 = he->source()->point();
			Point & p2 = he->target()->point();

			// compute the incident angle as weight
			double a = (p0 - p1).norm();
			double b = (p2 - p1).norm();
			double cAngles = acos(((p0 - p1)*(p2 - p1)) / (a * b));
			sumAngles += cAngles;
		}
		GuassCurv = 2 * PI - sumAngles;
		// add the vertex normal vector
		cGaussCurv.push_back(GuassCurv);
	}
}

void LocalMaxMinGC()
{
	for (MeshVertexIterator vit(cMesh); !vit.end(); ++vit)
	{
		Vertex *ver1 = *vit;
		int iMax = ver1->index();
		int iMin = ver1->index();
		double vMax = cGaussCurv[iMax];
		double vMin = cGaussCurv[iMin];
		vector<int> Index;
		for (VertexVertexIterator vvit(ver1); !vvit.end(); ++vvit)
		{
			Vertex *ver2 = *vvit;
			for (VertexVertexIterator vvit2(ver2); !vvit2.end(); ++vvit2)
			{
				Vertex *ver = *vvit2;
				Index.push_back(ver->index());
			}
		}

		// find the local Max and Min
		for (int i = 0; i < Index.size(); i++)
		{
			if (cGaussCurv[Index[i]] >= vMax)
			{
				vMax = cGaussCurv[Index[i]];
				iMax = Index[i];
			}

			if (cGaussCurv[Index[i]] <= vMin)
			{
				vMin = cGaussCurv[Index[i]];
				iMin = Index[i];
			}
		}

		LocalMaxGC.push_back(iMax);
		LocalMinGC.push_back(iMin);
	}
}

void ComputeShapeEdge()
{
	int count = 0;
	for (MeshEdgeIterator eit(cMesh); !eit.end(); ++eit)
	{
		Edge *e = *eit;
		if (e->boundary())
		{
			ShapeEdge.push_back(0);
		}
		else
		{
			Halfedge *he0 = e->he(0);
			Halfedge *he1 = e->he(1);
			//cout << he0->face()->index() << endl;
			//cout << he1->face()->index() << endl;
			Point & pt0 = cNormalFace[he0->face()->index()];
			Point & pt1 = cNormalFace[he1->face()->index()];
			if (0.5 > pt0*pt1)
			{
				ShapeEdge.push_back(1);
				count++;
			}
			else
			{
				ShapeEdge.push_back(0);
			}
		}
	}
	cout << count << endl;
}
/////////////////////////////////////////////////////////////////////////////////
// rendering code
char keyControl; // mouse control
char keyControl1; // normal vector control 
char keyControl2; // Gassian Curvature
char keyControl3; // ShapeEdge

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
	//Define Material Properties for the Objects in the Scene
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 100.0 };

	//The following color components specify the intensity for each type of lights.
	GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };

	GLfloat light_position[] = { 10.0, 10.0, 10.0, 0.0 };

	glClearColor(0, 0, 0, 1);
	keyControl = 0;
	keyControl1 = 0;
	keyControl2 = 0;
	keyControl3 = 0;
	glShadeModel(GL_SMOOTH);

	//Material properties determine how it reflects light
	//You can specify a material's ambient, diffuse, and specular colors and how shiny it is.
	//Here only the last two material properties are explicitly specified
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);


	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	//glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	//Specify the LIGHT0 position to be "light_position"
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);


	//This enables lighting calculations
	glEnable(GL_LIGHTING);

	//Remember to enable the light you just defined 
	glEnable(GL_LIGHT0);

	glEnable(GL_COLOR_MATERIAL);

	glEnable(GL_DEPTH_TEST);
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
	case 'f':
		keyControl1 = 'f';//face normal
		break;
	case 'v':
		keyControl1 = 'v';//vertex normal
		break;
	case 'k':
		keyControl2 = 'k';//Gaussian Curvature
		break;
	case 's':
		keyControl3 = 's';//Shape Edge
		break;
	case 'c':
		keyControl2 = 'c';
		keyControl3 = 'c';
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

//Begin Mesh
void Render_Mesh()
{
	glBegin(GL_TRIANGLES);
	for (MeshFaceIterator fit(cMesh); !fit.end(); ++fit)
	{
		Face * f = *fit;
		int ind = f->index();
		Halfedge *he1 = f->he();
		Halfedge *he2 = he1->next();
		Point &pt1 = he1->source()->point();
		Point &pt2 = he1->target()->point();
		Point &pt3 = he2->target()->point();
		int ind1 = he1->source()->index();
		int ind2 = he1->target()->index();
		int ind3 = he2->target()->index();

		// point 1
		glColor3f(1, 1, 1);
		glVertex3f(pt1.v[0], pt1.v[1], pt1[2]);

		// point 2
		glColor3f(1, 1, 1);
		glVertex3f(pt2.v[0], pt2.v[1], pt2[2]);

		// point 3
		glColor3f(1, 1, 1);
		glVertex3f(pt3.v[0], pt3.v[1], pt3[2]);
	}
	glEnd();
}

void Render_Mesh_f()
{
	glBegin(GL_TRIANGLES);
	for (MeshFaceIterator fit(cMesh); !fit.end(); ++fit)
	{
		Face * f = *fit;
		int ind = f->index();
		Halfedge *he1 = f->he();
		Halfedge *he2 = he1->next();
		Point &pt1 = he1->source()->point();
		Point &pt2 = he1->target()->point();
		Point &pt3 = he2->target()->point();
		int ind1 = he1->source()->index();
		int ind2 = he1->target()->index();
		int ind3 = he2->target()->index();
		Point p1 = cNormalFace[ind];
		Point p2 = cNormalFace[ind];
		Point p3 = cNormalFace[ind];

		// point 1
		glColor3f(1, 1, 1);
		glNormal3f(p1.v[0], p1.v[1], p1.v[2]);
		glVertex3f(pt1.v[0], pt1.v[1], pt1[2]);

		// point 2
		glColor3f(1, 1, 1);
		glNormal3f(p2.v[0], p2.v[1], p2.v[2]);
		glVertex3f(pt2.v[0], pt2.v[1], pt2[2]);

		// point 3
		glColor3f(1, 1, 1);
		glNormal3f(p3.v[0], p3.v[1], p3.v[2]);
		glVertex3f(pt3.v[0], pt3.v[1], pt3[2]);
	}
	glEnd();
}

void Render_Mesh_v()
{
	glBegin(GL_TRIANGLES);
	for (MeshFaceIterator fit(cMesh); !fit.end(); ++fit)
	{
		Face * f = *fit;
		int ind = f->index();
		Halfedge *he1 = f->he();
		Halfedge *he2 = he1->next();
		Point &pt1 = he1->source()->point();
		Point &pt2 = he1->target()->point();
		Point &pt3 = he2->target()->point();
		int ind1 = he1->source()->index();
		int ind2 = he1->target()->index();
		int ind3 = he2->target()->index();
		Point p1 = cNormalVertex[ind1];
		Point p2 = cNormalVertex[ind2];
		Point p3 = cNormalVertex[ind3];

		// point 1
		glColor3f(1, 1, 1);
		glNormal3f(p1.v[0], p1.v[1], p1.v[2]);
		glVertex3f(pt1.v[0], pt1.v[1], pt1[2]);

		// point 2
		glColor3f(1, 1, 1);
		glNormal3f(p2.v[0], p2.v[1], p2.v[2]);
		glVertex3f(pt2.v[0], pt2.v[1], pt2[2]);

		// point 3
		glColor3f(1, 1, 1);
		glNormal3f(p3.v[0], p3.v[1], p3.v[2]);
		glVertex3f(pt3.v[0], pt3.v[1], pt3[2]);
	}
	glEnd();
}

void Render_GaussCurv()
{
	for (MeshVertexIterator vit(cMesh); !vit.end(); ++vit)
	{
		Vertex *ver = *vit;
		int indMax = LocalMaxGC[ver->index()];
		int indMin = LocalMinGC[ver->index()];

		Vertex *vMax = cMesh->indVertex(indMax);
		Vertex *vMin = cMesh->indVertex(indMin);

		glPushMatrix();
		// render the local max
		glColor3f(1.0f, 0.0f, 0.0f);
		glTranslatef(vMax->point().v[0], vMax->point().v[1], vMax->point().v[2]);
		glutSolidSphere(0.005, 15, 15);
		glPopMatrix();

		glPushMatrix();
		// render the local min
		glColor3f(0.0f, 0.0f, 1.0f);
		glTranslatef(vMin->point().v[0], vMin->point().v[1], vMin->point().v[2]);
		glutSolidSphere(0.005, 15, 15);
		glPopMatrix();
	}
}

void Render_ShapeEdge()
{
	glBegin(GL_LINES);
	for (MeshEdgeIterator eit(cMesh); !eit.end(); ++eit)
	{
		Edge *e = *eit;
		if (ShapeEdge[e->index()] == 1)
		{
			Halfedge *he = e->he(0);
			Point pt0 = he->source()->point();
			Point pt1 = he->target()->point();
			glColor3f(0, 0, 1);
			glVertex3f(pt0.v[0], pt0.v[1], pt0.v[2]);
			glVertex3f(pt1.v[0], pt1.v[1], pt1.v[2]);
		}
	}
	glEnd();
}

void display()
{
	// clean the window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(x_cam, y_cam, z_cam, x_BCenter, y_BCenter, z_BCenter, 0, 1, 0);

	// Draw the mesh
	glPushMatrix();
	glTranslatef(+x_BCenter, +y_BCenter, +z_BCenter);
	glRotatef(obj_angle_x, 1, 0, 0);
	glRotatef(obj_angle_y, 0, 1, 0);
	glTranslatef(obj_trans[0], obj_trans[1], 0);
	glScalef(camera_zoom, camera_zoom, camera_zoom);
	glTranslatef(-x_BCenter, -y_BCenter, -z_BCenter);
	if (keyControl1 == 'f')
	{
		Render_Mesh_f();
	}
	else if (keyControl1 == 'v')
	{
		Render_Mesh_v();
	}
	else
	{
		Render_Mesh();
	}
	
	if (keyControl2 == 'k')
	{
		Render_GaussCurv();
	}
		
	if (keyControl3 == 's')
	{
		Render_ShapeEdge();
	}
	
	glPopMatrix();

	glutSwapBuffers();
}

int main(int argc, char** argv)
{
	const char ObjFileName[128] = "C:\\Users\\robin\\Desktop\\EE7755\\OBJ Meshes\\david.obj";
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

	cout << "Compute the Face Normal Vector ..." << endl;
	ComputerFaceNormal();
	cout << "Compute the Vertex Normal Vector ..." << endl;
	ComputeVertexNormal();

	cout << "Compute the Gaussian Curvature ..." << endl;
	ComputeGaussianCurv();
	cout << "Find the local feature ..." << endl;
	LocalMaxMinGC();

	cout << "Find the shape edge ..." << endl;
	ComputeShapeEdge();
	// Set the OpenGL for 3D rendering
	//Initialize GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(400, 400);
	glutInitWindowPosition(100, 100);

	//Create the window
	glutCreateWindow(ObjFileName);
	//glutCreateWindow(argv[1]);

	initRendering();

	// call the display function. The object will be render.
	glutDisplayFunc(display);

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