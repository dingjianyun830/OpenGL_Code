// This code is the homework1 of course EE7755
// Author: Yuqi Ding(yding18@lsu.edu)

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <stdlib.h>
#include <vector>
#include <GL/glut.h>

using namespace std;

struct VERTEX
{
	float x;
	float y;
	float z;
};

struct FACE
{
	int V1;
	int V2;
	int V3;
};

vector<VERTEX> Vertex;
vector<FACE> Face;

float x_BCenter = 0.0f;
float y_BCenter = 0.0f;
float z_BCenter = 0.0f;
float BoundingBoxDiagonalAxisLength = 0.0f;

bool ReadOBJFile(const char filename[])
{
	string name(filename);
	ifstream fobj(name);
	string str;
	string m_nLabelV("v");
	string m_nLabelF("f");
	while (getline(fobj, str))
	{
		istringstream input(str);
		string m_tS, m_nS1, m_nS2, m_nS3;
		input >> m_tS >> m_nS1 >> m_nS2 >> m_nS3;
		if (m_tS == m_nLabelV)
		{	
			// if vertex, add
			VERTEX v1;
			v1.x = stof(m_nS1);
			v1.y = stof(m_nS2);
			v1.z = stof(m_nS3);
			Vertex.push_back(v1);
			
		}
		else if (m_tS == m_nLabelF)
		{
			// if face, add
			FACE f;
			f.V1 = stoi(m_nS1)-1;
			f.V2 = stoi(m_nS2)-1;
			f.V3 = stoi(m_nS3)-1;
			Face.push_back(f);
		}
		else
		{
			cout << "Can not read the .obj file. Please check it." << endl;
			break;
		}
	}

	if ((Vertex.size() != 0) || (Face.size() != 0))
	{
		cout << "This .obj file has " << Vertex.size() << " vertexs" << endl;
		cout << "This .obj file has " << Face.size() << " faces" << endl;
	}

	return 1;
}

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
	for (int i = 0; i < Vertex.size(); i++)
	{
		if (Vertex[i].x >= xMax)
		{
			xMax = Vertex[i].x;
		}

		if (Vertex[i].x <= xMin)
		{
			xMin = Vertex[i].x;
		}

		if (Vertex[i].y >= yMax)
		{
			yMax = Vertex[i].y;
		}

		if (Vertex[i].y <= yMin)
		{
			yMin = Vertex[i].y;
		}

		if (Vertex[i].z >= zMax)
		{
			zMax = Vertex[i].z;
		}

		if (Vertex[i].z < zMin)
		{
			zMin = Vertex[i].z;
		}
	}

	x_BCenter = (xMax + xMin) / 2;
	y_BCenter = (yMax + yMin) / 2;
	z_BCenter = (zMax + zMin) / 2;
	BoundingBoxDiagonalAxisLength = ComputeDiag(xMax - xMin, yMax - yMin, zMax - zMin);

	return 1;
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
	gluPerspective(90, (double)w / (double)h, 1.0, 400.0);
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
	gluPerspective(90.0, 1.0, 1.0, 400.0);

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
	glColor3f(1.0, 0.0, 1.0);

	//Begin Mesh
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < Face.size(); i++)
	{
		glVertex3f(Vertex[Face[i].V1].x, Vertex[Face[i].V1].y, Vertex[Face[i].V1].z);

		glVertex3f(Vertex[Face[i].V2].x, Vertex[Face[i].V2].y, Vertex[Face[i].V2].z);

		glVertex3f(Vertex[Face[i].V3].x, Vertex[Face[i].V3].y, Vertex[Face[i].V3].z);	
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
	//char ObjFileName[] = "C:\\Users\\robin\\Desktop\\EE7755\\OBJ Meshes\\bunny.obj";
	// read the obj file
	ReadOBJFile(argv[1]);
	//ReadOBJFile(ObjFileName);
	cout << "Compute the Bounding Box..." << endl;
	ComputeBoundingBox();
	cout << "The Center(" << x_BCenter << "," << y_BCenter << "," << z_BCenter << ")" << endl;
	cout << "Set the camera position..." << endl;
	x_cam = x_BCenter;
	y_cam = y_BCenter;
	z_cam = z_BCenter + 1.5 * BoundingBoxDiagonalAxisLength;
	cout << "The Camera(" << x_cam << "," << y_cam << "," << y_cam << ")" << endl;

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