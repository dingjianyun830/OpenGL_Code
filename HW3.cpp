#include <iostream>
#include <math.h>
#include <GL/glut.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>


#include "Mesh.h"
#include "Iterators.h"
#include "Tex_CheckerBoard.h"
#include "CG_SPL\cg_spl.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

float whratio;
int win_height, win_width;
GLUquadricObj *obj;
CheckerBoard *CBTex = new CheckerBoard();// generate texture

/* Some variables to measure mouse movement	*/
int mousePositionX0 = 0, mousePositionY0 = 0;
int mouseButton = 0;

/* Some variables to describe the object dimension and camera placement	*/
float cameraPosition[3] = {0, 0, 2};	//default camera position
float objCenter[3] = {0, 0, 0};			//default object center
float boxMin[3] = { 0, 0, 0 };
float boxMax[3] = { 0, 0, 0 };
float axislen = 1.414;

/* Some variables to control interactive transformations	*/
float my_Rotation_x = 0, my_Rotation_y = 0;
float my_Translation[3] = { 0, 0, 0 };

int useNormal;
bool showSharpEdge;
bool showK_g;
bool showISO_u;
bool showISO_v;

void MyInit(void);

void ComputeBoundingBox();

//rendering
void SetCamera(void);
void Render_BoxAndAxes(void);
void Render_Mesh(void);
void Render_MinMaxK(void);
void Render_SharpEdge(void);

//glut functions
void mouseMove(int x, int y);
void mouseClick(int button , int state, int x, int y);
void reshape(int w, int h);

Mesh * myMesh;
std::vector<Point> fNList;
std::vector<double> cAngles;
std::vector<Point> vNList;
std::vector<double> kg;
std::vector<int> maxmin_kg;
std::vector<bool> sharpEdges;
std::vector<double> HSF; // harmonic scalar field
std::vector<Halfedge *> boundary;
std::vector<Halfedge *> boundaryLoop;
std::vector<double> arcLen;
std::vector<double> tex_U;
std::vector<double> tex_V;

void ComputeSharpEdges(){
	sharpEdges.resize(myMesh->numEdges(), false);
	for (MeshEdgeIterator eit(myMesh); !eit.end(); ++eit){
		Edge * e = *eit;
		if (e->boundary()) continue;
		Halfedge * he1 = e->he(0);
		Halfedge * he2 = e->he(1);
		Point & pFNorm1 = fNList[he1->face()->index()];
		Point & pFNorm2 = fNList[he2->face()->index()];
		double cdot = pFNorm1*pFNorm2;
		if (cdot < 0.5)
			sharpEdges[e->index()] = true;		
	}
}

void ComputeFaceNormals(){
	fNList.clear();
	for (MeshFaceIterator fit(myMesh); !fit.end(); ++fit){
		Face * f = *fit;
		Halfedge * he = f->he();
		Halfedge * nhe = he->next();
		Point p1 = he->target()->point() - he->source()->point();
		Point p2 = nhe->target()->point() - nhe->source()->point();
		Point np = p1 ^ p2;
		np /= np.norm();
		fNList.push_back(np);
	}
}

void ComputeCornerAngles(){
	cAngles.resize(myMesh->numEdges() * 2);
	for (MeshFaceIterator fit(myMesh); !fit.end(); ++fit){
		Face * f = *fit;
		Halfedge * he = f->he();
		Halfedge * nhe = he->next();
		Halfedge * phe = he->prev();
		Point & p1 = he->target()->point();
		Point & p2 = nhe->target()->point();
		Point & p3 = phe->target()->point();
		double l3 = (p2 - p1).norm();
		double l1 = (p3 - p2).norm();
		double l2 = (p1 - p3).norm();
		double theta1 = acos((l2*l2 + l3*l3 - l1*l1) / (2 * l2 * l3));
		double theta2 = acos((l1*l1 + l3*l3 - l2*l2) / (2 * l1 * l3));
		double theta3 = acos((l1*l1 + l2*l2 - l3*l3) / (2 * l1 * l2));
		cAngles[he->index()] = theta1;
		cAngles[nhe->index()] = theta2;
		cAngles[phe->index()] = theta3;
	}
}

void ComputeGaussianCurvature(){
	//double PI = 3.14159265358979;
	double twoPI = PI * 2.0;
	kg.resize(myMesh->numVertices());
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit){
		Vertex * v = *vit;
		double cK = twoPI;
		for (VertexInHalfedgeIterator hit(v); !hit.end(); ++hit){
			Halfedge * he = *hit;
			cK -= cAngles[he->index()];
		}
		kg[v->index()] = cK;
	}
}

void getMaxMinKG(){
	maxmin_kg.resize(myMesh->numVertices(), 0); //1 for max; -1 for min
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit){
		Vertex * v = *vit;
		double kv = kg[v->index()];
		if (abs(kv) < 0.1)
			continue;
		bool isMax = true;
		bool isMin = true;		
		for (VertexVertexIterator vvit(v); !vvit.end(); ++vvit){
			Vertex * vv = *vvit;
			double kvv = kg[vv->index()];
			if (kvv > kv)
				isMax = false;
			if (kvv < kv)
				isMin = false;
			if (!isMax && !isMin)
				break;
			for (VertexVertexIterator v2it(v); !v2it.end(); ++v2it){
				Vertex * v2 = *v2it;
				double kv2 = kg[v2->index()];
				if (kv2 > kv)
					isMax = false;
				if (kv2 < kv)
					isMin = false;
				if (!isMax && !isMin)
					break;
			}
			if (!isMax && !isMin)
				break;
		}
		if (isMax)
			maxmin_kg[v->index()] = 1;
		else if (isMin)
			maxmin_kg[v->index()] = -1;
	}
}

void ComputeVertNormals(){
	vNList.resize(myMesh->numVertices());
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit){
		Vertex * v = *vit;
		Point vNorm(0,0,0);
		for (VertexInHalfedgeIterator vhit(v); !vhit.end(); ++vhit){
			Halfedge * he = *vhit;
			Point & fN = fNList[he->face()->index()];
			double cAngle = cAngles[he->index()];
			vNorm += fN * cAngle;
		}
		vNList[v->index()] = vNorm / vNorm.norm();
	}
}

void IterEnergy()
{
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit)
	{
		Vertex * ver1 = *vit;
		if (ver1->index() >= 9 && ver1->index() < myMesh->numVertices() - 9)
		{
			double w = 0;
			double E1 = 0;
			// use 1-ring neighboor to pass the energy
			// use constant weight, also can use cAngles as weight
			for (VertexVertexIterator vvit(ver1); !vvit.end(); ++vvit)
			{
				Vertex *ver2 = *vvit;
				//w += cAngles[ver2->index()];	
				//E1 += cAngles[ver2->index()] * HSF[ver2->index()];
				w += 1;
				E1 += HSF[ver2->index()];
			}
			HSF[ver1->index()] = E1 / w;
		}
	}
}

double ComputeDirichletEngery()
{
	double DirichletEngery = 0;
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit)
	{
		Vertex * v = *vit;
		DirichletEngery += HSF[v->index()];
	}

	return DirichletEngery;
}

void ComputeHarmonicScalarField()
{
	HSF.resize(myMesh->numVertices(),0);
	tex_U.resize(myMesh->numVertices(), 0);
	tex_V.resize(myMesh->numVertices(), 0);

	// initialize the boundary value. Index 0-8=0, 80-88=1
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit)
	{
		Vertex * v = *vit;
		if (v->index() < 9)
		{
			HSF[v->index()] = 0;
		}
		else if(v->index() >= myMesh->numVertices()-9)
		{
			HSF[v->index()] = 1;
		}
		else
		{
			HSF[v->index()] = 0.5;
		}
	}

	// minimize the Dirichlet energy
	double oldEnergy = ComputeDirichletEngery();
	double delta = 1;
	double eplison = 1e-8;
	while (delta > eplison)
	{
		IterEnergy();
		double NewEnergy = ComputeDirichletEngery();
		delta = abs(NewEnergy - oldEnergy);
		oldEnergy = NewEnergy;
	}

	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit)
	{
		Vertex * v = *vit;
		tex_U[v->index()] = HSF[v->index()];
		tex_V[v->index()] = HSF[v->index()];
	}
}

void getBoundary()
{
	boundary.clear();
	// find the boundary edge
	for (MeshEdgeIterator eit(myMesh); !eit.end(); ++eit) 
	{
		Edge * e = *eit;
		if (e->boundary())
		{
			Halfedge * he0 = e->he(0);
			boundary.push_back(he0);	
		}
	}

	// select a random start point get the boundary loop
	// we select 0 as the start point
	Halfedge *he = boundary[0];
	for (int i = 0; i < boundary.size(); i++)
	{
		Halfedge *he1 = he->target()->most_clw_out_halfedge();
		boundaryLoop.push_back(he1);
		he = he1;
	}
}

void BoundaryPara()
{
	tex_U.assign(myMesh->numVertices(), 0);
	tex_V.assign(myMesh->numVertices(), 0);
	// compute the arclength for the boundary loop
	double aLength = 0;
	for (int i = 0; i < boundaryLoop.size(); i++)
	{
		Halfedge *he = boundaryLoop[i];
		Point &p1 = he->source()->point();
		Point &p2 = he->target()->point();
		
		aLength += (p2 - p1).norm();
		arcLen.push_back(aLength);
	}
	// use a unit disk to parameter
	for (int i = 0; i < boundaryLoop.size(); i++)
	{
		arcLen[i] = arcLen[i] / aLength;
		tex_U[boundaryLoop[i]->target()->index()] = cos(2 * PI * arcLen[i]);
		tex_V[boundaryLoop[i]->target()->index()] = sin(2 * PI * arcLen[i]);
	}
}

void HarmonicSolver()
{
	//initalize the adjmartix
	std::vector<int> rowIndex;
	std::vector<int> colIndex;
	std::vector<double> adjMatrix;
	
	//construct adjMartix => A
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit)
	{
		Vertex * ver1 = *vit;
		// for boundary, we just need to set current vertex
		if (ver1->boundary())
		{
			rowIndex.push_back(ver1->index());
			colIndex.push_back(ver1->index());
			adjMatrix.push_back(1);
			continue;
		}

		double tempSum = 0;
		for (VertexOutHalfedgeIterator vheit(ver1); !vheit.end(); ++vheit)
		{
			Halfedge *he = *vheit;
			Halfedge *he1 = he->twin()->next();
			Halfedge *he2 = he->next();

			double g1 = cAngles[he1->index()];
			double g2 = cAngles[he2->index()];
			double w = 1 / tan(g1) + 1 / tan(g2);
			tempSum += w;

			rowIndex.push_back(ver1->index());
			colIndex.push_back(he->target()->index());
			adjMatrix.push_back(-w);
			
		}
		rowIndex.push_back(ver1->index());
		colIndex.push_back(ver1->index());
		adjMatrix.push_back(tempSum);
	}

	// prepare the solver parameters
	double threshold = 1e-8;
	int num_maxIterations = 10000;
	int num_NonZeroCoef = rowIndex.size();
	int dim_x = myMesh->numVertices();
	VECTOR_double x = VECTOR_double(dim_x, 0.0);// create a vector with dimension "dim_x"
	int *row = new int[num_NonZeroCoef];    // because: num_NonZeroCoef = 6
	int *col = new int[num_NonZeroCoef];
	double *cVal = new double[num_NonZeroCoef];

	for (int i = 0; i < num_NonZeroCoef; i++)
	{
		row[i] = rowIndex[i];
		col[i] = colIndex[i];
		cVal[i] = adjMatrix[i];
	}

	double *bx = new double[dim_x]; 
	double *by = new double[dim_x];
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit)
	{
		Vertex * ver1 = *vit;
		bx[ver1->index()] = tex_U[ver1->index()];
		by[ver1->index()] = tex_V[ver1->index()];
	}

	// use the 3rd solver to solver the linear equation Ax = b
	CG_SPL theSolverX(row, col, cVal, bx, dim_x, num_NonZeroCoef, num_maxIterations, threshold);
	CG_SPL theSolverY(row, col, cVal, by, dim_x, num_NonZeroCoef, num_maxIterations, threshold);
	theSolverX.solve_GMRES("x.txt");
	theSolverY.solve_GMRES("y.txt");
	for (int i = 0; i < dim_x; ++i)
	{
		tex_U[i] = theSolverX.getx(i);
		tex_V[i] = theSolverY.getx(i);
	}		
}

void display (void)
{		
	SetCamera();
	Render_BoxAndAxes();
	Render_Mesh();
	if (showK_g)
		Render_MinMaxK();
	if (showSharpEdge)
		Render_SharpEdge();
	glutSwapBuffers();
}

void reshape(int w, int h)
{
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	win_height = h;
	win_width = w;
	glMatrixMode (GL_PROJECTION);     
	glLoadIdentity();
	whratio = (double)w / (double)h; 	//A Commonly Suggested Setting: set ratio in gluPerspective to the aspect ratio of the associated viewport
	gluPerspective(60, whratio, axislen*0.01, axislen*5);
	glMatrixMode (GL_MODELVIEW);	//change back to modelview
	glutPostRedisplay();
}

void SetCamera(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2], 
		objCenter[0], objCenter[1], objCenter[2], 0, 1, 0);

	glTranslatef(my_Translation[0], my_Translation[1], my_Translation[2]);

	glTranslatef(objCenter[0], objCenter[1], objCenter[2]);	//before doing rotation to the object, move the object center to the origin
	
	glRotatef(my_Rotation_y, 0.0, 1.0, 0.0);
	glRotatef(my_Rotation_x, 1.0, 0.0, 0.0);
			
	glTranslatef(-objCenter[0], -objCenter[1], -objCenter[2]);
}

void mouseMove(int x, int y)
{
	double movingScale = axislen / win_height;  // just a scaling factor to make the mouse moving not too sensitive
	/* rotation*/
	if (mouseButton == GLUT_LEFT_BUTTON ) 
	{
		////////////do something////////////////
		my_Rotation_y += x - mousePositionX0;
		my_Rotation_x += y - mousePositionY0;
	}

	/*xy translation */
	if (mouseButton == GLUT_MIDDLE_BUTTON)
	{
		////////////do something ////////////////
		my_Translation[0] += movingScale * (x - mousePositionX0);
		my_Translation[1] -= movingScale * (y - mousePositionY0);
	}

	/* zoom in and out */
	if (mouseButton == GLUT_RIGHT_BUTTON)
	{ // suppose we want to make moving up as zooming out
		my_Translation[2] += movingScale * (y - mousePositionY0);
	}
	mousePositionX0 = x;
	mousePositionY0 = y;
	glutPostRedisplay();
}

int count = 0;
void keyPressHandler(unsigned char key, int x, int y){	
	switch (key) {
	case 'f':
		useNormal = 1;
		glutPostRedisplay();
		break;
	case 'v':
		useNormal = 2;
		glutPostRedisplay();
		break;
	case 's':
		showSharpEdge = !showSharpEdge;
		glutPostRedisplay();
		break;
	case 'k':
		showK_g = !showK_g;
		glutPostRedisplay();
		break;
	case 't':
		count++;
		if (count % 3 == 1)
		{
			showISO_u = true;
			showISO_v = false;
		}
		else if (count % 3 == 2)
		{
			showISO_u = false;
			showISO_v = true;
		}
		else
		{
			showISO_u = false;
			showISO_v = false;
			count = 0;
		}	
		glutPostRedisplay();
		break;
	case 27: //Escape key
		exit(0);
	}
}

void mouseClick(int button , int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) 
		mouseButton = GLUT_LEFT_BUTTON;
	else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) 
		mouseButton = GLUT_MIDDLE_BUTTON;
	else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
		mouseButton = GLUT_RIGHT_BUTTON;

	mousePositionX0 = x;
	mousePositionY0 = y;
	return ;
}

//Rendering
void Render_BoxAndAxes() {

	float axiswidth = axislen / 100;

	glMatrixMode (GL_MODELVIEW);
	
	glColor3f (1,1,1);
	
	glPushMatrix();
	//bounding box
	glTranslatef(objCenter[0], objCenter[1], objCenter[2]);
	glutWireCube(axislen);
	glTranslatef(-axislen/2, -axislen/2, -axislen/2);
	glutSolidSphere(axiswidth*1.5, 10, 10);

	//x-axis
	glColor3f (1,0,0);
	glPushMatrix();
	glRotatef(90, 0, 1, 0);
	gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
	glPopMatrix();

	glPushMatrix();	
	glTranslatef(axislen,0,0);	
	glRotatef(90, 0, 1, 0);
	glutWireCone(axiswidth*1.5, axiswidth*3, 10, 10);
	glPopMatrix();
	

	//y-axis
	glColor3f (0,1,0);	
	glPushMatrix();
	glTranslatef(0,axislen,0);
	glRotatef(-90, 1, 0, 0);
	glutWireCone(axiswidth*1.5,axiswidth*3,10,10);
	glPopMatrix();

	glPushMatrix();
	glRotatef(-90, 1, 0, 0);
	gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
	glPopMatrix();

	//z-axis
	glColor3f (0,0,1);
	glPushMatrix();
	glTranslatef(0,0,axislen);	
	glRotatef(-90, 0, 0, 1);
	glutWireCone(axiswidth*1.5,axiswidth*3,10,10);
	glPopMatrix();

	glPushMatrix();	
	gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
	glPopMatrix();

	glPopMatrix();
}

void Render_Ball(double xyz[3], double radius, bool red_blue){ //fasle: red; true: blue
	glPushMatrix();
	glTranslated(xyz[0], xyz[1], xyz[2]);
	if (!red_blue)
		glColor3f(0.5, 0, 0);
	else
		glColor3f(0, 0, 0.5);
	glutSolidSphere(radius, 4, 4);
	glPopMatrix();
}

void Render_SharpEdge(){
	float currentColor[4];
	glGetFloatv(GL_CURRENT_COLOR, currentColor);
	glColor3f(0, 0, 0.8);	
	glLineWidth(5);
	glBegin(GL_LINES);
	for (int i = 0; i < sharpEdges.size(); ++i){
		if (sharpEdges[i]){
			Edge * e = myMesh->indEdge(i);
			Halfedge * he = e->he(0);
			Point & p1 = he->source()->point();
			Point & p2 = he->target()->point();
			glVertex3dv(p1.v);
			glVertex3dv(p2.v);
		}
	}
	glEnd();	
	glLineWidth(1);
	glColor3f(currentColor[0], currentColor[1], currentColor[2]);
}

void Render_MinMaxK(){
	float currentColor[4];
	glGetFloatv(GL_CURRENT_COLOR, currentColor);
	double ballSize = axislen / 100.0;
	for (int i = 0; i < maxmin_kg.size(); ++i){
		if (maxmin_kg[i] == 0)
			continue;
		else if (maxmin_kg[i] == 1) {//red ball
			Vertex * v = myMesh->indVertex(i);
			Render_Ball(v->point().v, ballSize, false);			
		}
		else if (maxmin_kg[i] == -1){//blue ball
			Vertex * v = myMesh->indVertex(i);
			Render_Ball(v->point().v, ballSize, true);
		}
	}
	glColor3f(currentColor[0], currentColor[1], currentColor[2]);
}

void Render_Mesh(){

	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glBindTexture(GL_TEXTURE_2D, CBTex->texName);

	glColor3f(0.7f, 0.7f, 0.7f);	
	//traverse all the face and draw them 
	glBegin(GL_TRIANGLES);
	for (MeshFaceIterator fit(myMesh); !fit.end(); ++fit){
		Face * f = *fit;
		if (useNormal == 1) {
			Point & pFNorm = fNList[f->index()];
			glNormal3f(pFNorm[0], pFNorm[1], pFNorm[2]);
		}
		for (FaceVertexIterator vit(f); !vit.end(); ++vit){
			Vertex * v = *vit;
			if (useNormal == 2){
				Point & pVNorm = vNList[v->index()];
				glNormal3f(pVNorm[0], pVNorm[1], pVNorm[2]);
			}			
			Point & pt = v->point();
			// add the texture to the mesh
			double isoU = tex_U[v->index()];
			double isoV = tex_V[v->index()];
			if (showISO_u)
			{
				glTexCoord2f(isoU,0); glVertex3f(pt[0], pt[1], pt[2]);
			}
			else if (showISO_v)
			{
				glTexCoord2f(0, isoV); glVertex3f(pt[0], pt[1], pt[2]);
			}
			else
			{
				glTexCoord2f(isoU, isoV); glVertex3f(pt[0], pt[1], pt[2]);
			}
			
		}
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);
}

void ComputeBoundingBox() {
	objCenter[0]=objCenter[1]=objCenter[2]=0;
	boxMin[0]=boxMin[1]=boxMin[2]=1e5;
	boxMax[0]=boxMax[1]=boxMax[2]=-1e5;
	int vNum = myMesh->numVertices();
	for (MeshVertexIterator vit(myMesh); !vit.end(); ++vit){
		Vertex * v = *vit;
		for (int j = 0; j < 3; ++j) {
			double value = v->point().v[j];
			objCenter[j] += value;
 			if (boxMax[j] < value)
				boxMax[j] = value;
			if (boxMin[j] > value)
				boxMin[j] = value;
		}
	}	
	axislen=sqrt((boxMax[2]-boxMin[2])*(boxMax[2]-boxMin[2])+(boxMax[1]-boxMin[1])*(boxMax[1]-boxMin[1])+(boxMax[0]-boxMin[0])*(boxMax[0]-boxMin[0]));

	objCenter[0] /= vNum;
	objCenter[1] /= vNum;
	objCenter[2] /= vNum;

	/*
	////You can also do a normalization here
	//for (int i = 0; i < vNum; ++i){
	//	for (int j = 0; j < 3; ++j)
	//		vert_xyz[j][i] -= objCenter[j];
	//}	
	//for (int i = 0; i < 3; ++i){
	//	boxMin[i] -= objCenter[i];
	//	boxMax[i] -= objCenter[i];
	//}
	//objCenter[0] = objCenter[1] = objCenter[2] = 0;
	*/

	cameraPosition[0]=objCenter[0];
	cameraPosition[1]=objCenter[1];
	cameraPosition[2]=objCenter[2] + axislen*1.5;

	std::cout << objCenter[0] << " " << objCenter[1] << " " <<objCenter[2] << " "
		<< cameraPosition[0] << " " << cameraPosition[1] << " " << cameraPosition[2] << "\n";
}

void MyInit()
{	
	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);      
	glEnable(GL_DEPTH_TEST);
	glClearColor(0,0,0,0);
	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT , GL_FILL);
	
	CBTex->GenTex();
	obj = gluNewQuadric();	//only for drawing spheres and cones		

	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHT0);

	// Create light components
	GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	GLfloat	diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
	GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat position[] = {cameraPosition[0], cameraPosition[1], cameraPosition[2], 1.0f }; // the light is on the camera position

	// Assign created components to GL_LIGHT0
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
	glLightfv(GL_LIGHT0, GL_POSITION, position);

	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
}

int main(int argc, char** argv)
{

	if (argc != 2) {
		std::cout << "HW3: Read a .obj file and render it in its boundingbox.\n";
		std::cout << "Usage: " << argv[0] << " input_mesh.obj\n";
		return -1;
	}

	myMesh = new Mesh;
	bool flag = myMesh->readOBJFile(argv[1]);
	if (!flag) {
		std::cerr << "Fail to read " << argv[1] << "\n";
		return -2;
	}

	showSharpEdge = false;
	showK_g = false;
	showISO_u = false;
	showISO_v = false;

	useNormal = 2;

	ComputeBoundingBox();
	ComputeFaceNormals();
	ComputeCornerAngles();
	ComputeVertNormals();
	ComputeGaussianCurvature();
	getMaxMinKG();
	ComputeSharpEdges();

	if (myMesh->numVertices() == 81)
	{// This is task 1 for Cyliner9.obj
		ComputeHarmonicScalarField();
	}
	else
	{// this is task 2 for Susan.obj
		getBoundary();
		BoundaryPara();
		HarmonicSolver();
	}

	//OpenGL Routines
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(700, 700);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Homework 3 submitted by Yuqi Ding");
	MyInit();

	//Register Callback Functions
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);
	glutKeyboardFunc(keyPressHandler);

	glutMainLoop();

	delete myMesh;
	return 0;

}