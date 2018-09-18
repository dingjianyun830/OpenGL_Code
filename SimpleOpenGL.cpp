// SimpleOpenGL.cpp : Defines the entry point for the console application.
// Author: Yuqi Ding
// Professor Jinwei Ye
// CSC4356
// Programing Assignment1
//
#include "GL\glut.h"
#include <math.h>
#include <time.h>   

int win_H, win_W;
time_t timer;
struct tm curr_time;

void reshape(int w, int h)
{
	glViewport(0, 0, w, h);       /* Establish viewing area to cover entire window. */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, -1, 1);
	glScalef(1, -1, 1);
	glTranslatef(0, -h, 0);
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	// Insert your own code here
	float Clock_Center[2] = { 0.0f, 0.0f };
	float radius = 150.0f;
	float Baseline[2] = { 0.0f + radius , 0.0f };
	float Hourhand1[2] = { 266.0f - 256.0f, 256.0f - 256.0f };
	float Hourhand2[2] = { 246.0f - 256.0f, 256.0f - 256.0f };
	float Hourhand3[2] = { 261.0f - 256.0f, 176.0f - 256.0f };
	float Hourhand4[2] = { 251.0f - 256.0f, 176.0f - 256.0f };

	float Minutehand1[2] = { 262.0f - 256.0f,256.0f - 256.0f };
	float Minutehand2[2] = { 250.0f - 256.0f,256.0f - 256.0f };
	float Minutehand3[2] = { 256.0f - 256.0f,196.0f - 256.0f };

	float Secondhand1[2] = { 256.0f - 256.0f, 256.0f - 256.0f };
	float Secondhand2[2] = { 256.0f - 256.0f, 146.0f - 256.0f };

	glPushMatrix();
	// translate all object to the new original position
	glTranslatef(256.0f, 256.0f, 0.0f);

	// Draw the Clock face
	float Baseline_R1[2] = { 0.0f, 0.0f };
	float Baseline_R2[2] = { 0.0f, 0.0f };
	// set the color of clock face
	glColor3f(1.0f, 1.0f, 1.0f);

	glBegin(GL_TRIANGLE_FAN);
	// from the x axis to rotate
	// a*pi/180 = r
	float PI = 3.1415926;
	for (int i = 0; i < 360; i++)
	{
		float theta1 = i * PI / 180.0;
		float theta2 = (i + 1) * PI / 180.0;
		Baseline_R1[0] = Baseline[0] * cos(theta1) - Baseline[1] * sin(theta1);
		Baseline_R1[1] = Baseline[0] * sin(theta1) + Baseline[1] * cos(theta1);

		Baseline_R2[0] = Baseline[0] * cos(theta2) - Baseline[1] * sin(theta2);
		Baseline_R2[1] = Baseline[0] * sin(theta2) + Baseline[1] * cos(theta2);

		glVertex2fv(Clock_Center);
		glVertex2fv(Baseline_R1);
		glVertex2fv(Baseline_R2);
	}
	glEnd();

	// Draw the Hourhand
	float hour_rotate = (float)curr_time.tm_hour*30.0f + (float)curr_time.tm_min / 2.0f + (float)curr_time.tm_sec / 120.0f;
	glPushMatrix();
	glRotatef(hour_rotate, 0, 0, 1);
	glColor3f(0.0f, 1.0f, 0.0f);
	glBegin(GL_TRIANGLE_STRIP);
	glVertex2fv(Hourhand1);
	glVertex2fv(Hourhand2);
	glVertex2fv(Hourhand3);
	glVertex2fv(Hourhand4);
	glEnd();
	glPopMatrix();

	// Draw the Minutehand
	float minute_rotate = (float)curr_time.tm_min*6.0f + (float)curr_time.tm_sec / 10.0f;
	glPushMatrix();
	glRotatef(minute_rotate, 0, 0, 1);
	glColor3f(0.0f, 0.0f, 1.0f);
	glBegin(GL_TRIANGLES);
	glVertex2fv(Minutehand1);
	glVertex2fv(Minutehand2);
	glVertex2fv(Minutehand3);
	glEnd();
	glPopMatrix();

	// Draw the Secondhand
	float second_rotate = (float)curr_time.tm_sec*6.0f;
	glPushMatrix();
	glRotatef(second_rotate, 0, 0, 1);
	glColor3f(0.1f, 0.2f, 0.3f);
	glBegin(GL_LINES);
	glVertex2fv(Secondhand1);
	glVertex2fv(Secondhand2);
	glEnd();
	glPopMatrix();

	glPopMatrix();
	// Your code ends here

	glutSwapBuffers(); // swap the back buffer to front


}

void TimeEvent(int time_val)
{
	time(&timer); // get the current date and time from system
	localtime_s(&curr_time, &timer);
	glutPostRedisplay();
	glutTimerFunc(30, TimeEvent, 1);// By using a timed event, your application should run at the same speed on any machine.
}

int main(int argc, char **argv)
{
	GLenum type;

	glutInit(&argc, argv);

	type = GLUT_DEPTH;
	type |= GLUT_RGB;
	type |= GLUT_DOUBLE;
	glutInitDisplayMode(type);

	time(&timer); // get current date and time
	localtime_s(&curr_time, &timer);

	// set window size and create a window for rendering
	win_W = 512;
	win_H = 512;
	glutInitWindowSize(win_H, win_W);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("My clock");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);

	glutTimerFunc(30, TimeEvent, 1);
	glutMainLoop();
	return 0;
}
