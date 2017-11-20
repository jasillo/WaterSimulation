#define GLUT_DISABLE_ATEXIT_HACK
#include <windows.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "particles.h"
#include <algorithm>

#include <GL\glut.h>
#include <vector>

using namespace std;

#define RED 0
#define GREEN 0
#define BLUE 0
#define ALPHA 1

GLvoid initGL();
GLvoid window_display();
GLvoid window_reshape(GLsizei width, GLsizei height);
GLvoid window_key(unsigned char key, int x, int y);
//function called on each frame
GLvoid window_idle();

Grid grid(9.8, 50, 50, 1);
Particles particles(grid);


using namespace std;

float fluidphi(Grid &grid, float x, float y)
{
	//return y-0.5*grid.ly;
	//return min(sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.625*grid.ly))-0.02*grid.ly, y-0.6*grid.ly);
	//return min(sqrt(sqr(x-0.3333*grid.lx)+sqr(y-0.71*grid.ly))-0.3*grid.ly, y-0.2*grid.ly);
	//return max(y-0.8*grid.ly, -sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.2*grid.ly))+0.1*grid.lx);
	//return sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.75*grid.ly))-0.15*grid.lx;
	
	return min(y - 0.05f*grid.ly, sqrt(pow(x - 0.5f*grid.lx, 2.0f) + pow(y - 0.5f*grid.ly, 2.0f)) - 0.05f*grid.lx);
	//return 0.75*grid.lx-x;
	//return max(x-0.75*grid.lx, 0.25*grid.lx-x, y-0.75*grid.ly, 0.25*grid.ly-y);
}

void project(Grid &grid, float &x, float &y, float current, float target)
{
	float dpdx = (fluidphi(grid, x + 1e-4, y) - fluidphi(grid, x - 1e-4, y)) / 2e-4;
	float dpdy = (fluidphi(grid, x, y + 1e-4) - fluidphi(grid, x, y - 1e-4)) / 2e-4;
	float scale = (target - current) / sqrt(dpdx*dpdx + dpdy*dpdy);
	x += scale*dpdx;
	y += scale*dpdy;
}

void init_water_drop(Grid &grid, Particles &particles, int na, int nb)
{
	int i, j, a, b;
	float x, y, phi;

	for (i = 1; i<grid.marker.nx - 1; ++i) {
		for (j = 1; j<grid.marker.ny - 1; ++j) {
			for (a = 0; a<na; ++a) {
				for (b = 0; b<nb; ++b) {
					x = (i + (a + 0.1 + 0.8*rand() / (double)RAND_MAX) / na)*grid.h;
					y = (j + (b + 0.1 + 0.8*rand() / (double)RAND_MAX) / nb)*grid.h;
					phi = fluidphi(grid, x, y);
					if (phi>-0.25*grid.h / na)
						continue;
					else if (phi>-1.5*grid.h / na) {
						project(grid, x, y, phi, -0.75*grid.h / na);
						phi = fluidphi(grid, x, y);
						project(grid, x, y, phi, -0.75*grid.h / na);
						phi = fluidphi(grid, x, y);
					}
					particles.add_particle(Vec2f(x, y), Vec2f(0, 0));
				}
			}
		}
	}
}

void advance_one_step(Grid &grid, Particles &particles, double dt)
{
	for (int i = 0; i<5; ++i)
		particles.move_particles_in_grid(0.2*dt);
	particles.transfer_to_grid();
	grid.save_velocities();
	grid.add_gravity(dt);
	grid.compute_distance_to_fluid();
	grid.extend_velocity();
	grid.apply_boundary_conditions();
	grid.make_incompressible();
	grid.extend_velocity();
	grid.get_velocity_update();
	particles.update_from_grid();
}

void advance_one_frame(Grid &grid, Particles &particles, double frametime)
{
	double t = 0;
	double dt;
	bool finished = false;
	while (!finished) {
		dt = 2 * grid.CFL();
		if (t + dt >= frametime) {
			dt = frametime - t;
			finished = true;
		}
		else if (t + 1.5*dt >= frametime)
			dt = 0.5*(frametime - t);
		advance_one_step(grid, particles, dt);
		t += dt;
	}
}




int main(int argc, char **argv)
{
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);


	glutInitWindowSize(600, 600);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("TP 2 : Transformaciones");


	initGL();
	glutDisplayFunc(&window_display);
	glutReshapeFunc(&window_reshape);
	glutKeyboardFunc(&window_key);
	//function called on each frame
	glutIdleFunc(&window_idle);

	glutMainLoop();

	return 1;
}



GLvoid initGL()
{
	GLfloat position[] = { 0.0f, 5.0f, 10.0f, 0.0 };

	//enable light : try without it
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glEnable(GL_LIGHTING);
	//light 0 "on": try without it
	glEnable(GL_LIGHT0);

	//shading model : try GL_FLAT
	glShadeModel(GL_SMOOTH);

	glEnable(GL_DEPTH_TEST);

	//enable material : try without it
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glClearColor(RED, GREEN, BLUE, ALPHA);
	glPointSize(3.0f);	

	init_water_drop(grid, particles, 2, 2);
}




GLvoid window_display()
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(-25.0f, 25.0f, -25.0f, 25.0f, -25.0f, 25.0f);
	gluPerspective(60.0, 1.0, 20, 40);
	glTranslated(-10.0, 0.0, -32.0);

	advance_one_frame(grid, particles, 1. / 30);

	glBegin(GL_POINTS);
	for (size_t i = 0; i < particles.x.size(); i++)
	{
		glVertex2f(particles.x[i].v[0]*20, particles.x[i].v[1]*20);
	}
	glEnd();



	glutSwapBuffers();

	glFlush();
}

GLvoid window_reshape(GLsizei width, GLsizei height)
{
	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, 1.0, 10, 50);

	glMatrixMode(GL_MODELVIEW);
}

GLvoid window_key(unsigned char key, int x, int y)
{
	switch (key) {
	case 'r':
		break;
	default:
		printf("La touche %d non active.\n", key);
		break;
	}
}


//function called on each frame
GLvoid window_idle()
{
	glutPostRedisplay();
}
