#define GLUT_DISABLE_ATEXIT_HACK
#include <windows.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <GL\glut.h>
#include <glm\glm.hpp>
#include <glm\gtc\constants.hpp>
#include <vector>

using namespace std;

#define RED 0
#define GREEN 0
#define BLUE 0
#define ALPHA 1
#define Xmax 10.0f
#define Ymax 25.0f
#define Zmax 10.0f

#define SPACE .5f //espacio entre particulas

const static float H = 4.0f; // area de smoothness
const static float Mass = 6500.0f; 
const static float Visc = 250.0f; //indice de viscocidad
const static float REST_DENS = 200.0f; //resistencia (densidad) natural del agua
const static float GAS_CONST = 100.0f; //
const static float VISC_LAP = 45.0f / (glm::pi<float>() * glm::pow(H, 6.0));
const static float BOUND = -0.5; //fuerza con que se refleja en bordes
const static float H2 = H*H;
const static glm::vec3 Gravity(0.f,-400.0f,0.f);

static float kernelConstant = 315.0f / (64.0f * glm::pi<float>() * glm::pow(H, 9.0f)) ;
//static float gradientConstant = -945 / (32 * glm::pi<float>() * glm::pow(H, 9.0f));
static float gradientConstant = -45.0f / (glm::pi<float>() * glm::pow(H, 6.0f));
static float LaplaceConstant = -945.0f / (32.0f * glm::pi<float>() * glm::pow(H, 9.0f));
static float h2 = glm::pow(H, 2);
float currentTime = 0.0;
float newTime = 0.0;

GLvoid initGL();
GLvoid window_display();
GLvoid window_reshape(GLsizei width, GLsizei height);
GLvoid window_key(unsigned char key, int x, int y);
//function called on each frame
GLvoid window_idle();

int times = 5;

struct Particle{
	glm::vec3 x; //posision
	glm::vec3 u; //velocidad
	glm::vec3 f; //fuerza
	float rho; //densidad
	float p; //presion

	Particle(float posX, float posY, float posZ): 
		x(glm::vec3(posX,posY,posZ)),
		u(glm::vec3(0.0f, 0.0f, 0.0f)),
		f(glm::vec3(0.0f,0.0f,0.0f))
	{
		rho = 0;
		p = 0;
	}
};

vector<Particle> particles;

void createParticles(int size)
{
	float limit = 10.0;
	for (float y = 21.0; y < Ymax; y += SPACE) {
		for (float x = 1.0 ; x <= 1.0+ 4.0*SPACE; x += SPACE) {
			for (float z = 1.0; z <= 1.0+ 4.0*SPACE; z+= SPACE) {
				if (particles.size() >= size)
					return;
				particles.push_back(struct Particle(x,y,z));
			}			
		}
	}	
}

float kernel(float r)
{
	//cout << r << " ";
	if (r > H)
		return 0.0f;
	float res =  kernelConstant * glm::pow(H2 - glm::pow(r,2.0f), 3.0f);
	//cout << res << " ";
	if (res < 0.0f)
		return 0.0f;
	if (res > 1.0f)
		return 1.0f;
	return res;
}

glm::vec3 gradient(glm::vec3 r, float d)
{
	if (d == 0)
		return glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 res =  gradientConstant * r / d * glm::pow(H - d, 2.0f);
	//cout << res.x << " " << res.y << " " << res.z << endl;
	return res;
}

float laplace(float r)
{
	return LaplaceConstant * (H2 - glm::pow(r, 2.0f)) * (3 * H2 * - 7 * glm::pow(r,2.0f));
}

void calculateDensityPresure() 
{
	for (auto &pi : particles)
	{
		pi.rho = 0.0f;
		for (auto &pj : particles)
		{
			float r;
			if (&pi.x == &pj.x)
				r = 0.0f;
			else
				r = glm::distance(pj.x, pi.x);
			// cout << r << " ";
			if (r < H)
			{
				pi.rho += kernel(r);				
			}
		}
		//cout << pi.rho << " ";
		pi.rho *= Mass;
		
		pi.p = GAS_CONST * (pi.rho - REST_DENS);
	}
}

void calculateForces() 
{
	for (auto &pi : particles)
	{
		glm::vec3 fpress (0.0f, 0.0f, 0.0f);
		glm::vec3 fvisc (0.0f, 0.0f, 0.0f);
		for (auto &pj : particles)
		{
			if (&pi.x == &pj.x)
				continue;
			float r = glm::distance(pj.x, pi.x);
			if (r < H)
			{
				fpress += (pi.p + pj.p) / (2.0f * pj.rho) * gradient(pi.x - pj.x, r);
				fvisc += (pj.u - pi.u) / pj.rho * laplace(r);
			}
		}
		fvisc *= Visc * Mass;
		fpress *= -Mass;
		//cout << fpress.x << " " << fpress.y << endl;
		glm::vec3 fgrav = Gravity * pi.rho;
		//cout << fgrav.x << " " << fgrav.y << endl;
		pi.f = fpress + fgrav;
	} 
}

void integrate()
{
	//times--;
	//if (times < 0)
		//return;
	newTime = glutGet(GLUT_ELAPSED_TIME) /10000.0f; // time elapsed
	float DT =  newTime - currentTime;
	
	for (auto &pi : particles)
	{
		//cout << pi.f.x << " " << pi.f.y << " " << pi.f.z << endl;
		//cout << pi.rho << " ";
		pi.u += DT * pi.f / pi.rho;
		//cout << pi.u.x << " " << pi.u.y << " " << pi.u.z << endl;
		pi.x += DT * pi.u;
		
		if (pi.x.x < 0.0)
		{
			pi.u.x *= -0.5;
			pi.x.x = 0.0;
		}
		if (pi.x.x > Xmax)
		{
			pi.u.x *= -0.5;
			pi.x.x = Xmax;
		}
		if (pi.x.y < 0.0)
		{
			pi.u.y *= -0.5;
			pi.x.y = 0.0;
		}
		if (pi.x.y > Ymax)
		{
			pi.u.y *= -0.5;
			pi.x.y = Ymax;
		}
		if (pi.x.z < 0.0)
		{
			pi.u.z *= -0.5;
			pi.x.z = 0.0;
		}
		if (pi.x.z > Zmax)
		{
			pi.u.z *= -0.5;
			pi.x.z = Zmax;
		}
	}
	currentTime = newTime;
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
	glPointSize(8.0f);
	createParticles(70);
	/*for (size_t i = 0; i < particles.size(); i++)
	{
		cout << particles[i].x.x << " " << particles[i].x.y << " " << particles[i].x.z << endl;
	}*/
	cout << particles.size()<<endl;
}




GLvoid window_display()
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	 
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(-25.0f, 25.0f, -25.0f, 25.0f, -25.0f, 25.0f);
	gluPerspective(60.0, 1.0, 20, 40);
	glTranslated(-Xmax/2.0, -Ymax/2.0, -32.0);
	
	calculateDensityPresure();
	calculateForces();
	integrate(); 
	
	//glutSolidSphere(1, 10, 10);
	glBegin(GL_LINE_STRIP);
	glVertex3f(0,0, 0);	
	glVertex3f(Xmax, 0, 0);
	glVertex3f(Xmax, 0, Zmax);
	glVertex3f(0,0 ,Zmax );
	glVertex3f(0, 0, 0);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, Ymax, 0);
	glVertex3f(Xmax, 0, 0);
	glVertex3f(Xmax, Ymax, 0);
	glVertex3f(Xmax, 0, Zmax);
	glVertex3f(Xmax, Ymax, Zmax);
	glVertex3f(0, 0, Zmax);
	glVertex3f(0, Ymax, Zmax);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex3f(0, Ymax, 0);
	glVertex3f(Xmax, Ymax, 0);
	glVertex3f(Xmax, Ymax, Zmax);
	glVertex3f(0, Ymax, Zmax);
	glVertex3f(0, Ymax, 0);
	glEnd();

	for (auto &Pi : particles) {
		glPushMatrix();
		glTranslated(Pi.x.x, Pi.x.y, Pi.x.z);
		glutSolidSphere(0.25, 7, 7);
		glPopMatrix();
		//cout << Pi.x.x <<" "<< Pi.x.y <<" "<< Pi.x.z << endl;
		//glBegin(GL_POINTS);
		//glVertex3f(Pi.x.x, Pi.x.y, Pi.x.z);
		//glEnd();
	}

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
