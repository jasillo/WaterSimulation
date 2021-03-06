#include <Windows.h>
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <glm\glm.hpp>

#define M_PI 3.14159265358979323846 

using namespace std;
using namespace glm;

// "Particle-Based Fluid Simulation for Interactive Applications"
// solver parameters
const static vec2 G(0.f, 12000 * -9.8f); // external (gravitational) forces
const static float REST_DENS = 1000.f; // rest density
const static float GAS_CONST = 2000.f; // const for equation of state
const static float H = 16.f; // kernel radius
const static float HSQ = H*H; // radius^2 for optimization
const static float MASS = 65.f; // assume all particles have the same mass
const static float VISC = 250.f; // viscosity constant
const static float DT = 0.001f; // integration timestep

								 // smoothing kernels defined in M�ller and their gradients
const static float POLY6 = 315.f / (65.f*M_PI*pow(H, 9.f));
const static float SPIKY_GRAD = -45.f / (M_PI*pow(H, 6.f));
const static float VISC_LAP = 45.f / (M_PI*pow(H, 6.f));

// simulation parameters
const static float EPS = H; // boundary epsilon
const static float BOUND_DAMPING = -0.5f;

float currentTime;
float oldTime;
int numOfFrames;

// particle data structure
// stores position, velocity, and force for integration
// stores density (rho) and pressure values for SPH
struct Particle {
	Particle(float _x, float _y) : x(_x, _y), v(0.f, 0.f), f(0.f, 0.f), rho(0), p(0.f) {}
	vec2 x, v, f;
	float rho, p;
};

// solver data
vector<Particle> particles;
vector<int> grilla[50][50];
vector<vector<int>> distanciaH;
vector<vector<int>> distancia2H;

void limpiargrilla() {
	for (size_t i = 0; i < 50; i++)
	{
		for (size_t j = 0; j < 50; j++)
		{
			grilla[i][j].clear();
		}
	}
}

void llenargrilla() {
	int x, y;
	for (size_t i = 0; i < particles.size(); i++)
	{
		x = particles[i].x.x / 16;
		y = particles[i].x.y / 16;
		grilla[x][y].push_back(i);
	}
}

// interaction
const static int MAX_PARTICLES = 2500;
const static int DAM_PARTICLES = 100;
const static int BLOCK_PARTICLES = 250;

// rendering projection parameters
const static int WINDOW_WIDTH = 600;
const static int WINDOW_HEIGHT = 600;
const static double VIEW_WIDTH = 1.5*800.f;
const static double VIEW_HEIGHT = 1.5*600.f;

void InitSPH(void)
{	
	
	float x = 300.0, y = 400.0;
	float temp = 1;
	particles.push_back(Particle(x, y));

	for (float i = 1.0; i < 4.0; i += 1.0)
	{
		temp *= 2;
		for (float j = 0.0; j < 2 * M_PI; j += M_PI/temp)
		{
			//cout << j << endl;
			particles.push_back(Particle(x + i*H*cos(j), y + i*H*sin(j)));
			//cout << x + i*H*cos(j) << " - " << y + i*H*sin(j) << endl;
		}
	}

	for (float y = EPS; y <= 5 * EPS; y += H)	
	{
		for (float x = EPS; x < 800 - EPS; x += H)
		{
			if (particles.size() < 250)
				particles.push_back(Particle(x, y));
			else
				break;
		}
	}
	//cout << particles.size() << endl;
	oldTime = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
	llenargrilla();
}

void Integrate(void)
{
	for (auto &p : particles)
	{
		// forward Euler integration
		p.v += DT*p.f / p.rho;
		p.x += DT*p.v;

		// enforce boundary conditions
		if (p.x.x - EPS < 0.0f)
		{
			p.v.x *= BOUND_DAMPING;
			p.x.x = EPS;
		}
		if (p.x.x + EPS > 800)
		{
			p.v.x *= BOUND_DAMPING;
			p.x.x = 800 - EPS;
		}
		if (p.x.y - EPS < 0.0f)
		{
			p.v.y *= BOUND_DAMPING;
			p.x.y = EPS;
		}
		if (p.x.y + EPS > 600)
		{
			p.v.y *= BOUND_DAMPING;
			p.x.y = 600 - EPS;
		}
	}
	limpiargrilla();
	llenargrilla();
}

void ComputeDensityPressure(void)
{
	/*
	for (auto &pi : particles)
	{
		pi.rho = 0.f;
		for (auto &pj : particles)
		{
			vec2 rij = pj.x - pi.x;
			float r2 = pow(distance(pj.x, pi.x),2.0f);

			if (r2 < HSQ)
			{
				// this computation is symmetric
				pi.rho += MASS*POLY6*pow(HSQ - r2, 3.f);
			}
		}
		pi.p = GAS_CONST*(pi.rho - REST_DENS);
	}**/
	
	int x, y;
	for (size_t i = 0; i < particles.size(); i++)
	{
		x = particles[i].x.x / 16;
		y = particles[i].x.y / 16;
		particles[i].rho = 0.f;

		for (size_t m = x -1; m <= x + 1; m++)
		{
			for (size_t n = y-1; n <= y + 1; n++)
			{
				if (m >= 0 && n >= 0 && m < 50 && n < 50) {
					for (size_t j = 0; j < grilla[m][n].size(); j++)
					{
						float r2 = pow(distance(particles[(grilla[m][n])[j]].x, particles[i].x), 2.0f);						
						if (r2 < HSQ)
						{
							particles[i].rho += MASS*POLY6*pow(HSQ - r2, 3.f);
						}
					}
				}					
			}
		}
		particles[i].p = GAS_CONST*(particles[i].rho - REST_DENS);
	}
}

void ComputeForces(void)
{
	/*
	for (auto &pi : particles)
	{
		vec2 fpress(0.f, 0.f);
		vec2 fvisc(0.f, 0.f);
		for (auto &pj : particles)
		{
			if (&pi == &pj)
				continue;

			vec2 rij = pj.x - pi.x;
			float r = distance(pj.x, pi.x);

			if (r < H)
			{
				// compute pressure force contribution
				fpress += -normalize(rij)*(pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD*pow(H - r, 2.f);
				// compute viscosity force contribution
				fvisc += VISC*MASS*(pj.v - pi.v) / pj.rho * VISC_LAP*(H - r);
			}
		}
		fpress *= MASS;
		vec2 fgrav = G * pi.rho;
		pi.f = fpress + fvisc + fgrav;
	}*/
	int x, y;
	for (size_t i = 0; i < particles.size(); i++)
	{
		vec2 fpress(0.f, 0.f);
		vec2 fvisc(0.f, 0.f);
		x = particles[i].x.x / 16;
		y = particles[i].x.y / 16;
		
		for (size_t m = x - 1; m <= x + 1; m++)
		{
			for (size_t n = y - 1; n <= y + 1; n++)
			{
				if (m >= 0 && n >= 0 && m < 50 && n < 50) {
					for (size_t j = 0; j < grilla[m][n].size(); j++)
					{
						if (i == grilla[m][n][j]) {
							continue;
						}
						vec2 rij = particles[(grilla[m][n])[j]].x - particles[i].x;
						float r = distance(particles[(grilla[m][n])[j]].x, particles[i].x);

						if (r < H)
						{
							fpress += -normalize(rij)*(particles[i].p + particles[(grilla[m][n])[j]].p) / (2.f * particles[(grilla[m][n])[j]].rho) * SPIKY_GRAD*pow(H - r, 2.f);							
							fvisc += VISC*MASS*(particles[(grilla[m][n])[j]].v - particles[i].v) / particles[(grilla[m][n])[j]].rho * VISC_LAP*(H - r);
						}
					}
				}
			}
		}
		fpress *= MASS;
		vec2 fgrav = G * particles[i].rho;
		particles[i].f = fpress + fvisc + fgrav;
	}
}

void Update(void)
{
	glutPostRedisplay();
}

void InitGL(void)
{
	glClearColor(1.0f, 1.0f, 1.0f, 1);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(4);
	glMatrixMode(GL_PROJECTION);
	//oldTime = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
	numOfFrames = 0;
}

void Render(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	currentTime = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
	glLoadIdentity();
	glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);
	
	//for (size_t i = 0; i < 3; i++)
	//{
	ComputeDensityPressure();
	ComputeForces();
	Integrate();
	//}
	
	glTranslated(100,100,0);
	glColor4f(0.f, 0.3f, 0.7f, 1);
	/*glBegin(GL_POINTS);
	for (auto &p : particles)
		glVertex2f(p.x.x, p.x.y);
	glEnd();*/
	//for (auto &p : particles)
		//cout << p.x.x << endl;
	//cout<< glutGet(GLUT_ELAPSED_TIME) / 1000.0 <<endl;
	if (currentTime - oldTime <= 10.0) //10seg	
		numOfFrames++;	
	else {
		cout << numOfFrames << endl;
	}
		

	//oldTime = currentTime;
	glutSwapBuffers();
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

int main(int argc, char** argv)
{
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInit(&argc, argv);
	glutCreateWindow("SPH");
	glutDisplayFunc(Render);
	glutIdleFunc(Update);
	glutKeyboardFunc(window_key);

	InitGL();
	InitSPH();

	glutMainLoop();
	return 0;
}