#pragma once
#include <glm\glm.hpp>
#include <vector>

using namespace std;

struct Node
{
	vector<glm::vec2*> points;
	vector<Node> sons;
	glm::vec2 ini;
	glm::vec2 fin;
};

class Quadtree
{
private:
	glm::vec2 maximo, minimo;
	Node root;
	float h;

public:
	Quadtree(float size_, float h_)  
	{
		h = h_;
		root.ini = glm::vec2(0.0, 0.0);
		root.fin = glm::vec2(size_, size_);
		createNodes(&root, size_/2.0f);
	};

	void createNodes(Node *n, float h);
	~Quadtree();
};

