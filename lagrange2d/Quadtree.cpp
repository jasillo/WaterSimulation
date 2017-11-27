#include "Quadtree.h"




void Quadtree::createNodes(Node * n, float size)
{
	float size = n->fin.x - n->ini.x;
	if (size <= h)
		return;
	for (size_t i = 0; i < 4; i++)
	{
		n->sons.push_back(struct Node());
	}

	n->sons[0].ini = n->ini;
	n->sons[0].fin = n->ini + glm::vec2(size, size);

	n->sons[1].ini = n->fin - glm::vec2(size, size);
	n->sons[1].fin = n->fin;

	n->sons[2].ini = n->ini + glm::vec2(0.0f,size);
	n->sons[2].fin = n->fin - glm::vec2(size,0.0f);

	n->sons[3].ini = n->ini + glm::vec2(size,0.0f);
	n->sons[3].fin = n->fin - glm::vec2(0.0f, size);

	for (size_t i = 0; i < 4; i++)
	{
		createNodes(&(n->sons[i]), size / 2.0f);
	}
}

Quadtree::~Quadtree()
{
}
