#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <vector>
#include <iostream>

#include "pugixml/pugixml.hpp"
#include "atom.h"

class Geometry {
public:
	std::vector<Atom> atoms;

	int mode;
	float step;
	float energy;

	Geometry();

	void write(std::ostream &, int indent);
	void load(pugi::xml_node node);
};

#endif
