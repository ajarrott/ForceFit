#ifndef _ATOM_H
#define _ATOM_H

#include <string>
#include <iostream>
#include "pugixml/pugixml.hpp"

using std::string;

class Atom {
public:
	//These should be straight from the periodic table
	int number;
	std::string symbol;
	float mass;

	//This will be provided by the user
	float charge;

	//Position
	float x, y, z;

	//Force for this atom
	float forcex, forcey, forcez;

	Atom();

	void write(std::ostream &, int indent);
	void load(pugi::xml_node node);
};

#endif
