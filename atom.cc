#include <iostream>
#include <cstdlib>

#include "atom.h"

Atom::Atom() {
	charge = 0;
}

void Atom::write(std::ostream & out, int indent){
	for(int i = 0; i < indent; i++){
		out << "\t";
	}
	out << "<atom><number>" << number << "</number><symbol>" << symbol << "</symbol><mass>" << mass << "</mass><charge>" << charge << 
		"</charge><x>" << x << "</x><y>" << y << "</y><z>" << z << "</z><forcex>" << forcex << "</forcex><forcey>" << forcey << 
		"</forcey><forcez>" << forcez << "</forcez></atom>" << std::endl;
}

void Atom::load(pugi::xml_node node){	
	number = atoi(node.child_value("number"));
	symbol = node.child_value("symbol");
	mass = atof(node.child_value("mass"));
	charge = atof(node.child_value("charge"));
	x = atof(node.child_value("x"));
	y = atof(node.child_value("y"));
	z = atof(node.child_value("z"));
	forcex = atof(node.child_value("forcex"));
	forcey = atof(node.child_value("forcey"));
	forcez = atof(node.child_value("forcez"));
}
