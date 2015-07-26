#include <cstdlib>

#include "geometry.h"

Geometry::Geometry(){
	mode = 0;
	step = 0;
	energy = 0;
}

void Geometry::write(std::ostream & out, int indent){
	for(int i = 0; i < indent; i++){
		out << "\t";
	}
	out << "<geometry>" << std::endl;
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "<mode>" << mode << "</mode><step>" << step << "</step><energy>" << energy << "</energy>" << std::endl;
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "<atoms>" << std::endl;
	for(std::vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); it++){
		it->write(out, indent+2);
	}
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "</atoms>" << std::endl;
	for(int i = 0; i < indent; i++){
		out << "\t";
	}
	out << "</geometry>" << std::endl;
}

void Geometry::load(pugi::xml_node node){
	mode = atof(node.child_value("mode"));
	step = atof(node.child_value("step"));
	energy = atof(node.child_value("energy"));
	
	atoms.clear();
	for(pugi::xml_node atomNode = node.child("atoms").child("atom"); atomNode; atomNode = atomNode.next_sibling("atom")){
		Atom atom;
		atom.load(atomNode);
		atoms.push_back(atom);
	}
}
