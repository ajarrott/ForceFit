#ifndef _GEOMETRYSET_H
#define _GEOMETRYSET_H

#include "pugixml/pugixml.hpp"

#include "geometry.h"

class GeometrySet {
public:
	GeometrySet();
	
	std::vector<Geometry> geometries;

	std::vector<float> frequencies;
	std::vector< std::vector <std::vector<float> > > normals;
	std::string name;
	float transformationMatrix[3][3];
	float a,b,c,alpha,beta,gamma,volume;
	float energyZero;

	void write(std::ostream &, int indent);
	void load(pugi::xml_node node);
};

#endif
