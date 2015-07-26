#include "geometrySet.h"

#include <sstream>
#include <cstdlib>

GeometrySet::GeometrySet(){
	a = 0;
	b = 0;
	c = 0;
	alpha = 0;
	beta = 0;
	gamma = 0;
	volume = 0;
	name = "";
}

void GeometrySet::write(std::ostream & out, int indent){
	for(int i = 0; i < indent; i++){
		out << "\t";
	}
	out << "<geometrySet>" << std::endl;
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "<name>" << name << "</name>" << std::endl;
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "<a>" << a << "</a><b>" << b << "</b><c>" << c << "</c><alpha>" << alpha << "</alpha><beta>" << beta << "</beta><gamma>" << gamma << "</gamma><volume>" << volume << "</volume><energyZero>" << energyZero << "</energyZero>" << std::endl;
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "<transformationMatrix>";
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			out << transformationMatrix[i][j];
			if(i != 2 || j != 2)
				out << ",";
		}
	}
	out << "</transformationMatrix>" << std::endl;
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "<frequencies>" << std::endl;
	for(std::vector<float>::iterator it = frequencies.begin(); it != frequencies.end(); it++){
		for(int i = 0; i < indent+2; i++){
			out << "\t";
		}
		out << "<frequency>" << (*it) << "</frequency>" << std::endl;
	}
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "</frequencies>" << std::endl;
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "<normals>" << std::endl;
	for(std::vector< std::vector< std::vector<float> > >::iterator xit = normals.begin(); xit != normals.end(); xit++){
		for(int i = 0; i < indent+2; i++){
			out << "\t";
		}
		out << "<xrow>" << std::endl;
		for(std::vector< std::vector<float> >::iterator yit = xit->begin(); yit != xit->end(); yit++){
			for(int i = 0; i < indent+3; i++){
				out << "\t";
			}
			out << "<yrow>" << std::endl;
			for(std::vector<float>::iterator zit = yit->begin(); zit != yit->end(); zit++){
				for(int i = 0; i < indent+4; i++){
					out << "\t";
				}
				out << "<normal>" << *zit << "</normal>" << std::endl;
			}
			for(int i = 0; i < indent+3; i++){
				out << "\t";
			}
			out << "</yrow>" << std::endl;
		}
		for(int i = 0; i < indent+2; i++){
			out << "\t";
		}
		out << "</xrow>" << std::endl;
	}
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "</normals>" << std::endl;
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}

	out << "<geometries>" << std::endl;
	for(std::vector<Geometry>::iterator it = geometries.begin(); it != geometries.end(); it++){
		it->write(out, indent+2);
	}
	for(int i = 0; i < indent+1; i++){
		out << "\t";
	}
	out << "</geometries>" << std::endl;
	
	for(int i = 0; i < indent; i++){
		out << "\t";
	}
	out << "</geometrySet>" << std::endl;
}

void GeometrySet::load(pugi::xml_node node){
	name = node.child_value("name");
	a = atof(node.child_value("a"));
	b = atof(node.child_value("b"));
	c = atof(node.child_value("c"));
	alpha = atof(node.child_value("alpha"));
	beta = atof(node.child_value("beta"));
	gamma = atof(node.child_value("gamma"));
	volume = atof(node.child_value("volume"));
	energyZero = atof(node.child_value("energyZero"));

	{
		std::istringstream iss(node.child_value("transformationMatrix"));
		char buffer[256];
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				iss.getline(buffer, 255, ',');
				transformationMatrix[i][j] = atof(buffer);
			}
		}
	}

	frequencies.clear();
	for(pugi::xml_node frequency = node.child("frequencies").child("frequency"); frequency; frequency = frequency.next_sibling("frequency")){
		frequencies.push_back(atof(frequency.child_value()));
	}

	normals.clear();
	for(pugi::xml_node xrow = node.child("normals").child("xrow"); xrow; xrow = xrow.next_sibling("xrow")){
		std::vector< std::vector<float> > xrowVec;
		for(pugi::xml_node yrow = xrow.child("yrow"); yrow; yrow = yrow.next_sibling("yrow")){
			std::vector<float> yrowVec;
			for(pugi::xml_node normal = yrow.child("normal"); normal; normal = normal.next_sibling("normal")){
				yrowVec.push_back(atof(normal.child_value()));
			}
			xrowVec.push_back(yrowVec);
		}
		normals.push_back(xrowVec);
	}

	geometries.clear();
	for(pugi::xml_node geometryNode = node.child("geometries").child("geometry"); geometryNode; geometryNode = geometryNode.next_sibling("geometry")){
		Geometry geometry;
		geometry.load(geometryNode);
		geometries.push_back(geometry);
	}
}
