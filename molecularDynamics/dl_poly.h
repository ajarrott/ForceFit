#ifndef _DL_POLY_H
#define _DL_POLY_H

#include "molecularDynamics.h"

struct potential {
	std::string atom1, atom2, type;
	std::vector<std::string> parameters;
};

class Dl_Poly : public MolecularDynamics {
public:
	virtual std::vector<struct mdQuestion> getQuestions();
	virtual void prepareMD(std::vector<GeometrySet> geoSets, std::vector<double> variables, std::vector<struct mdQuestion> questions);
	virtual void runMD();
	static std::string const name;
	virtual std::string const getName();
private:
	void generateField(std::string directory, Geometry geom);
	void generateConfig(std::string directory, Geometry geom);
	void generateConfig(std::string directory, Geometry geom, float transformationMatrix[][3]);
	std::string controlFile, fieldFile, dlPolyExe;
	void readForces(std::string);
};

#endif
