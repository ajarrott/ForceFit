#ifndef _TINKER_H
#define _TINKER_H

#include "molecularDynamics.h"

class Tinker : public MolecularDynamics {
public:
	virtual std::vector<struct mdQuestion> getQuestions();
	virtual void prepareMD(std::vector<GeometrySet> geoSets, std::vector<double> variables, std::vector<struct mdQuestion> questions);
	virtual void runMD();
	static std::string const name;
	virtual std::string const getName();
private:
	std::string xyzFile, prmFile, keyFile, tinkerExe;

	void generateXyz(std::string directory, Geometry geom);
	void generatePrm(std::string directory, Geometry geom);
	void readForces(std::string);
};
        bool RemoveDirectory(string);

#endif
