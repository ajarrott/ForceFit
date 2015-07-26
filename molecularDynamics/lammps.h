#ifndef _LAMMPS_H
#define _LAMMPS_H

#include "molecularDynamics.h"

class Lammps : public MolecularDynamics {
public:
	Lammps();
	virtual std::vector<struct mdQuestion> getQuestions();
	virtual void prepareMD(std::vector<GeometrySet> geoSets, std::vector<double> variables, std::vector<struct mdQuestion> questions);
	virtual void runMD();
	static std::string const name;
	virtual std::string const getName();
private:
	void generateData(Geometry geom, std::string directory);
	void generateInput(Geometry geom, std::string directory);
	void checkForError(std::string filename);
	void readForces(std::string);

	std::string inputFile;
	std::string dataFile;
};

#endif
