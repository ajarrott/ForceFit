#ifndef _MOLECULARDYNAMICS_H
#define _MOLECULARDYNAMICS_H

#include <unistd.h>
#include <string.h>
#include "../geometrySet.h"

// Our code will provide the interface with a vector of questions
// The interface will then respond with answers by passing the same struct back
enum mdAnswerType { MDBOOL = 0, MDSTRING, MDFILENAME };

struct mdQuestion {
	enum mdAnswerType type;

	std::string phrase;

	// Boolean answers are filled in in the boolAnswer field
	bool boolAnswer;
	// String and filename answers are filled in in the stringAnswer field
	std::string stringAnswer;
};

class MolecularDynamics {
public:
	MolecularDynamics();
	virtual std::vector<struct mdQuestion> getQuestions() = 0;
	virtual void prepareMD(std::vector<GeometrySet> geoSets, std::vector<double> variables, std::vector<struct mdQuestion> questions) = 0;
	virtual void runMD() = 0;
	virtual int getIterations();
	virtual std::vector< std::vector< std::vector<float> > > getForces();
	virtual std::vector<float> getEnergies();
	const std::vector<double> & getVariables();
	void setVariables(const std::vector<double> & values);
	static std::string const name;
	virtual std::string const getName() = 0;

protected:
	std::vector<GeometrySet> geoSets;
	void resetRun();
	std::string mdPath;
	std::vector<double> variables;
	int iterations;
	std::vector< std::vector< std::vector<float> > > forces;
	std::vector<float> energies;
};

class MolecularDynamicsException : public std::exception {
public:
	MolecularDynamicsException(const char* message);

	virtual const char * what() const throw();
private:
	const char* whatMessage;
};

#endif
