#ifndef _MINIMIZER_H
#define _MINIMIZER_H

#include "../geometrySet.h"
#include "../molecularDynamics/molecularDynamics.h"
#include <stdlib.h>
#include <iterator>
#include <algorithm>

// Our code will provide the interface with a vector of questions
// The interface will then respond with answers by passing the same struct back
enum minAnswerType { MINBOOL = 0, MINSTRING, MINFILENAME };

struct minQuestion {
	enum minAnswerType type;

	std::string phrase;

	// Boolean answers are filled in in the boolAnswer field
	bool boolAnswer;
	// String and filename answers are filled in in the stringAnswer field
	std::string stringAnswer;
};

class Minimizer {
public:
	virtual std::vector<struct minQuestion> getQuestions() = 0;
	virtual void prepareMin(const std::vector<GeometrySet> & geoSets, std::vector<struct minQuestion> questions) = 0;
	virtual void runMin(MolecularDynamics * md) = 0;
	static std::string const name;
	virtual std::string const getName() = 0;

protected:
	std::vector<GeometrySet> geoSets;
};

#endif
