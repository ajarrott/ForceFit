#ifndef _POWELL1965_H
#define _POWELL1965_H

#include "minimizer.h"
#include "../molecularDynamics/molecularDynamics.h"

class Powell1965 : public Minimizer {
public:
	virtual std::vector<struct minQuestion> getQuestions();
	virtual void prepareMin(const std::vector<GeometrySet> & geoSets, std::vector<struct minQuestion> questions);
	virtual void runMin(MolecularDynamics * md);
	static std::string const name;
	virtual std::string const getName();
private:
	std::vector<double> calculateUForce(struct dataPoint data);
	std::vector<double> calculateUEnergy(struct dataPoint data);
	struct dataPoint runValues(MolecularDynamics * md, std::vector<double> values);

	MolecularDynamics * md;
	
	double linmin(std::vector<double> & p, std::vector<double> & xi);
	double fidim(std::vector<double> & p, std::vector<double> & xi, double x);
	double brent(double ax, double bx, double cx, std::vector<double> & p, std::vector<double> & xi, double tol, double *xmin);
	void mnbrak(double * ax, double * bx, double * cx, double * fa, double * fb, double *fc, std::vector<double> & p, std::vector<double> & xi);

	double epsilon;
	bool useForce, useEnergy;
	std::vector<double> finalVariables;
};

#endif
