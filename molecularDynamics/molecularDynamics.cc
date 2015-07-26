#include <iostream>
#include <fstream>
#include <sstream>

#include "molecularDynamics.h"

std::string const MolecularDynamics::name = "Generic Molecular Dynamics";

MolecularDynamics::MolecularDynamics() : iterations(0) {
}

int MolecularDynamics::getIterations(){
	return iterations;
}

std::vector< std::vector< std::vector<float> > > MolecularDynamics::getForces(){
	return forces;
}

std::vector<float> MolecularDynamics::getEnergies(){
	return energies;
}

const std::vector<double> & MolecularDynamics::getVariables(){
	return variables;
}

void MolecularDynamics::setVariables(const std::vector<double> & values){
	if(variables.size() != values.size())
		return;
	
	variables.clear();
	for(std::vector<double>::const_iterator it = values.begin(); it != values.end(); it++)
		variables.push_back(*it);
}

void MolecularDynamics::resetRun(){
	iterations++;
	forces.clear();
	energies.clear();
}

MolecularDynamicsException::MolecularDynamicsException(const char * message){
	whatMessage = message;
}

const char * MolecularDynamicsException::what() const throw(){
	return whatMessage;
}
