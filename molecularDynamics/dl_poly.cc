#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cstdlib>
#include <sstream>

#include "dl_poly.h"
	
std::string const Dl_Poly::name = "DL_POLY";

extern bool RemoveDirectory(string);

std::string const Dl_Poly::getName(){
	return "DL_POLY";
}

void Dl_Poly::generateConfig(std::string directory, Geometry geom){
	mkdir(directory.c_str(),0700);
	std::ofstream configfile;
	configfile.open((directory + "/CONFIG").c_str());
	configfile << "newffm Generated Original Geometry" << std::endl;
	
	configfile << "0          0" << std::endl;
	
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		configfile << " " << it->symbol << "             " << (it-geom.atoms.begin()+1) << std::endl;
		configfile << "  ";
		configfile.width(16);
		configfile.precision(7);
		configfile << it->x;
		configfile << "  ";
		configfile.width(16);
		configfile.precision(7);
		configfile << it->y;
		configfile << "  ";
		configfile.width(16);
		configfile.precision(7);
		configfile << it->z;
		configfile << "  " << std::endl;
	}

	configfile.close();
}

void Dl_Poly::generateConfig(std::string directory, Geometry geom, float transformationMatrix[][3]){
	mkdir(directory.c_str(),0700);
	std::ofstream configfile;
	configfile.open((directory + "/CONFIG").c_str());
	configfile << "newffm Generated Original Geometry" << std::endl;

	configfile << "0          3" << std::endl;
	configfile.precision(7);
	configfile.setf(std::ios_base::fixed,std::ios::floatfield);
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			configfile << "  ";
			configfile.width(16);
			configfile << transformationMatrix[i][j];
			configfile << "  ";
		}
		configfile << std::endl;
	}
	
	int i = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		configfile << " " << it->symbol << "             " << (i+1) << std::endl;
		configfile << "  ";
		configfile.width(16);
		configfile.precision(7);
		configfile << it->x;
		configfile << "  ";
		configfile.width(16);
		configfile.precision(7);
		configfile << it->y;
		configfile << "  ";
		configfile.width(16);
		configfile.precision(7);
		configfile << it->z;
		configfile << "  " << std::endl;

		i++;
	}

	configfile.close();
}

void Dl_Poly::generateField(std::string directory, Geometry geom){
	std::ofstream outField;
	std::ifstream inField(fieldFile.c_str());
	std::string line;
	
	std::cout << directory << std::endl;
	outField.open((directory + "/FIELD").c_str());
	while(std::getline(inField, line)){
		for(std::vector<double>::iterator it = variables.begin(); it != variables.end(); it++){
			std::string needle = "#_#";
			needle[1] = (it - variables.begin()) + 65;
			std::string::size_type pos = 0;
			while((pos = line.find(needle, pos)) != std::string::npos){
				std::ostringstream variable(std::ostringstream::out);
				variable << (*it);
				line.replace(pos, needle.size(), variable.str());
			}
		}
		outField << line << std::endl;
	}
	inField.close();
	outField.close();
}

void Dl_Poly::readForces(std::string directory){
	std::ifstream inFile((directory + "/OUTPUT").c_str());
	std::string line;
	float energy, temp, x, y, z;
	int step, atom;
	bool finished = false;
	
	while(std::getline(inFile, line)){
		if(line.find("run terminated") == std::string::npos){
			finished = true;
			break;
		}
	}

	if(!finished){
		MolecularDynamicsException excep("DL_POLY failed. Check dl_polyRun/ for errors and try again.");
		throw excep;
	}

	while(std::getline(inFile, line) && line.find("eng_tot") == std::string::npos);

	//Next four lines aren't what we want
	std::getline(inFile, line);
	std::getline(inFile, line);
	std::getline(inFile, line);
	std::getline(inFile, line);

	std::getline(inFile, line);
	std::istringstream iss(line);
	iss >> step >> energy;

	if(step != 1){
		MolecularDynamicsException excep("Couldn't find step 1 in OUTPUT. Check dl_polyRun/ for errors and try again.");
		throw excep;
	}
	
	energies.push_back(energy);

	while(std::getline(inFile, line) && line.find("sample of final configuration") == std::string::npos);

	//Next four lines aren't what we want
	std::getline(inFile, line);
	std::getline(inFile, line);
	std::getline(inFile, line);
	std::getline(inFile, line);

	std::vector< std::vector<float> > geom;
	while(std::getline(inFile, line)){
		std::istringstream iss2(line);

		if(!(iss2 >> atom))
			break;

		std::vector<float> atomForce;

		iss2 >> temp >> temp >> temp >> temp >> temp >> temp >> x >> y >> z;

		//		std::cout << x << ":" << y << ":" << z << std::endl;
		// Convert forces from DLPOLY OUTPUT file from internal units to kcalmol-1A-1
		x=x/418.4;
		y=y/418.4;
		z=z/418.4;
		//	        std::cout << x << ":" << y << ":" << z << std::endl;
		//		x=1/0;

		atomForce.push_back(x);
		atomForce.push_back(y);
		atomForce.push_back(z);

		geom.push_back(atomForce);
	}

	forces.push_back(geom);
	inFile.close();
}

std::vector<struct mdQuestion> Dl_Poly::getQuestions(){
	std::vector<struct mdQuestion> result;

	struct mdQuestion qControlFile;
	qControlFile.phrase = "CONTROL file";
	qControlFile.type = MDFILENAME;
	
	result.push_back(qControlFile);

	struct mdQuestion qFieldFile;
	qFieldFile.phrase = "FIELD file";
	qFieldFile.type = MDFILENAME;
	
	result.push_back(qFieldFile);

	struct mdQuestion qDlPolyExe;
	qDlPolyExe.phrase = "DL_POLY executable";
	qDlPolyExe.type = MDFILENAME;

	result.push_back(qDlPolyExe);

	return result;
}

void Dl_Poly::prepareMD(std::vector<GeometrySet> geoSets, std::vector<double> variables, std::vector<struct mdQuestion> questions){
	this->geoSets = geoSets;
	this->variables = variables;

	controlFile = questions[0].stringAnswer;
	fieldFile = questions[1].stringAnswer;
	dlPolyExe = questions[2].stringAnswer;

        RemoveDirectory("dl_polyRun");
	mkdir("dl_polyRun",0700);
}

void Dl_Poly::runMD(){
	std::string line;


	resetRun();
	
	int set = 1;
	int i;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ostringstream oss(std::ostringstream::out);
			oss << "dl_polyRun/set_" << set << "_geometry_" << i;
			mkdir(oss.str().c_str(),0700);

			generateField(oss.str(), *geomit);
			
			i++;
		}
		set++;
	}
	
	set = 1;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ostringstream oss(std::ostringstream::out);
			
			oss << "dl_polyRun/set_" << set << "_geometry_" << i;
			
			if(setit->a == 0 && setit->b == 0 && setit->c == 0) //We don't have a transformation matrix
				generateConfig(oss.str(), *geomit);
			else
				generateConfig(oss.str(), *geomit, setit->transformationMatrix);

			i++;
		}
		set++;
	}

	std::ifstream inControl;
	std::ofstream outControl;
	
	set = 1;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ifstream gInControl(controlFile.c_str());
			std::ostringstream oss(std::ostringstream::out);
			
			oss << "dl_polyRun/set_" << set << "_geometry_" << i << "/CONTROL";
			outControl.open(oss.str().c_str());
			while(std::getline(gInControl, line))
				outControl << line << std::endl;
			i++;
			gInControl.close();
			outControl.close();
		}
		set++;
	}

	set = 1;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ostringstream oss(std::ostringstream::out);
			
			oss << "dl_polyRun/set_" << set << "_geometry_" << i;
			chdir(oss.str().c_str());
			std::ostringstream execString;
			execString << dlPolyExe << " >& stdout";
			system(execString.str().c_str());
			chdir("..");
			chdir("..");
			readForces(oss.str());

			if(geomit->energy == 0.0) //Kind of a hack. Undoes reading energy from the output
				energies[energies.size()-1] = 0.0;
			
			i++;
		}
		set++;
	}

}
