#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <sstream>

#include "lammps.h"

Lammps::Lammps(){
}
	
std::string const Lammps::name = "LAMMPS";

extern bool RemoveDirectory(string);

std::string const Lammps::getName(){
	return "LAMMPS";
}

void Lammps::generateInput(Geometry geom, std::string directory){
	mkdir(directory.c_str(),0700);
	std::ofstream outInput;
	outInput.setf(std::ios_base::fixed,std::ios::floatfield);
	outInput.open((directory + "/input").c_str());

	if(directory[0] != '/'){
		char cwd[256];
		getcwd(cwd, 256);
		directory.insert(0, "/");
		directory.insert(0, cwd);
	}
	
	outInput << "shell           cd " << directory << std::endl;
	outInput << "" << std::endl;
	
	std::ifstream inInput(inputFile.c_str());
	std::string line;
	
	while(std::getline(inInput, line)){
		if(line.find("read_data") != string::npos){
			outInput << "read_data   data" << std::endl;
			continue;
		} else if(line.find("dump") != string::npos){
			outInput << "dump            1 all custom 1 geometry.out fx fy fz" << std::endl;
			continue;
		}

		for(std::vector<double>::iterator it = variables.begin(); it != variables.end(); it++){
			std::string needle = "#_#";
			needle[1] = (it - variables.begin()) + 65;
			std::string::size_type pos = 0;
			while((pos = line.find(needle, pos)) != std::string::npos){
				std::ostringstream variable(std::ostringstream::out);
				variable.width(7);
				variable << (*it);
				line.replace(pos, needle.size(), variable.str());
			}
		}
		outInput << line << std::endl;
	}


	inInput.close();
	outInput.close();
}

void Lammps::generateData(Geometry geom, std::string directory){
	mkdir(directory.c_str(),0700);
	std::ofstream outData((directory + "/data").c_str());
	outData.setf(std::ios_base::fixed,std::ios::floatfield);
	std::ifstream inData(dataFile.c_str());
	std::string line;
	
	while(std::getline(inData, line)){
		if(line.find("Atoms") != string::npos){
			outData << line;
			outData << std::endl;
			outData << std::endl;

			std::vector<std::string> types;

			bool found;
			int hydrogens = 0;
			for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
				found = false;
				
				if(it->symbol == "H")
					hydrogens++;
				
				for(std::vector<std::string>::iterator findit = types.begin(); findit != types.end(); findit++){
					if(*findit == it->symbol){
						found = true;
						break;
					}
				}
				
				if(!found){
					types.push_back(it->symbol);
				}
			}

			int i = 0;
			for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
				outData.width(7);
				outData << i+1;
				outData.width(7);
				outData << 0;
				outData.width(4);
				int typeNum = 1;
				for(std::vector<std::string>::iterator typeit = types.begin(); typeit != types.end(); typeit++){
					if(*typeit == it->symbol){
						outData << typeNum;
						break;
					}
					typeNum++;
				}
				outData.width(10);
				outData.precision(6);
				outData << it->charge;
				
				outData.width(16);
				outData.precision(9);
				outData << it->x;
				outData << "  ";
				outData.width(16);
				outData.precision(9);
				outData << it->y;
				outData << "  ";
				outData.width(16);
				outData.precision(9);
				outData << it->z;
				outData << "  " << std::endl;

				i++;
			}

			outData << std::endl;

			while(std::getline(inData, line)){
				float temp;
				std::istringstream iss(line);
				if(line.length() != 0 && !(iss >> temp))
					break;
			}
		}

		for(std::vector<double>::iterator it = variables.begin(); it != variables.end(); it++){
			std::string needle = "#_#";
			needle[1] = (it - variables.begin()) + 65;
			std::string::size_type pos = 0;
			while((pos = line.find(needle, pos)) != std::string::npos){
				std::ostringstream variable(std::ostringstream::out);
				variable.width(7);
				variable << (*it);
				line.replace(pos, needle.size(), variable.str());
			}
		}
		outData << line << std::endl;
	}

	inData.close();
	outData.close();
}

std::vector<struct mdQuestion> Lammps::getQuestions(){
	std::vector<struct mdQuestion> result;

	struct mdQuestion qMixed;

	struct mdQuestion qInputPath;
	qInputPath.phrase = "Input File";
	qInputPath.type = MDFILENAME;

	result.push_back(qInputPath);

	struct mdQuestion qDataPath;
	qDataPath.phrase = "Data File";
	qDataPath.type = MDFILENAME;

	result.push_back(qDataPath);

	struct mdQuestion qLammpsPath;
	qLammpsPath.phrase = "Lammps Executable";
	qLammpsPath.type = MDFILENAME;

	result.push_back(qLammpsPath);

	return result;
}

void Lammps::prepareMD(std::vector<GeometrySet> geoSets, std::vector<double> variables, std::vector<struct mdQuestion> questions){
	this->geoSets = geoSets;
	this->variables = variables;

        RemoveDirectory("lammpsRun");
	mkdir("lammpsRun",0700);

	std::string line;
	
	inputFile = questions[0].stringAnswer;
	dataFile = questions[1].stringAnswer;
	mdPath = questions[2].stringAnswer;
}

void Lammps::readForces(std::string directory){
	std::ifstream inFile((directory + "/output").c_str());
	std::string line;
	float energy, temp, x, y, z;
	int step, atom;
	bool finished = false;
	
	while(std::getline(inFile, line)){
		if(line.find("TotEng") != std::string::npos){
			finished = true;
			break;
		}
	}

	if(!finished){
		MolecularDynamicsException excep("Lammps failed. Check lammpsRun/ for errors and try again.");
		throw excep;
	}

	//Skip step 0
	std::getline(inFile, line);

	std::getline(inFile, line);
	std::istringstream iss(line);
	iss >> step >> temp >> temp >> temp >> energy;

	if(step != 1){
		MolecularDynamicsException excep("Couldn't find step 1 in output. Check lammpsRun/ for errors and try again.");
		throw excep;
	}
	
	energies.push_back(energy);

	inFile.close();

	std::ifstream in2File((directory + "/geometry.out").c_str());

	step = -1;

	while(step != 1){
		while(std::getline(in2File, line) && line.find("TIMESTEP") == std::string::npos);

		std::getline(in2File, line);

		std::istringstream iss2(line);
		iss2 >> step;

		if(step != 1)
			continue;
		
		while(std::getline(in2File, line) && line.find("ITEM: ATOMS") == std::string::npos);

		std::vector< std::vector<float> > geom;
		while(std::getline(in2File, line)){
			std::istringstream iss2(line);

			std::vector<float> atomForce;

			if(!(iss2 >> x >> y >> z))
					break;

			std::cout << x << " " << y << " " << z << std::endl;

			atomForce.push_back(x);
			atomForce.push_back(y);
			atomForce.push_back(z);

			geom.push_back(atomForce);
		}

		forces.push_back(geom);

		break;
	}
	in2File.close();
}

void Lammps::runMD(){
	int i;
	std::string line;

	resetRun();
	
	int set = 1;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ostringstream oss(std::ostringstream::out);
			oss << "lammpsRun/set_" << set << "_geometry_" << i;
			generateData(*geomit, oss.str());
			generateInput(*geomit, oss.str());
			
			system((mdPath + " < " + oss.str() + "/input >& " + oss.str() + "/output").c_str());
			readForces(oss.str());
			i++;
		}
		set++;
	}
}
