#include "crystalScan.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#define PI 3.14159265

std::string const CrystalScan::name = "Crystal Output Type";

std::string const CrystalScan::getName(){
	return name;
}

GeometrySet CrystalScan::getGeometries(){
	
	/*
	   Data to be gathered
	
	float a,b,c,alpha,beta,gamma,volume;
	float transformationMatrix[3][3];
	//Normals indexed by: Mode, Atom, Axis (0-X 1-Y 2-Z)
	std::vector< std::vector <std::vector<float> > > normals;
	std::vector<float> frequencies;
	std::vector<Atom> originalAtoms;
	std::vector<Geometry> geometries;
	float energyZero;
	*/

	GeometrySet geoSet;
	geoSet.name = "Crystal Set from ";
	bool comma = false;
	
	for(std::vector<std::string>::iterator it = inputFiles.begin(); it != inputFiles.end(); it++){
		if(comma)
			geoSet.name += ", ";
		comma = true;
		geoSet.name += *it;
		
		std::string line;
		std::ifstream inputFile(it->c_str());
		while(!inputFile.eof()){
			std::getline(inputFile, line);
			bool points = false;
	
			std::vector<float> newFreqs;

			if(line.find("LATTICE PARAMETERS  (ANGSTROMS AND DEGREES) - PRIMITIVE CELL") != std::string::npos){
				//Skip "A B C Alpha Beta Gamma Volume"
				std::getline(inputFile, line);
				
				//Reading and parsing A, B, C, Alpha, Beta, Gamma
				std::getline(inputFile, line);
				std::istringstream iss(line);
				iss >> geoSet.a >> geoSet.b >> geoSet.c >> geoSet.alpha >> geoSet.beta >> geoSet.gamma >> geoSet.volume;
			}

			if(line.find("COORDINATES OF THE EQUIVALENT ATOMS (FRACTIONARY UNITS)") != std::string::npos){
				while(!inputFile.eof()){
					std::getline(inputFile, line);

					if(line.length() == 0)
						continue;

					int num, atomicNum;
					std::string symbol;
					float x,y,z;
					std::string temp;

					std::istringstream iss(line);

					iss >> symbol;

					if(symbol == "N.")
						continue;

					if(atoi(symbol.c_str()) == 0)
						break;

					iss >> num >> temp >> atomicNum >> symbol >> x >> y >> z;

					Geometry newGeom;

					if(newGeom.atoms.size() < num){
						Atom newAtom;
						newAtom.number = atomicNum;
						newAtom.symbol = symbol;
						newAtom.x = x;
						newAtom.y = y;
						newAtom.z = z;
						newGeom.atoms.push_back(newAtom);
						newGeom.mode = 0;
						newGeom.step = 0;
					} else {
						newGeom.atoms[num-1].number = atomicNum;
						newGeom.atoms[num-1].symbol = symbol;
						newGeom.atoms[num-1].x = x;
						newGeom.atoms[num-1].y = y;
						newGeom.atoms[num-1].z = z;
					}
					
					geoSet.geometries.push_back(newGeom);
				}
			}

			//Finding, reading, and parsing the input transformation matrix
			//Note that we generate our own transformation matrix later, in order to perform
			//   the converion from fractional to cartesian coordinates
			if(line.find("DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)") != std::string::npos){
				//Skip "X  Y  Z"
				std::getline(inputFile, line);
				
				for(int i = 0; i < 3; i++){
					std::getline(inputFile, line);
					std::istringstream iss(line);
					iss >> geoSet.transformationMatrix[i][0] >> geoSet.transformationMatrix[i][1] >> geoSet.transformationMatrix[i][2];
				}
			}
			
			if(line.find("NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES") != std::string::npos){
				int currentAtom = -1;
				int currentMode = 0;
				int startMode = 0;
				std::string temp;

				while(!inputFile.eof()){
					std::getline(inputFile, line);

					if(line.find("*******************************************************************************") != std::string::npos)
						break;

					if(line.length() == 0)
						continue;

					std::istringstream iss(line);
					iss >> temp;

					if(temp == "FREQ(CM**-1)"){
						startMode = currentMode;
						float num;

						while(!iss.eof()){
							iss >> num;
							newFreqs.push_back(num);
						}
					}

					if(temp == "AT." || temp == "Y" || temp == "Z"){
						currentMode = startMode;
						
						if(temp == "AT."){
							iss >> currentAtom;
							currentAtom--;
							//Symbol
							iss >> temp;
							//X
							iss >> temp;
						}

						float num;
						while(!iss.eof()){
							while(geoSet.normals.size() <= currentMode){
								std::vector<std::vector<float> > newVec;
								geoSet.normals.push_back(newVec);
							}

							while(geoSet.normals[currentMode].size() <= currentAtom){
								std::vector<float> newVec;
								geoSet.normals[currentMode].push_back(newVec);
							}

							iss >> num;
							geoSet.normals[currentMode][currentAtom].push_back(num);

							currentMode++;
						}
					}
				}
			}
			
			if(line.find("ATOMS ISOTOPIC MASS (AMU) FOR FREQUENCY CALCULATION") != std::string::npos){
				std::string line;
				while(!inputFile.eof()){
					std::getline(inputFile, line);

					if(line.length() == 0)
						continue;

					int num;
					std::string symbol;
					float mass;

					std::istringstream iss(line);
					while(!iss.eof()){
						iss >> num >> symbol >> mass;
						if(num < 0 || symbol.length() > 2 || symbol.length() == 0)
							break;

						if(geoSet.geometries[0].atoms.size() < num){
							Atom newAtom;
							newAtom.symbol = symbol;
							newAtom.mass = mass;
							geoSet.geometries[0].atoms.push_back(newAtom);
						} else {
							geoSet.geometries[0].atoms[num-1].symbol = symbol;
							geoSet.geometries[0].atoms[num-1].mass = mass;
						}
					}

					if(num < 0 || symbol.length() > 2 || symbol.length() == 0)
						break;
				}
			}

			if(line.find("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS") != std::string::npos){
				int currentMode;
				float currentStep;
				int startLineCount = 0;
				std::string temp;

				while(std::getline(inputFile, line)){
					if(line.length() == 0)
						continue;

					//We're in the energies file! Let the next if statement catch this.
					if(line.find("MODE(CM-1) DISPLAC   TOTAL ENE(DFT)(AU)  CLASSICAL HARM ENE(AU)  NCYC    DE") != std::string::npos){
						break;
					}
					
					//This is our first hint at being in the points file and NOT in the energies file.
					points = true;

					std::istringstream iss(line);
					iss >> temp;

					if(temp == "MODE"){
						iss >> currentMode;
					}

					if(temp == "GEOMETRY"){
						startLineCount = 0;
						iss >> temp >> temp >> currentStep;
					}

					if(temp == "*******************************************************************************"){
						startLineCount++;

						if(startLineCount == 3){
							Geometry newGeom;
							newGeom.mode = currentMode;
							newGeom.step = currentStep;
							while(std::getline(inputFile, line)){
								if(line.length() == 0)
									break;

								std::istringstream iss2(line);

								iss2 >> temp;

								//Make sure the first word in the line is a positive number
								if(atoi(temp.c_str()) <= 0)
										break;

								Atom newAtom;

								//Skip the 'T' and the number
								iss2 >> temp >> temp >> newAtom.symbol >> newAtom.x >> newAtom.y >> newAtom.z;

								newAtom.mass = geoSet.geometries[0].atoms[newGeom.atoms.size()].mass;
								newGeom.atoms.push_back(newAtom);
							}

							geoSet.geometries.push_back(newGeom);
						}
					}
				}
			}
			
			if(line.find("MODE(CM-1) DISPLAC   TOTAL ENE(DFT)(AU)  CLASSICAL HARM ENE(AU)  NCYC    DE") != std::string::npos) {
				int currentMode = 0;

				while(!inputFile.eof()){
					std::getline(inputFile, line);

					if(line.length() == 0)
						break;

					if(line.find_first_of('(') != std::string::npos){
						std::istringstream iss(line.substr(0,line.find_first_of('(')));
						iss >> currentMode;
					} else {
						std::istringstream iss(line);

						float step, energy;
						std::string temp;

						iss >> step >> temp >> energy;

						for(std::vector<Geometry>::iterator it = geoSet.geometries.begin(); it != geoSet.geometries.end(); it++){
							if(it->mode == currentMode && it->step == step)
								it->energy = energy;
						}

						if(step == 0.0){
							geoSet.energyZero = energy;
						}
					}
				}
			}

			//The frequencies that we get from the energies file are wrong.
			if(points)
				geoSet.frequencies = newFreqs;
		}
		inputFile.close();
	}

	
	float realTransformationMatrix[3][3];
	realTransformationMatrix[0][0] = geoSet.a;
	realTransformationMatrix[0][1] = geoSet.b*cos(geoSet.gamma*PI/180);
	realTransformationMatrix[0][2] = geoSet.c*cos(geoSet.beta*PI/180);
	realTransformationMatrix[1][0] = 0.0;
	realTransformationMatrix[1][1] = geoSet.b*sin(geoSet.gamma*PI/180);
	realTransformationMatrix[1][2] = geoSet.c*(cos(geoSet.alpha*PI/180)-cos(geoSet.beta*PI/180)*cos(geoSet.gamma*PI/180))/sin(geoSet.gamma*PI/180);
	realTransformationMatrix[2][0] = 0.0;
	realTransformationMatrix[2][1] = 0.0;
	realTransformationMatrix[2][2] = geoSet.volume / (geoSet.a*geoSet.b*sin(geoSet.gamma*PI/180));

	float determinant = realTransformationMatrix[0][0]*(realTransformationMatrix[2][2]*realTransformationMatrix[1][1] - realTransformationMatrix[2][1]*realTransformationMatrix[1][2]) -
		realTransformationMatrix[1][0]*(realTransformationMatrix[2][2]*realTransformationMatrix[0][1] - realTransformationMatrix[2][1]*realTransformationMatrix[0][2]) +
		realTransformationMatrix[2][0]*(realTransformationMatrix[1][2]*realTransformationMatrix[0][1] - realTransformationMatrix[1][1]*realTransformationMatrix[0][2]);

	float inverseTransformationMatrix[3][3];

	inverseTransformationMatrix[0][0] = 1/determinant*(realTransformationMatrix[2][2]*realTransformationMatrix[1][1] - realTransformationMatrix[2][1]*realTransformationMatrix[1][2]);
	inverseTransformationMatrix[0][1] = 1/determinant*(0-(realTransformationMatrix[2][2]*realTransformationMatrix[0][1] - realTransformationMatrix[2][1]*realTransformationMatrix[0][2]));
	inverseTransformationMatrix[0][2] = 1/determinant*(realTransformationMatrix[1][2]*realTransformationMatrix[0][1] - realTransformationMatrix[1][1]*realTransformationMatrix[0][2]);
	inverseTransformationMatrix[1][0] = 1/determinant*(0-(realTransformationMatrix[2][2]*realTransformationMatrix[1][0] - realTransformationMatrix[2][0]*realTransformationMatrix[1][2]));
	inverseTransformationMatrix[1][1] = 1/determinant*(realTransformationMatrix[2][2]*realTransformationMatrix[0][0] - realTransformationMatrix[2][0]*realTransformationMatrix[0][2]);
	inverseTransformationMatrix[1][2] = 1/determinant*(0-(realTransformationMatrix[1][2]*realTransformationMatrix[0][0] - realTransformationMatrix[1][0]*realTransformationMatrix[0][2]));
	inverseTransformationMatrix[2][0] = 1/determinant*(realTransformationMatrix[2][1]*realTransformationMatrix[1][0] - realTransformationMatrix[2][0]*realTransformationMatrix[1][1]);
	inverseTransformationMatrix[2][1] = 1/determinant*(0-(realTransformationMatrix[2][1]*realTransformationMatrix[0][0] - realTransformationMatrix[2][0]*realTransformationMatrix[0][1]));
	inverseTransformationMatrix[2][2] = 1/determinant*(realTransformationMatrix[1][1]*realTransformationMatrix[0][0] - realTransformationMatrix[1][0]*realTransformationMatrix[0][1]);

	for(std::vector<Atom>::iterator it = geoSet.geometries[0].atoms.begin(); it != geoSet.geometries[0].atoms.end(); it++){
		if(it->x < 0.0)
			it->x += 1.0;

		if(it->y < 0.0)
			it->y += 1.0;

		if(it->z < 0.0)
			it->z += 1.0;
	}

	for(std::vector<Geometry>::iterator geomit = geoSet.geometries.begin(); geomit != geoSet.geometries.end(); geomit++){
		int atom = 0;
		for(std::vector<Atom>::iterator atomit = geomit->atoms.begin(); atomit != geomit->atoms.end(); atomit++){
			if(atomit->x - geoSet.geometries[0].atoms[atom].x > 0.5)
				atomit->x -= 1.0;

			if(atomit->x - geoSet.geometries[0].atoms[atom].x < -0.5)
				atomit->x += 1.0;

			if(atomit->y - geoSet.geometries[0].atoms[atom].y > 0.5)
				atomit->y -= 1.0;

			if(atomit->y - geoSet.geometries[0].atoms[atom].y < -0.5)
				atomit->y += 1.0;

			if(atomit->z - geoSet.geometries[0].atoms[atom].z > 0.5)
				atomit->z -= 1.0;

			if(atomit->z - geoSet.geometries[0].atoms[atom].z < -0.5)
				atomit->z += 1.0;

			float x, y, z;
			
			x = atomit->x * realTransformationMatrix[0][0] + atomit->y * realTransformationMatrix[0][1] + atomit->z * realTransformationMatrix[0][2];
			y = atomit->x * realTransformationMatrix[1][0] + atomit->y * realTransformationMatrix[1][1] + atomit->z * realTransformationMatrix[1][2];
			z = atomit->x * realTransformationMatrix[2][0] + atomit->y * realTransformationMatrix[2][1] + atomit->z * realTransformationMatrix[2][2];

			atomit->x = x;
			atomit->y = y;
			atomit->z = z;

			atom++;
		}
	}

	float x, y, z;
	
	for(std::vector<Atom>::iterator it = geoSet.geometries[0].atoms.begin(); it != geoSet.geometries[0].atoms.end(); it++){
		x = it->x * realTransformationMatrix[0][0] + it->y * realTransformationMatrix[0][1] + it->z * realTransformationMatrix[0][2];
		y = it->x * realTransformationMatrix[1][0] + it->y * realTransformationMatrix[1][1] + it->z * realTransformationMatrix[1][2];
		z = it->x * realTransformationMatrix[2][0] + it->y * realTransformationMatrix[2][1] + it->z * realTransformationMatrix[2][2];
		
		it->x = x;
		it->y = y;
		it->z = z;
	}

	return geoSet;
}
