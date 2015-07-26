#include "gaussianScan.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#define PI 3.14159265

std::string const GaussianScan::name = "Gaussian Output Type";

std::string const GaussianScan::getName(){
	return name;
}

GeometrySet GaussianScan::getGeometries(){
	GeometrySet geoSet;
	
	double stepSize;
	int originalCount = 0;

	geoSet.name = "Gaussian Set from ";

	bool comma = false;
	for(std::vector<std::string>::iterator it = inputFiles.begin(); it != inputFiles.end(); it++){
		if(comma)
			geoSet.name += ", ";
		comma = true;
		geoSet.name += *it;
		
		std::ifstream inputFile(it->c_str());

		bool nosymmFound = false, hpmodesFound = false, freqFound = false;

		bool needForces = true;

		int nextDisplacement = 0;
		std::string line;

		while(std::getline(inputFile, line)){
			if(line.find(" # ") == 0 || line.find("# ") == 0){
				while(line.find("--------")  == std::string::npos){
					std::transform(line.begin(), line.end(), line.begin(), ::tolower);
					if(line.find("nosymm") != std::string::npos)
						nosymmFound = true;

					if(line.find("hpmodes") != std::string::npos)
						hpmodesFound = true;

					if(line.find("freq") != std::string::npos)
						freqFound = true;

					if(line.find("opt") != std::string::npos)
						needForces = false;

					if(line.find("--------") != std::string::npos){
						//if(!hpmodesFound)
						//	std::cout << "Warning: HPMODES not found." << std::endl;

						//if(!nosymmFound)
						//	std::cout << "Error: nosymm not found." << std::endl;

						if(!freqFound)
							std::cout << "Error: FREQ not found." << std::endl;

						if(!needForces)
							std::cout << "Error: opt found." << std::endl;

						if(!freqFound || !needForces){
							ScanFileFormatException excep("Too many errors. Aborting.");
							throw excep;
						}

						break;
					}
					if(!std::getline(inputFile, line)){
						ScanFileFormatException excep("Couldn't read past keywords. Aborting.");
						throw excep;
					}
				} 
			}
			
			if(line.find("Input orientation:") != std::string::npos){
				originalCount++;
				Geometry newGeom;
				// Skip "--------------------------------------------------------------------"
				std::getline(inputFile, line);
				// Skip "Center     Atomic     Atomic              Coordinates (Angstroms)"
				std::getline(inputFile, line);
				// Skip "Number     Number      Type              X           Y           Z"
				std::getline(inputFile, line);
				// Skip "--------------------------------------------------------------------"
				std::getline(inputFile, line);

				while(std::getline(inputFile, line)){
					Atom newAtom;
					std::string temp;
					int num;

					std::istringstream iss(line);

					if(!(iss >> num))
						break;

					if(newGeom.atoms.size() > num-1)
						break;

					iss >> newAtom.number >> temp >> newAtom.x >> newAtom.y >> newAtom.z;

					newAtom.forcex = 0.0;
					newAtom.forcey = 0.0;
					newAtom.forcez = 0.0;

					newGeom.atoms.push_back(newAtom);
				}

				while(line.find("Distance matrix (angstroms):") == std::string::npos){
					if(!std::getline(inputFile, line)){
						ScanFileFormatException excep("Couldn't find 'Distance matrix (angstroms):'. Aborting.");
						throw excep;
					}
				}
				
				while(std::getline(inputFile, line)){
					int num;
					std::string symbol;

					std::istringstream iss(line);

					if(!(iss >> num))
						break;

					iss >> symbol;

					if(symbol[0] >= '0' && symbol[0] <= '9')
						continue;

					newGeom.atoms[num-1].symbol = symbol;
				}
				
				while(line.find("SCF Done") == std::string::npos && std::getline(inputFile, line));

				if(true) {
					std::istringstream iss(line);
					std::string temp;

					iss >> temp >> temp >> temp >> temp >> newGeom.energy;
					newGeom.energy *= 627.509;
				}
				
				while(line.find("Mulliken atomic charges:") == std::string::npos){
					if(!std::getline(inputFile, line)){
						ScanFileFormatException excep("Couldn't find 'Mulliken atomic charges'. Aborting.");
						throw excep;
					}
				}
				//Skip first line after "Mulliken atomic charges"
				std::getline(inputFile, line);

				while(std::getline(inputFile, line)){
					std::istringstream iss(line);
					int num;
					std::string symbol;
					float charge;

					if(!(iss >> num >> symbol >> charge))
						break;

					newGeom.atoms[num-1].charge = charge;


				}
				
				if(hpmodesFound){
					std::cout << "'hpmodes' found. " << std::endl;
					
					std::getline(inputFile, line);
					geoSet.frequencies.clear();
					geoSet.normals.clear();
					
				
					while(line.find("- Thermochemistry -") == std::string::npos){
					
						if(line.find("Force constants ---") != std::string::npos){
							float forceConstant;
							std::string temp;
							std::istringstream iss(line);
						
							// Skip "Force constants ---"
							iss >> temp >> temp >> temp;

							while(iss >> forceConstant){
								geoSet.frequencies.push_back(forceConstant);
								//std::cout << forceConstant << std::endl;
							}
					
						}
						
						if(line.find("Coord Atom Element:") != std::string::npos){
							while(std::getline(inputFile, line)){
								int newDisplacement;

								int axis, atom;
								float displacement;
								std::string temp;

								std::istringstream iss(line);

								iss >> axis >> atom >> temp;

								if(!iss || axis > 3 || axis < 1){ //If all we got was one field, we're done.
									nextDisplacement = newDisplacement;
									break;
								}

								newDisplacement = nextDisplacement;

								while(iss >> displacement){
									while(geoSet.normals.size() <= newDisplacement){
										std::vector<std::vector<float> > newVec;
										geoSet.normals.push_back(newVec);
									}

									while(geoSet.normals[newDisplacement].size() < atom){
										std::vector<float> newVec;
										geoSet.normals[newDisplacement].push_back(newVec);
									}

									geoSet.normals[newDisplacement][atom-1].push_back(displacement);
								
									newDisplacement++;
								}
							}
						}
						if(!std::getline(inputFile, line)){
						ScanFileFormatException excep("Couldn't find '- Thermochemistry -'. Aborting.");
						throw excep;
						}
					}	
					
					
				} else if(!hpmodesFound){
					std::cout << "'hpmodes' not found. " << std::endl;				
				
					std::getline(inputFile, line);
					geoSet.frequencies.clear();
					geoSet.normals.clear();
					
					while(line.find("- Thermochemistry -") == std::string::npos){
					
						if(line.find("Frc consts  --") != std::string::npos){ 
							float forceConstant;
							std::string temp;
							std::istringstream iss(line);

							// Skip "Frc consts  --"
							iss >> temp >> temp >> temp;

							while(iss >> forceConstant){
								geoSet.frequencies.push_back(forceConstant);
								//std::cout << forceConstant << std::endl;
							}
						}
					
						
						
					
						if(line.find("Atom AN") != std::string::npos){
							while(std::getline(inputFile, line)){
								int newDisplacement;
	
								int atom;
								float displacementx, displacementy, displacementz;
								std::string temp;

								std::istringstream iss(line);
								
								iss >> atom >> temp;
								
								if(!iss){ 
									nextDisplacement = newDisplacement+3;
									break;
								}

								newDisplacement = nextDisplacement;

								while(iss >> displacementx >> displacementy >> displacementz){
									while(geoSet.normals.size() <= newDisplacement){
										std::vector<std::vector<float> > newVec;
										geoSet.normals.push_back(newVec);
									}
	
									while(geoSet.normals[newDisplacement].size() < atom){
										std::vector<float> newVec;
										geoSet.normals[newDisplacement].push_back(newVec);
									}

									geoSet.normals[newDisplacement][atom-1].push_back(displacementx);
									geoSet.normals[newDisplacement][atom-1].push_back(displacementy);
									geoSet.normals[newDisplacement][atom-1].push_back(displacementz);
								
									newDisplacement++;
								}
							}
						}
					
						if(!std::getline(inputFile, line)){
							ScanFileFormatException excep("Couldn't find '- Thermochemistry -'. Aborting.");
							throw excep;
						}
					}
				}
					
				//if we want to see displacements, we just need to un-comment section below
				
				//for(std::vector< std::vector< std::vector<float> > >::iterator it3 = geoSet.normals.begin(); it3 != geoSet.normals.end(); it3++){
				//	for(std::vector< std::vector<float> >::iterator it2 = it3->begin(); it2 != it3->end(); it2++){
				//		std::cout << "(" << (*it2)[0] << ", " << (*it2)[1] << ", " << (*it2)[2] << ")\t";
				//	}
				//	std::cout << std::endl;
				//}
				
				
				while(std::getline(inputFile, line)){
					std::string temp;
					int atomNum, number;
					float mass;

					std::istringstream iss(line);

					iss >> temp;

					if(temp == "Molecular")
						break;
					
					if(temp == "Atom"){
						iss >> atomNum >> temp >> temp >> temp >> number >> temp >> temp >> mass;
						newGeom.atoms[atomNum-1].number = number;
						newGeom.atoms[atomNum-1].mass = mass;

						for(std::vector<Geometry>::iterator geomit = geoSet.geometries.begin(); geomit != geoSet.geometries.end(); geomit++){
							geomit->atoms[atomNum-1].number = number;
							geomit->atoms[atomNum-1].mass = mass;
						}
					}
				}
				
				while(line.find("Forces (Hartrees/Bohr)") == std::string::npos){
					if(!std::getline(inputFile, line)){
						ScanFileFormatException excep("Couldn't find 'Forces (Hartrees/Bohr)'. Aborting.");
						throw excep;
					}
				}
				
				//Skip "Number     Number              X              Y              Z"
				std::getline(inputFile, line);
				//Skip "-------------------------------------------------------------------"
				std::getline(inputFile, line);

				while(std::getline(inputFile, line)){
					std::istringstream iss(line);

					int number, atomNum;

					if(!(iss >> number))
						break;
					
					if(!(iss >> atomNum >> newGeom.atoms[number-1].forcex >> 
						newGeom.atoms[number-1].forcey >> newGeom.atoms[number-1].forcez))
						break;
					
					newGeom.atoms[number-1].forcex *= 627.509/0.529177249;
					newGeom.atoms[number-1].forcey *= 627.509/0.529177249;
					newGeom.atoms[number-1].forcez *= 627.509/0.529177249;
				}

				geoSet.geometries.push_back(newGeom);
			
				if(geoSet.normals.size() > 0){
					std::cout << "Please enter a step size: ";
					std::cin >> stepSize;
					std::getline(std::cin, line);
					
					std::vector<bool> usingNormal(geoSet.normals.size(), false);
					std::cout << "List all modes to be included (1 through " << geoSet.normals.size() << ") for geometry " << originalCount << " in " << *it << ", separated by spaces:" << std::endl;
					std::getline(std::cin, line);
					std::istringstream iss2(line);
					int num;
					while(iss2 >> num){
						if(num > 0 && num <= geoSet.normals.size())
							usingNormal[num-1] = true;
					} 
							
					Geometry latest = geoSet.geometries[geoSet.geometries.size()-1];
					int geomNum = 0;
					for(std::vector<std::vector<std::vector<float> > >::iterator modeit = geoSet.normals.begin(); modeit != geoSet.normals.end(); modeit++ , geomNum++){
						if(!usingNormal[geomNum])
							continue;

						Geometry newGeom;

						int atomNum = 0;
						for(std::vector<std::vector<float> >::iterator atomit = modeit->begin(); atomit != modeit->end(); atomit++){
							Atom newAtom;

							newAtom.x = latest.atoms[atomNum].x + (*atomit)[0] * stepSize;
							newAtom.y = latest.atoms[atomNum].y + (*atomit)[1] * stepSize;
							newAtom.z = latest.atoms[atomNum].z + (*atomit)[2] * stepSize;

							if(geomNum < geoSet.frequencies.size()){
								newAtom.forcex = (latest.atoms[atomNum].forcex - stepSize * (*atomit)[0] * geoSet.frequencies[geomNum] * 143.8878);
								newAtom.forcey = (latest.atoms[atomNum].forcey - stepSize * (*atomit)[1] * geoSet.frequencies[geomNum] * 143.8878);
								newAtom.forcez = (latest.atoms[atomNum].forcez - stepSize * (*atomit)[2] * geoSet.frequencies[geomNum] * 143.8878);
							}
							
							//std::cout << newAtom.forcex-latest.atoms[atomNum].forcex << " " << newAtom.forcey-latest.atoms[atomNum].forcey << " " << newAtom.forcez-latest.atoms[atomNum].forcez << std::endl;

							newAtom.number = latest.atoms[atomNum].number;
							newAtom.symbol = latest.atoms[atomNum].symbol;
							newAtom.mass = latest.atoms[atomNum].mass;
							newAtom.charge = latest.atoms[atomNum].charge;
							newGeom.atoms.push_back(newAtom);
							atomNum++;
						}

						geoSet.geometries.push_back(newGeom);
					}
					
					geomNum = 0;
					for(std::vector<std::vector<std::vector<float> > >::iterator modeit = geoSet.normals.begin(); modeit != geoSet.normals.end(); modeit++ , geomNum++){
						if(!usingNormal[geomNum])
							continue;
						
						Geometry newGeom;

						int atomNum = 0;
						for(std::vector<std::vector<float> >::iterator atomit = modeit->begin(); atomit != modeit->end(); atomit++){
							Atom newAtom;

							newAtom.x = latest.atoms[atomNum].x - (*atomit)[0] * stepSize;
							newAtom.y = latest.atoms[atomNum].y - (*atomit)[1] * stepSize;
							newAtom.z = latest.atoms[atomNum].z - (*atomit)[2] * stepSize;

							if(geomNum < geoSet.frequencies.size()){
								newAtom.forcex = (latest.atoms[atomNum].forcex - stepSize * (*atomit)[0] * geoSet.frequencies[geomNum] * 143.8878);
								newAtom.forcey = (latest.atoms[atomNum].forcey - stepSize * (*atomit)[1] * geoSet.frequencies[geomNum] * 143.8878);
								newAtom.forcez = (latest.atoms[atomNum].forcez - stepSize * (*atomit)[2] * geoSet.frequencies[geomNum] * 143.8878);
							}

							newAtom.number = latest.atoms[atomNum].number;
							newAtom.symbol = latest.atoms[atomNum].symbol;
							newAtom.mass = latest.atoms[atomNum].mass;
							newAtom.charge = latest.atoms[atomNum].charge;
							newGeom.atoms.push_back(newAtom);
							atomNum++;
						}

						geoSet.geometries.push_back(newGeom);
					}
				}
			}			
		}
	}


	return geoSet;
}
