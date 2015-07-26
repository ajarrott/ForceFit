#include "nwchemScan.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <utility>

#define PI 3.14159265
#define c 29979245800 //in cm/s

std::string const NWChemScan::name = "NWChem Output Type";

std::string const NWChemScan::getName(){
	return name;
}

GeometrySet NWChemScan::getGeometries(){
	GeometrySet geoSet;

	double stepSize;
	int originalCount = 0;

	geoSet.name = "NWChem Set from ";

	bool comma = false;
        for(std::vector<std::string>::iterator it = inputFiles.begin(); it != inputFiles.end(); it++){
                if(comma)
                        geoSet.name += ", ";
                comma = true;
                geoSet.name += *it;
	
		std::ifstream inputFile(it->c_str());

		bool freqFound = false;
		bool gradientFound = false;

		std::string line;

		while(std::getline(inputFile, line)){

			if(line.find("task dft freq") != std::string::npos){
				freqFound = true;
				std::cout << "Frequency calculations" << std::endl;
			} else if(line.find("task dft gradient") != std::string::npos) {
				gradientFound = true;
				std::cout << "Gradient calculations" << std::endl;
			} else if(line.find("task mp2 gradient") != std::string::npos) {
				gradientFound = true;
				std::cout << "Gradient calculations" << std::endl;
			}


			if(line.find("Geometry \"complete\"") != std::string::npos || line.find("Geometry \"geometry\"") != std::string::npos){
				originalCount++;
	
				while(std::getline(inputFile, line) && line.find("  No.       Tag          Charge          X              Y              Z") == std::string::npos);
	
				if(line.find("  No.       Tag          Charge          X              Y              Z") == std::string::npos){
					ScanFileFormatException excep("Couldn't find geometries. Aborting.");
					throw excep;
				}
				
				std::getline(inputFile, line);    //Skip " ---- ---------------- ---------- -------------- -------------- --------------"
	
				Geometry newGeom;
				int atomicNumber = 0;
				Atom newAtom;
				
				if(freqFound){
					while(std::getline(inputFile, line)){
					
						int num;
					
						std::string symbol;
					
						std::istringstream iss(line);

						if(!(iss >> num >> symbol >> atomicNumber >> newAtom.x >> newAtom.y >> newAtom.z))
							break;

						if(newGeom.atoms.size() > num-1) //Make sure we don't skip a number
							break;

						newAtom.number = (int)atomicNumber;
						newAtom.symbol = symbol;
					
						newAtom.forcex = 0.0;
						newAtom.forcey = 0.0;
						newAtom.forcey = 0.0;

						for(std::vector<Geometry>::iterator geomit = geoSet.geometries.begin(); geomit != geoSet.geometries.end(); geomit++){
							geomit->atoms[num-1].number = (int)atomicNumber;
						}
					
						newGeom.atoms.push_back(newAtom);
					}
				
					std::getline(inputFile, line);
				
					
					geoSet.geometries.push_back(newGeom);

		
					while(line.find("NWChem DFT Module") == std::string::npos){
						if(!std::getline(inputFile, line)){
							ScanFileFormatException excep("Couldn't find 'NWChem DFT Module'. Aborting.");
							throw excep;
						} 
					}
				
					std::getline(inputFile, line);
					
					while(line.find("Total Density - Mulliken Population Analysis") == std::string::npos){
						if(!std::getline(inputFile, line)){
							ScanFileFormatException excep("Couldn't find 'Total Density - Mulliken Population Analysis'.");
							throw excep;
						}
					}
				
					std::getline(inputFile, line); //We're making sure, that we extract second set of charges
				
					while(line.find("Total Density - Mulliken Population Analysis") == std::string::npos){
						if(!std::getline(inputFile, line)){
							ScanFileFormatException excep("Couldn't find 'Total Density - Mulliken Population Analysis'.");
							throw excep;
						}
					}
				
				
					std::getline(inputFile, line);    //Skip --------------------------------------------- 
					std::getline(inputFile, line);    //Skip blank line
					std::getline(inputFile, line);    //Skip "Atom       Charge   Shell Charges"
					std::getline(inputFile, line);    //Skip -----------   ------   -------------------------------------------------------
					
					while(std::getline(inputFile, line)){
						std::istringstream iss(line);
						int num;
						std::string symbol;
						float atomNum;
						float totalCharge;
	
						if(!(iss >> num >> symbol >> atomNum >> totalCharge))
							break;
		
						if(num-1 < newGeom.atoms.size())
							newGeom.atoms[num-1].charge = atomNum-totalCharge;
						
						for(std::vector<Geometry>::iterator geomit = geoSet.geometries.begin(); geomit != geoSet.geometries.end(); geomit++){
							geomit->atoms[num-1].charge = atomNum-totalCharge;	
						}
					
						//std::cout << atomNum-totalCharge << std::endl;
						
					}
					
					while(line.find("finite difference hessian delta") == std::string::npos){
				
						Geometry newGeom;
						Atom newAtom;	
					
						if(line.find("NWChem DFT Module") != std::string::npos){
						
						
							std::getline(inputFile, line);    //Skip -----------------
							std::getline(inputFile, line);    //Skip blank line
							std::getline(inputFile, line);    //Skip blank line
							std::getline(inputFile, line);    //Skip line with title
							std::getline(inputFile, line);    //Skip blank line
							std::getline(inputFile, line);    //Skip blank line
							std::getline(inputFile, line);    //Skip "Caching 1-el integrals"
							std::getline(inputFile, line);    //Skip "Time after variat. SCF:"
							std::getline(inputFile, line);    //Skip blank line
							std::getline(inputFile, line);    //Skip ---------------------------------------------
							std::getline(inputFile, line);    //Skip "Total Density - Mulliken Population Analysis"
							std::getline(inputFile, line);    //Skip ---------------------------------------------
							std::getline(inputFile, line);    //Skip blank line
							std::getline(inputFile, line);    //Skip "Atom       Charge   Shell Charges"
							std::getline(inputFile, line);    //Skip -----------   ------   -------------------------------------------------------
						
							while(std::getline(inputFile, line)){
								std::istringstream iss(line);
								int num;
								std::string symbol;
								float atomNum;
								float totalCharge;
	
								if(!(iss >> num >> symbol >> atomNum >> totalCharge))
									break;
		
								if(num-1 < newGeom.atoms.size())
									newGeom.atoms[num-1].charge = atomNum-totalCharge;
							
								newAtom.charge = atomNum-totalCharge;
						
								for(std::vector<Geometry>::iterator geomit = geoSet.geometries.begin(); geomit != geoSet.geometries.end(); geomit++){
									geomit->atoms[num-1].charge = atomNum-totalCharge;
								}
							
								newGeom.atoms.push_back(newAtom);
						
							}
						
						}
					
						if(line.find("DFT ENERGY GRADIENTS") != std::string::npos){
					
						
							std::getline(inputFile, line);    //Skip blank line
							std::getline(inputFile, line);    //Skip "atom               coordinates                        gradient"
							std::getline(inputFile, line);    //Skip "x          y          z           x          y          z"

							while(std::getline(inputFile, line)){
								std::istringstream iss(line);
								std::string temp;
								std::string symbol;
								float x, y, z, forcex, forcey, forcez;
								
								iss >> temp >> symbol >> x >> y >> z >> forcex >> forcey >> forcez;
								if(!(iss))
									break;
							
								newAtom.symbol = symbol;	
								newAtom.x = x*0.529177249;
								newAtom.y = y*0.529177249;
								newAtom.z = z*0.529177249;
								newAtom.forcex = forcex*627.509/0.529177249;
								newAtom.forcey = forcey*627.509/0.529177249;
								newAtom.forcez = forcez*627.509/0.529177249;
							

								newGeom.atoms.push_back(newAtom);
									
							}
					
						geoSet.geometries.push_back(newGeom);
					
						}
					
						
						if(!std::getline(inputFile, line))
							break;
					}

					while(line.find("Atom information") == std::string::npos){
						if(!std::getline(inputFile, line)){
							ScanFileFormatException excep("Couldn't find 'Atom information'.Aborting.");
							throw excep;
						}
					}
				
					std::getline(inputFile, line);    //Skip "atom    #        X              Y              Z            mass"
					std::getline(inputFile, line);    //Skip "--------------------------------------------------------------------------"

					while(std::getline(inputFile, line)){
						float mass;
						std::string massString;
						std::istringstream iss(line);
						std::string symbol;
						std::string temp;
					
						int num;
			
						if(!(iss >> symbol >> num >> temp >> temp >> temp >> massString))
							break;


						if(massString.find("D") != std::string::npos)
							massString[massString.find("D")] = 'E';
						
						std::istringstream massIss(massString);
						massIss >> mass;
	
						newGeom.atoms[num-1].symbol = symbol;
						newGeom.atoms[num-1].mass = mass;
					
						for(std::vector<Geometry>::iterator geomit = geoSet.geometries.begin(); geomit != geoSet.geometries.end(); geomit++){
							geomit->atoms[num-1].mass = mass;
						}
						
						//std::cout << symbol << " " << mass << std::endl;
						
						newGeom.atoms.push_back(newAtom);

					}
					
				} else if(gradientFound){

					using namespace std;
					typedef map<string, float> map1;
					map1 atomicNumber;
					std::string sym;
					float atomNumber;

					while(std::getline(inputFile, line)){
						
						std::istringstream iss(line);
						std::string temp;
						float atomNumber;

						iss >> temp >> sym >> atomNumber >> temp >> temp >> temp;

						atomicNumber.insert(pair<string, float>(sym,atomNumber));

						if(!iss)
							break;

						//std::cout << atomicNumber[sym] << std::endl;
					}
					
					
					std::getline(inputFile, line);
					
					while(line.find("Atomic Mass") == std::string::npos){
						if(!std::getline(inputFile, line)){
							ScanFileFormatException excep("Couldn't find 'Atomic Mass'.");
							throw excep;
						}
					}
					
					std::getline(inputFile, line); // Skip -----------
					std::getline(inputFile, line); // Skip blank line 
					
					using namespace std;
					typedef map<string, float> map2;
					map2 atomicMass;
					//map2::iterator aM = atomicMass.begin();
					std::string symbol;
					float mass;
					
					while(std::getline(inputFile, line)){

						std::istringstream iss(line);
						
						iss >> symbol >> mass;
						atomicMass.insert(pair<string, float>(symbol,mass));
						
						if(!iss)
							break;
					
						//std::cout << atomicMass[symbol] << std::endl;
					}

// Commented out Dan 5/16/12
//                    
//					while(line.find("NWChem DFT Module") == std::string::npos){
//						if(!std::getline(inputFile, line)){
//							ScanFileFormatException excep("Couldn't find 'NWChem DFT Module'. Aborting.");
//							throw excep;
//						} 
//					}
//				
//					std::getline(inputFile, line);
//					
//					while(line.find("DFT Final Molecular Orbital Analysis") == std::string::npos && line.find("DFT Final Beta Molecular Orbital Analysis") == std::string::npos){
//						if(!std::getline(inputFile, line)){
//							ScanFileFormatException excep("Couldn't find 'DFT Final Molecular Orbital Analysis' nor 'DFT Final Beta Molecular Orbital Analysis'.");
//							throw excep;
//						}
//					}
//					
//					std::getline(inputFile, line); //We're making sure, that we extract second set of charges
//				
//					while(line.find("Total Density - Mulliken Population Analysis") == std::string::npos){
//						if(!std::getline(inputFile, line)){
//							ScanFileFormatException excep("Couldn't find 'Total Density - Mulliken Population Analysis'.");
//							throw excep;
//						}
//					}
//					
//					std::getline(inputFile, line);    //Skip --------------------------------------------- 
//					std::getline(inputFile, line);    //Skip blank line
//					std::getline(inputFile, line);    //Skip "Atom       Charge   Shell Charges"
//					std::getline(inputFile, line);    //Skip -----------   ------   -------------------------------------------------------
//					
//					using namespace std;
//					typedef map<int, float> map3;
//					map3 charges;
//					
//					float atomNum;
//					float totalCharge;
//					int num;
//					
//					while(std::getline(inputFile, line)){
//						std::istringstream iss(line);
//						
//						std::string symbol;
//						
//						iss >> num >> symbol >> atomNum >> totalCharge;
//						
//						charges.insert(pair<int, float>(num,(atomNum-totalCharge)));
//						
//						if(!iss)
//							break;
//						
//						//std::cout << charges[num] << std::endl;
//						//if(num-1 < newGeom.atoms.size())
//						//	newGeom.atoms[num-1].charge = atomNum-totalCharge;
//						
//						//for(std::vector<Geometry>::iterator geomit = geoSet.geometries.begin(); geomit != geoSet.geometries.end(); geomit++){
//						//	geomit->atoms[num-1].charge = atomNum-totalCharge;	
//						//}
//					
//						//std::cout << atomNum-totalCharge << std::endl;
//						
//						//newGeom.atoms.push_back(newAtom);
//						
//					}

					while(line.find("ACKNOWLEDGEMENT") == std::string::npos){
				
						Geometry newGeom;
						Atom newAtom;	
					
						if(line.find("ENERGY GRADIENTS") != std::string::npos){
						
							std::getline(inputFile, line);    //Skip blank line
							std::getline(inputFile, line);    //Skip "atom               coordinates                        gradient"
							std::getline(inputFile, line);    //Skip "x          y          z           x          y          z"

							while(std::getline(inputFile, line)){
								std::istringstream iss(line);
								std::string temp;
								std::string symbol;
								int num;
								float x, y, z, forcex, forcey, forcez, mass;
								
								iss >> num >> symbol >> x >> y >> z >> forcex >> forcey >> forcez;
								if(!(iss))
									break;

								newAtom.symbol = symbol;	
								newAtom.x = x*0.529177249;
								newAtom.y = y*0.529177249;
								newAtom.z = z*0.529177249;
								newAtom.forcex = forcex*627.509/0.529177249;
								newAtom.forcey = forcey*627.509/0.529177249;
								newAtom.forcez = forcez*627.509/0.529177249;
								
								for(map1::const_iterator aN = atomicNumber.begin(); aN != atomicNumber.end(); ++aN){
									if(newAtom.symbol == (*aN).first){
										newAtom.number = (*aN).second;
									}
								}

								//std::cout << newAtom.number << std::endl;

								for(map2::const_iterator aM = atomicMass.begin(); aM != atomicMass.end(); ++aM){
									if(newAtom.symbol == (*aM).first){
										newAtom.mass = (*aM).second;
									}
								}
								
//								for(map3::const_iterator ch = charges.begin(); ch != charges.end(); ++ch){
//									if(num == (*ch).first){
//										newAtom.charge = (*ch).second;
//									}
//								}
							
								newGeom.atoms.push_back(newAtom);
									
							}
					
						geoSet.geometries.push_back(newGeom);
					
						}
					
						
						if(!std::getline(inputFile, line))
							break;
					}
				}	
					
			}
		}
	}
	return geoSet;

}
