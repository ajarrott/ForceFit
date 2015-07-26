#include "crystalGradient.h"

#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

#define PI 3.14159265
	
std::string const CrystalGradient::name = "Crystal Gradient Creator";

std::string const CrystalGradient::getName(){
	return name;
}

void writeGradient(std::string filename, std::string beforeText, std::string afterText, std::string str){
	std::ofstream fp(filename.c_str());

	fp << beforeText << str << afterText;

	fp.close();
}

void CrystalGradient::createGradients(GeometrySet geoSet, std::string templateFile, float distance, int steps, int startMode, int endMode){
	std::string line;
	std::ostringstream bfOss, afOss;

	std::ifstream templ(templateFile.c_str());
	
	std::getline(templ, line); //Title text
	bfOss << line << std::endl;
	
	std::getline(templ, line);
	bfOss << line << std::endl;
	if(line.find("CRYSTAL") != std::string::npos){
		std::getline(templ, line);
		bfOss << line << std::endl;
	}
	std::getline(templ, line);
	bfOss << line << std::endl;

	bfOss << geoSet.a << ", " << geoSet.b << ", " << geoSet.c << ", ";
	bfOss << geoSet.alpha << ", " << geoSet.beta << ", " << geoSet.gamma << std::endl;

	// bfOss << geoSet.originalGeometry.atoms.size() << std::endl; FIXME

	while(line.find("OPTGEOM") == std::string::npos && !templ.eof())
		std::getline(templ, line);

	while(!templ.eof()){
		afOss << line << std::endl;
		std::getline(templ, line);
	}
	afOss << line << std::endl;

	templ.close();
	
	if(startMode < 1)
		startMode = 1;
	
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

	for(int currentStep = 1; currentStep <= steps; currentStep++){
		for(int currentMode = startMode; currentMode <= endMode && currentMode <= geoSet.normals.size(); currentMode++){
			std::ostringstream ossNeg;
			std::ostringstream ossPos;

			for(int currentAtom = 0; currentAtom < geoSet.normals[currentMode-1].size(); currentAtom++){
				float x,y,z;
				float posX, posY, posZ, negX, negY, negZ;
				float fracPosX, fracPosY, fracPosZ, fracNegX, fracNegY, fracNegZ;

				/* FIXME
				x = geoSet.originalGeometry.atoms[currentAtom].x;
				y = geoSet.originalGeometry.atoms[currentAtom].y;
				z = geoSet.originalGeometry.atoms[currentAtom].z;
				*/

				posX = x + currentStep * (distance/0.001) * geoSet.normals[currentMode-1][currentAtom][0] * 0.5291772;
				posY = y + currentStep * (distance/0.001) * geoSet.normals[currentMode-1][currentAtom][1] * 0.5291772;
				posZ = z + currentStep * (distance/0.001) * geoSet.normals[currentMode-1][currentAtom][2] * 0.5291772;
				
				negX = x - currentStep * (distance/0.001) * geoSet.normals[currentMode-1][currentAtom][0] * 0.5291772;
				negY = y - currentStep * (distance/0.001) * geoSet.normals[currentMode-1][currentAtom][1] * 0.5291772;
				negZ = z - currentStep * (distance/0.001) * geoSet.normals[currentMode-1][currentAtom][2] * 0.5291772;

				fracPosX = posX * inverseTransformationMatrix[0][0] + 
					posY * inverseTransformationMatrix[0][1] +
					posZ * inverseTransformationMatrix[0][2];
				fracPosY = posX * inverseTransformationMatrix[1][0] + 
					posY * inverseTransformationMatrix[1][1] +
					posZ * inverseTransformationMatrix[1][2];
				fracPosZ = posX * inverseTransformationMatrix[2][0] + 
					posY * inverseTransformationMatrix[2][1] +
					posZ * inverseTransformationMatrix[2][2];

				fracNegX = negX * inverseTransformationMatrix[0][0] + 
					negY * inverseTransformationMatrix[0][1] +
					negZ * inverseTransformationMatrix[0][2];
				fracNegY = negX * inverseTransformationMatrix[1][0] + 
					negY * inverseTransformationMatrix[1][1] +
					negZ * inverseTransformationMatrix[1][2];
				fracNegZ = negX * inverseTransformationMatrix[2][0] + 
					negY * inverseTransformationMatrix[2][1] +
					negZ * inverseTransformationMatrix[2][2];

				ossPos.width(3);
				// ossPos << geoSet.originalGeometry.atoms[currentAtom].number << "   "; FIXME
				ossPos.width(19);
				ossPos << fracPosX << " ";
				ossPos.width(19);
				ossPos << fracPosY << " ";
				ossPos.width(19);
				ossPos << fracPosZ << std::endl;

				ossNeg.width(3);
				// ossNeg << geoSet.originalGeometry.atoms[currentAtom].number << "   "; FIXME
				ossNeg.width(19);
				ossNeg << fracNegX << " ";
				ossNeg.width(19);
				ossNeg << fracNegY << " ";
				ossNeg.width(19);
				ossNeg << fracNegZ << std::endl;
			}
		
			std::ostringstream nameOssPos(std::ostringstream::out);
			nameOssPos << "sio2gradm";
			nameOssPos << currentMode;
			nameOssPos << "s";
			nameOssPos << currentStep;
			nameOssPos << "p.d12";
			writeGradient(nameOssPos.str(), bfOss.str(), afOss.str(), ossPos.str());
			
			std::ostringstream nameOssNeg(std::ostringstream::out);
			nameOssNeg << "sio2gradm";
			nameOssNeg << currentMode;
			nameOssNeg << "s";
			nameOssNeg << currentStep;
			nameOssNeg << "n.d12";
			writeGradient(nameOssNeg.str(), bfOss.str(), afOss.str(), ossNeg.str());
		}
	}
}

GeometrySet CrystalGradient::readGradients() {
	GeometrySet geoSet;

	return geoSet;
}

