#ifndef _CRYSTALGRADIENT_H
#define _CRYSTALGRADIENT_H

#include "gradientCreator.h"

class CrystalGradient : public GradientCreator {
public:
	void createGradients(GeometrySet geoSet, std::string templateFile, float distance, int steps, int startMode, int endMode);
	GeometrySet readGradients();
	static std::string const name;
	virtual std::string const getName();
};

#endif
