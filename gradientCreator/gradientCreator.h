#ifndef _GRADIENTCREATOR_H
#define _GRADIENTCREATOR_H

#include "../geometrySet.h"

class GradientCreator {
public:
	virtual void createGradients(GeometrySet geoSet, std::string templateFile, float distance, int steps, int startMode, int endMode) = 0;
	virtual GeometrySet readGradients() = 0;
	static std::string const name;
	virtual std::string const getName() = 0;
};

#endif
