#ifndef _GAUSSIANSCAN_H
#define _GAUSSIANSCAN_H

#include <string>

#include "../geometrySet.h"
#include "scanReader.h"

class GaussianScan : public ScanReader {
public:
	virtual GeometrySet getGeometries();
	virtual std::string const getName();
	static std::string const name;
};

#endif
