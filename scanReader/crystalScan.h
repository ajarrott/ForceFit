#ifndef _CRYSTALSCAN_H
#define _CRYSTALSCAN_H

#include <string>

#include "../geometrySet.h"
#include "scanReader.h"

class CrystalScan : public ScanReader {
public:
	virtual GeometrySet getGeometries();
	virtual std::string const getName();
	static std::string const name;
};

#endif
