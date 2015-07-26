#ifndef _NWCHEMSCAN_H
#define _NWCHEMSCAN_H

#include <string>

#include "../geometrySet.h"
#include "scanReader.h"

class NWChemScan : public ScanReader {
public:
	virtual GeometrySet getGeometries();
	virtual std::string const getName();
	static std::string const name;
};

#endif
