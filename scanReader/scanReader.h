#ifndef _SCANREADER_H
#define _SCANREADER_H

#include <string>
#include <vector>
#include <exception>

#include "../geometrySet.h"

class ScanReader {
public:
	~ScanReader();

	virtual GeometrySet getGeometries() = 0;
	
	std::vector<std::string> inputFiles;

	static std::string const name;
	virtual std::string const getName() = 0;
};

class ScanFileFormatException : public std::exception {
public:
	ScanFileFormatException(const char* message);

	virtual const char * what() const throw();
private:
	const char* whatMessage;
};

#endif
