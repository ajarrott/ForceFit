#include "scanReader.h"

std::string const name = "Generic Scan Reader";

ScanReader::~ScanReader(){
}

ScanFileFormatException::ScanFileFormatException(const char * message){
	whatMessage = message;
}

const char * ScanFileFormatException::what() const throw(){
	return whatMessage;
}
