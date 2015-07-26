#ifndef _nogui_h_
#define _nogui_h_

#include "parser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>

//#include "mdQuestionWindow.h"
//#include "minimizerQuestionWindow.h"
#include "geometrySet.h"
#include "geometry.h"

#include "pugixml/pugixml.hpp"

#define ScanReaderInclude
#define MolecularDynamicsInclude
#define GradientCreatorInclude
#define MinimizerInclude
#include "classes.h"
#undef MolecularDynamicsInclude
#undef GradientCreatorInclude
#undef MinimizerInclude
#undef ScanReaderInclude

class nogui
{
	public:
		/* constructors */
		nogui(string);
		~nogui();

		void run();

	private:
		parser p;
		vector<string> qmFilenames;
		vector<string> mdFilenames;
		vector<GeometrySet> geoSets;
		GeometrySet set;
		//vector<GeometrySet *> setPtrs;
		vector<double> variables;
		int qmsw, QType;
		vector<struct mdQuestion> mdQuestions;
		vector<struct minQuestion> minQuestions;
		ScanReader * reader;
		MolecularDynamics * md;
		Minimizer * min;

		void fillMinAnswers();
		void fillMDAnswers();
		void setGeometrySets();
};

#endif
