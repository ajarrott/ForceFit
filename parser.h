/**********************
 *     Created by     *
 *   Anthony Arrott   *
 *      3/14/2015     *
 **********************/

#ifndef _parser_h_
#define _parser_h_

#define GAUSSIAN 1
#define NWCHEM 2
#define CRYSTAL 3

#define DL_POLY 4
#define LAMMPS 5
#define TINKER 6

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

class parser
{
	public:
		/* constructors */
		parser();
		~parser();

		void createNewTextFile();
		void updateFilename(string);
		bool run();
		int getQType();
		int getQMSW();
		vector<string> getQMFiles();
		vector<double> getGuess();
		vector<string> getMDFiles();
		bool getForces();
		bool getEnergies();
		string getStepSize();

	private:
		/* member variables */
		ifstream inputFile;			// initial input file
		ifstream temp;
		string filename;
		string line;
		unsigned int i, QType, lineCount;

		vector<string> qmFiles;		// quantum mechanic files (geometry sets)
		int qmsw;					// gaussian, NWChem, or crystal
		vector<double> guess;		// guess parameters, e.g. a = 2.1032....
									// later in the program the value names are not used
									// therefore I just need a vector of the values
		vector<string> md_files;	// md files
		bool forces;				// use forces
		bool energies;				// use energies
		string step_size;

		/* functions */
		bool openFile();			// open input file designated/default
		bool checkFile(string);		// check to see if file exists
		void findNextLine();		// helper function to get to next non blank line
		bool checkQ(string);		// better user experence
		bool mdHelper(string, string);
		bool skipToQuestion(string);
		bool lastQuestions(string, bool);
		bool setQMFiles();			// sets the geometries for the specific files
		bool setQMSW();				// sets quantum mechanics software to gaussian/NWChem
		bool setGuessParam();		// set Molecular Dynamics software to dl_poly/lammps/tinker
		bool setMDSim();			// set md_file
		bool finalQuestions();

		
		
		
};

#endif
