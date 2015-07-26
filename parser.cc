#include "parser.h"

// PRIVATE FUNCTIONS
bool parser::checkFile(string name)
{
	struct stat st;

	stat(name.c_str(), &st);

	if (st.st_mode & S_IFDIR )
	{
		cout << name << " is a directory, not a file." << endl
			<< "Please add the filename to the end of the path" << endl
			<< "Exiting ForceFit." << endl;
		exit(-1);
	}
	cout << "Checking " << name;
	ifstream f(name.c_str());

	if (f.good()){
		f.close();
		cout << "...OK!" << endl;
		return true;
	}
	else
	{
		f.close();
		cout << "...Failed!" << endl;
		return false;
	}
}

// helper function to reduce code
bool parser::checkQ(string question)
{
	int qSize = question.length();

	if ( !line.compare(0,qSize,question) )
		return true;
	return false;
	
}

void parser::findNextLine()
{
	getline(inputFile, line);
	lineCount++;

	// for comments
	while ( line[0] == '#' )
	{
		getline(inputFile, line);
		lineCount++;
	}
}

bool parser::setQMFiles()
{
	// START QM Files
	string listOfQM = "List of QM files:";
	findNextLine();

	while (line.empty())
	{
		findNextLine();
	}

	// get first part of files, if line does not compare, will not return 0, will fail

	if (!checkQ(listOfQM))
	{
		cout << "Invalid input file, line " << lineCount << " should be \"" << listOfQM << "\"" << endl
			<< "without quotes." << endl;
		return false;
	}


	findNextLine();

	while (line.length())
	{
		if ( checkFile(line) )
		{
			qmFiles.push_back(line);
			findNextLine();
		}
		else
		{
			cout << "Error, file from path: " << endl
			<< line << endl
			<< "at line " << lineCount << " does not exist." << endl;
			return false;
		}
	}

	if (!qmFiles.size())
	{
		cout << "No QM files given at line " << lineCount 
		<< ", please add some to your input file." 
		<< endl;
		return false;
	}

	cout << "QM Files read successfully." << endl;
	// END QM Files
	return true;
}

bool parser::setQMSW()
{
	// START QM Software
	string QMSW = "QM software used to create these files:";

	// get next line, after blank space(s)
	findNextLine();
	

	while (line.empty())
	{
		findNextLine();
	}

	if (!checkQ(QMSW))
	{
		cout << "Invalid input file, line " << lineCount << " should be \"" << QMSW << "\"" << endl
			<< "without quotes." << endl;
		return false;
	}

	findNextLine();

	// make line tolower for different user inputs choices
	for (i = 0; i < line.length(); i++)
	{
		line[i] = tolower(line[i]);
	}

	if (checkQ("guassian"))
	{
		qmsw = GAUSSIAN;
	}
	else if (checkQ("nwchem"))
	{
		qmsw = NWCHEM;
	}
	else if (checkQ("crystal"))
	{
		qmsw = CRYSTAL;
	}
	else
	{
		cout << "Invalid software type, line should be \"NWChem\", \"Gaussian\", or \"Crystal\"" << endl
			<< "without quotes at line " << lineCount << "." << endl;
		return false;
	}

	cout << "QM Software loaded successfully." << endl;
	// END QM Software
	return true;
}

bool parser::setGuessParam()
{
	// START Guess Parameters
	string guessParam = "Guess parameters for force field:";
	findNextLine();

	while (line.empty())
	{
		findNextLine();
	}

	if (!checkQ(guessParam))
	{
		cout << "Invalid input file, line" << lineCount << " should be \"" << guessParam << "\"" << endl
			<< "without quotes at line " << lineCount << "." << endl;
		return false;
	}

	findNextLine();


	while (!line.empty())
	{
		char* dup = strdup(line.c_str());
		char* token = NULL;
		string newStr;
		double newVal;
		char* context = NULL;
		char delimiters[] = "\t= ";
		// set a max loops in case the user tries to break strtok
		int maxLoops = 20;

		token = strtok_r(dup, delimiters, &context);
		newStr = token;

		// get guess value
		try
		{
			token = strtok_r(NULL, delimiters, &context);

			// don't parse possible null tokens 
			while ( !token && maxLoops)
			{
				token = strtok_r(NULL, delimiters, &context);
				maxLoops--;
			}

			if ( !maxLoops )
			{
				cout << "ERROR parsing guess variables input file at line " << lineCount
					<< "." << endl
					<< "Make sure your assignment is \"(variable) = (value)\"" << endl;
				exit(-1);
			}
			newVal = atof(token);
		}
		catch (...)
		{
			cout << "An error occured while parsing line " << lineCount << ":" << endl
				<< line << endl
				<< "Make sure the formatting is (Variable)(Space or Tab)(Number)"
				<< "without parenthesis." << endl;
			return false;
		}

		// add items to vector
		guess.push_back(newVal);

		cout << "Added guess value " << newVal << "." << endl;

		findNextLine();
	}

	if (!guess.size())
	{
		cout << "No guess variables declared on line " << lineCount << endl
			<< "Please add some." << endl;

		return false;
	}

	cout << "Guess variables successfully added." << endl;

	// END Guess Parameters
	return true;
}

bool parser::setMDSim()
{
	// START MD Simulation
	string MDSim = "MD simulations software to be used during the fit:";
	string dl_poly = "dl_poly";
	string lammps = "lammps";
	string tinker = "tinker";

	findNextLine();

	// get to next non-empty line
	while (line.empty())
	{
		findNextLine();
	}

	if (!checkQ(MDSim))
	{
		cout << "Invalid input file, line" << lineCount << " should be \"" << MDSim << "\"" << endl
			<< "without quotes at line " << lineCount << "." << endl;
		return false;
	}

	findNextLine();

	// make line tolower for different user inputs choices
	for (i = 0; i < line.length(); i++)
	{
		line[i] = tolower(line[i]);
	}

	if (checkQ(dl_poly))
	{
		cout << "Setting software to DL_Poly." << endl;
		QType = DL_POLY;
	}
	else if (checkQ(lammps))
	{
		cout << "Setting software to lammps" << endl;
		QType = LAMMPS;
	}
	else if (checkQ(tinker))
	{
		cout << "Setting software to tinker." << endl;
		QType = TINKER;
	}
	else
	{
		cout << "Invalid input file, line " << lineCount << " should be: " << endl
			<< "\"" << dl_poly << "\", " << "\"" << lammps << "\", or " << "\"" << tinker << "\"" << endl
			<< "without quotes at line " << lineCount << "." << endl;

		return false;
	}

	// END MD Simulation

	cout << "MD Simulation finished successfully." << endl;
	return true;
}

bool parser::mdHelper(string question1, string question2)
{
	// fix bug where skipping questions goes one past where it should be
	if ( !checkQ(question1) )
	{
		findNextLine();
		while(line.empty()) { findNextLine(); }
		if (!checkQ(question1))
		{
			cout << "Invalid input file" << endl
			<< "Line " << lineCount << " should be \"" << question1 << "\"" << endl;
			return false;
		}
	}

	// next line should be the inputfile
	findNextLine();

	while(line.empty())
	{
		findNextLine();
	}

	if (checkQ(question2))
	{
		cout << "No input file given for Input at line" << lineCount - 2 << endl;
		return false;
	}

	if ( checkFile(line) )
	{
		md_files.push_back(line);
		findNextLine();
	}
	else
	{
		cout << "Error, file from path: " << endl
		<< line << endl
		<< "at line " << lineCount << " does not exist." << endl;
		return false;
	}

	
	return true;
}

bool parser::skipToQuestion(string question)
{
	// get to next question set (LAMMPS)
	while (!inputFile.eof())
	{
		findNextLine();

		if (checkQ(question))
			break;
	}

	if (inputFile.eof())
	{
		cout << "Unexpected end of input file at line " << lineCount << "." << endl
		<< "Leave all questions in the input file with blank entries. " << endl;
		return false;
	}

	return true;
}

bool parser::lastQuestions(string question, bool type)
{
	findNextLine();

	while (line.empty())
	{
		findNextLine();
	}

	if (checkQ(question))
	{
		findNextLine();

		while (line.empty())
		{
			findNextLine();
		}

		for (i = 0; i < line.length(); i++)
		{
			line[i] = tolower(line[i]);
		}

		if (checkQ("y") || checkQ("yes"))
		{
			type = true;
		}
		else if (checkQ("n") || checkQ("no"))
		{
			type = false;
		}
		else
		{
			cout << "Invalid input file, line" << lineCount << " should be \"Y\" or \"N\"" << endl
				<< "without quotes at line " << lineCount << "." << endl;
			return false;
		}
	}

	return true;
}

bool parser::finalQuestions()
{
	// START Final Questions
	string dlQ1 = "For Dl_Poly provide path to Control:";
	string dlQ2 = "For Dl_Poly provide path to Field:";
	string dlQ3 = "For Dl_poly provide path to Dl_Poly executable:";

	string lammpsQ1 = "For Lammps provide path to Input:";
	string lammpsQ2 = "For Lammps provide path to Data:";
	string lammpsQ3 = "For Lammps provide path to Lammps executable:";

	string tinkerQ1 = "For Tinker provide path to XYZ Template:";
	string tinkerQ2 = "For Tinker provide path to PRM:";
	string tinkerQ3 = "For Tinker provide path to Key:";
	string tinkerQ4 = "For Tinker provide path to Tinker executable:";

	string forceQ = "Are forces to be used in fitting:";
	string energyQ = "Are energies to be used in fitting:";
	string stepQ = "Step size:";

/////////////// DL_POLY
	if (QType == DL_POLY )
	{
		// control
		if (!mdHelper(dlQ1, dlQ2))
			return false;
		// field
		if (!mdHelper(dlQ2, dlQ3))
			return false;
		// executable
		if (!mdHelper(dlQ3, lammpsQ1))
			return false;

		cout << "Read DL_Poly information successfully." << endl;
	}
	else
	{
		if(!skipToQuestion(lammpsQ1))
			return false;
	}

/////////////// LAMMPS
	if (QType == LAMMPS)
	{
		// input
		if (!mdHelper(lammpsQ1, lammpsQ2))
			return false;
		// data
		if (!mdHelper(lammpsQ2, lammpsQ3))
			return false;
		// executable
		if (!mdHelper(lammpsQ3, tinkerQ1))
			return false;

		cout << "Read LAMMPS information successfully." << endl;
	}
	else
	{
		if(!skipToQuestion(tinkerQ1))
			return false;
	}
///////////////// TINKER
	if (QType == TINKER)
	{
		// XYZ
		if (!mdHelper(tinkerQ1, tinkerQ2))
			return false;
		// PRM
		if (!mdHelper(tinkerQ2, tinkerQ3))
			return false;
		// Key
		if (!mdHelper(tinkerQ3, tinkerQ4))
			return false;
		// executable
		if (!mdHelper(tinkerQ4, forceQ))
			return false;

		cout << "Read Tinker information successfully." << endl;
	}
	else
	{
		if(!skipToQuestion(forceQ))
			return false;
	}

	// make sure a file is open, if it is, then we're good to move forward
	if (!md_files.size())
	{
		cout << "Invalid input file, line" << lineCount << "." << endl
			<< "No input file's given for DL_Poly, Lammps, or Tinker." << endl
			<< "Or wrong simulation software designated." << endl;
		return false;

	}

/////// forces Q
	if (!lastQuestions(forceQ, forces))
		return false;

	cout << "Forces value read successfully." << endl;

/////// energies Q
	if(!lastQuestions(energyQ, energies))
		return false;

	cout << "Energies value read successfully." << endl;

	findNextLine();

	while (line.empty())
	{
		findNextLine();
	}

	if (checkQ(stepQ))
	{
		findNextLine();

		while (line.empty())
		{
			findNextLine();
		}

		try
		{
			// if it can't parse this as a double, then it will cause the error to call;
			double temp = atof(line.c_str());
			step_size = line;
		}
		catch (...)
		{
			cout << "Invalid input file, line" << lineCount << " should be a percentage" << endl
				<< "e.g. \"10\" is 10%, \".1\" is 0.1%" << endl
				<< "without quotes at line " << lineCount << "." << endl;
			return false;
		}
	}
	else
	{
		cout << "Invalid input file, line" << lineCount << " should be \"" << stepQ << "\"" << endl
			<< "without quotes at line " << lineCount << "." << endl;
		return false;
	}

	return true;
}

bool parser::run()
{
	if (!openFile())
		return false;
	if (!setQMFiles())
		return false;
	if (!setQMSW())
		return false;
	if (!setGuessParam())
		return false;
	if (!setMDSim())
		return false;
	if (!finalQuestions())
		return false;

	cout << "ForceFit parser completed successfully." << endl;
	return true;
}

bool parser::openFile()
{
	cout << "Opening " << filename << "." << endl;
	inputFile.open(filename.c_str(), ifstream::in);		// open input for read

	if (!inputFile.is_open())
	{
		cout << "Error opening input file, please check the name." << endl
		<< "You can also run \"ForceFit create\", without quotes, to create a new input file" << endl;
		return false;
	}

	cout << filename << " opened successfully." << endl;
	return true;
}

void parser::updateFilename(string fn)
{
	filename = fn;
}

// PUBLIC functions

int parser::getQMSW()
{
	return qmsw;
}

int parser::getQType()
{
	return QType;
}

vector<string> parser::getQMFiles()
{
	return qmFiles;
}

vector<double> parser::getGuess()
{
	return guess;
}

vector<string> parser::getMDFiles()
{
	return md_files;
}

bool parser::getForces()
{
	return forces;
}

bool parser::getEnergies()
{
	return energies;
}

string parser::getStepSize()
{
	return step_size;
}

void parser::createNewTextFile()
{
	string filename;
	string cont;
	ofstream outfile;
	
	cout << "Enter Filename for the inputfile: " << endl;
	getline(cin, filename);

	if (!filename.length())
		filename = "forcefit.txt";

	while ( 1 )
	{
		cout << "If " << filename << " exists in the current directory it's contents will be emptied, " 
			<< endl << "Continue? (y/n) : ";
		getline(cin, cont);
		if ( !cont.compare("y") ) 
			break;
		else if ( !cont.compare("n") )
			break;
		else
			cout << "Invalid input, either lowercase y or lowercase n" << endl;
	}

	if ( !cont.compare("n") )
	{
		cout << "Did not create a new file." << endl
			<< "Exiting ForceFit." << endl;
		exit(1);
	}

	// truncate file if it exists
	outfile.open(filename.c_str());
	
	if (!outfile.is_open())
	{
		cout << "Error opening " << filename << "." << endl;
		exit(-1);
	}

	outfile << "List of QM files:\n\n\n"
	<< "QM software used to create these files:\n\n\n"
	<< "Guess parameters for force field:\n\n\n"
	<< "MD simulations software to be used during the fit:\n\n\n"
	<< "For Dl_Poly provide path to Control:\n\n\n"
	<< "For Dl_Poly provide path to Field:\n\n\n"
	<< "For Dl_poly provide path to Dl_Poly executable:\n\n\n"
	<< "For Lammps provide path to Input:\n\n\n"
	<< "For Lammps provide path to Data:\n\n\n"
	<< "For Lammps provide path to Lammps executable:\n\n\n"
	<< "For Tinker provide path to XYZ Template:\n\n\n"
	<< "For Tinker provide path to PRM:\n\n\n"
	<< "For Tinker provide path to Key:\n\n\n"
	<< "For Tinker provide path to Tinker executable:\n\n\n"
	<< "Are forces to be used in fitting:\n\n\n"
	<< "Are energies to be used in fitting:\n\n\n"
	<< "Step size:\n\n\n";

	cout << "Finished creating new input file \"" << filename
	<< "\" in ForceFit's directory." << endl << endl;

	cout << "After inputting information in \"" << filename
	<< "\" run forcefit as: " << endl
	<< "\"ForceFit " << filename << "\" without quotes." << endl
	<< "Exiting." << endl;

	outfile.close();

	exit(1);
}

parser::parser()
{
	filename = "forcefit.txt";
	lineCount = 0;
	QType = 0;
}

parser::~parser()
{
	ifstream *temp;
	// make sure to close all files we have opened, if we opened any
	if (inputFile.is_open())
	{
		inputFile.close();
	}
}
