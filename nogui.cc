#include "nogui.h"

// get each file with geometry sets into geoSets
void nogui::setGeometrySets()
{
	// read geometries from qmFilenames given
	for(vector<string>::iterator it = qmFilenames.begin();
	it != qmFilenames.end();
	it++)
	{
		reader->inputFiles.push_back(*it);
	}

	// get geometries
	set = reader->getGeometries();
	geoSets.push_back(reader->getGeometries());

	int i = 1;
	// give user output for how many geometries were loaded
	for (vector<Geometry>::iterator geomit = set.geometries.begin();
	geomit != set.geometries.end();
	geomit++)
	{
		cout << "Geometry " << i++ << " loaded." << endl; 
	}
}

void nogui::fillMDAnswers()
{
// get molecular dynamics questions
	int i = 0;
	mdQuestions = md->getQuestions();

	// fill in answers for mdQuestions
	// all answers will be MDFILENAME type, aka strings
	// dl_poly && lammps == 3 questions, tinker == 4 questions

	for ( vector<string>::iterator it = mdFilenames.begin();
	it != mdFilenames.end(); 
	it++ )
	{
		mdQuestions[i].stringAnswer = *it;
		i++;;
	}

}

void nogui::fillMinAnswers()
{
	// get minimizer questions
	minQuestions = min->getQuestions();

	// questions are always the same order, string, bool, bool
	(minQuestions)[0].stringAnswer = p.getStepSize();
	(minQuestions)[1].boolAnswer = p.getForces();
	(minQuestions)[2].boolAnswer = p.getEnergies();
}




// PUBLIC functions
void nogui::run()
{

	if (qmsw == GAUSSIAN)
	{
		cout << "Setting scanner to Gaussian." << endl;
		reader = new GaussianScan;
	}
	else if ( qmsw == CRYSTAL )
	{
		cout << "Setting scanner to Crystal." << endl;
		reader = new CrystalScan;
	}
	else if ( qmsw == NWCHEM )
	{
		cout << "Setting scanner to NWChem." << endl;
		reader = new NWChemScan;  
	}
	else
	{
		cout << "Could not find QM software type!" << endl;
		return;
	}

	if ( QType == DL_POLY )
	{
		cout << "MD is Dl_Poly" << endl;
		md = new Dl_Poly;
	}
	else if ( QType == LAMMPS )
	{
		cout << "MD is Lammps" << endl;
		md = new Lammps;
	}
	else if ( QType == TINKER )
	{
		cout << "MD is Tinker" << endl;
		md = new Tinker;
	}
	else
	{
		cout << "Parsing error, no qm software set!" << endl;
		return;
	}

	min = new Powell1965;

	// create GeometrySets (1)
	setGeometrySets();

	// create guess variables (2)
	// already done in the constructor

	// create mdQuestions (3)
	fillMDAnswers();

	// create minQuestions (4)
	fillMinAnswers();

	try{
		// prepare molecular dynamics
		// using the above (1-3)
		md->prepareMD(geoSets, variables, mdQuestions);
		// prepare minimizer
		// using the above (1, 4)
		min->prepareMin(geoSets, minQuestions);
		// run the program
		min->runMin(md);
	}
	catch(const MolecularDynamicsException & e)
	{
		cout << "Molecular Dynamics Error!" << endl;
		cout << e.what() << endl;
	}	//md->prepareMD(geoSets, variables, mdQuestions);

	// free memory
	delete reader;
	delete min;
	delete md;	
}
nogui::nogui(string s)
{
	if (!strlen(s.c_str()))
	{
		cout << "No filename given, defaulting to \"forcefit.txt\"" << endl;
	}
	else if (!s.compare("create"))
	{
		p.createNewTextFile();
	}
	else
	{
		p.updateFilename(s);
	}

	// run parser as soon as nogui is instantiated
	if (p.run())
	{
		variables = p.getGuess();
		qmFilenames = p.getQMFiles();
		mdFilenames = p.getMDFiles();
		qmsw = p.getQMSW();
		QType = p.getQType();
	}
	// parser failed, exit program
	else
	{
		cout << "Parser failed, exiting ForceFit." << endl;
		exit(-1);
	}
}

nogui::~nogui()
{

}
