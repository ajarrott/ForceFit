--------------
  7/25/2015
--------------
On error stopped program from deleting the directory that had useful output files in it.

--------------
  3/14/2015
--------------
This new version of ForceFit has no GUI, so there should no longer be any problems with compiling the graphics libraries.

--------------
 INSTALLATION
--------------
Open a terminal and type make.  If you have all the necessary scientific libraries there should be no issues compiling.

--------------
    USAGE
--------------
Initally to use the new ForceFit create a new input file by using the command "ForceFit create".  You can then designate a filename and it will open that file and input all the information you need to fill out.  Be careful though, it will truncate (overwrite) a file that already exists with the same name in ForceFit's directory.

In this new file DO NOT delete questions that are unneccesary to you, even though they are unused the program still ensures the input file is valid from these questions.  Leave it all the same and just input information where it is appropriate for your use.  Leave at least one empty line after an answer and before a question

For example:
List of QM files:
/home/user/path/to/gradient.out

QM software used to create these files:
NWChem

Guess parameters for force field:
...

If you want to easily comment out some guess variables or QM files just use the character '#' at the very beginning of the line.  For example:

# This line is commented out
 # This line is not

When designating files within the input file use the full path to ensure that the correct file is being chosen.
For Example:

For Tinker provide path to XYZ Template:
/home/user/path/to/myfile.xyz

For guess variables only the number after an '=' will be used, the variable name does not matter as long as it contains only characters.  Spaces or tabs around the '=' will be removed
For Example:

Guess parameters for force field:
A=.504444
B = -2.71092
Cat	 = 		1.24325

The above are all valid.

DO NOT just type a value, it does need a designated variable name and an equals sign.  The program will alert you of this error.

--------------
  3/20/2015
--------------
The following are updates made to Forcefit.

Made user experience nicer by allowing user to have characters after the question be ignored.
Still make sure to write your information on the next line after the question.

For example:
List of QM files:   aflksdj;f salk;dsjfoiejoijoij98374
forcefit.txt

The above is valid, but it still needs to have the correct question so this is NOT valid for lines 1 and 2:
QM files:
forcefit.txt

--------------
 USING CREATE
--------------
Now you can now just press return when using "ForceFit create" it will default a text file to forcefit.txt

--------------
  DEBUG HELP
--------------
If you recieve output like this after everything has been parsed:

RMSD: -nan
X: -nan Y: -nan Z: -nan
Center Forces Us:
Center Energy Us:
Force	q	s:
Energy	q	s:
U: 0 Variables: (variable1), (variable2), .... // where these are your designated variables
Segmentation fault

There is likely a problem with your input files, the file may exist but it is not the correct type for the MD simulation software you are using.  Make sure you have the correct type of MD Simulation and the correct filenames.

--------------
  DIRECTORY
    CHECK
--------------
The program now checks to make sure that the file is NOT a directory, to prevent possibility of accidentally using directories as filenames.

