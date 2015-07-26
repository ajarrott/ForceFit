#include <iostream>
#include "nogui.h"

int main(int argc, char* argv[])
{
	if ( argc == 1 || argc > 2 )
	{
		cout << "ForceFit usage: " << endl;
		cout << "ForceFit \"inputfilename\" without quotes." << endl;
		cout << "OR" << endl;
		cout << "ForceFit create" << endl;
		exit(1);
	}
	nogui n(argv[1]);

	n.run();

	return 0;
}
