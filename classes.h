#ifdef ScanReaderInclude
#include "scanReader/crystalScan.h"
#include "scanReader/gaussianScan.h"
#include "scanReader/nwchemScan.h"
#endif

#ifdef ScanReaderClass
ScanReaderClass(CrystalScan)
ScanReaderClass(GaussianScan)
ScanReaderClass(NWChemScan)
#endif
	
#ifdef GradientCreatorInclude
#include "gradientCreator/crystalGradient.h"
#endif

#ifdef GradientCreatorClass
GradientCreatorClass(CrystalGradient)
#endif

#ifdef MolecularDynamicsInclude
//#include "molecularDynamics/amber.h"
#include "molecularDynamics/dl_poly.h"
#include "molecularDynamics/lammps.h"
#include "molecularDynamics/tinker.h"
#endif

#ifdef MolecularDynamicsClass
//MolecularDynamicsClass(Amber)
MolecularDynamicsClass(Dl_Poly)
MolecularDynamicsClass(Lammps)
MolecularDynamicsClass(Tinker)
#endif

#ifdef MinimizerInclude
#include "minimizer/powell1965.h"
#endif

#ifdef MinimizerClass
MinimizerClass(Powell1965)
#endif
