CXX=g++
CXXFLAGS=-g
LDFLAGS=-g

TARGET=ForceFit

SRCS=atom.cc main.cc nogui.cc parser.cc geometry.cc geometrySet.cc scanReader/scanReader.cc scanReader/crystalScan.cc scanReader/gaussianScan.cc scanReader/nwchemScan.cc gradientCreator/crystalGradient.cc molecularDynamics/molecularDynamics.cc molecularDynamics/dl_poly.cc molecularDynamics/lammps.cc molecularDynamics/tinker.cc minimizer/minimizer.cc minimizer/powell1965.cc pugixml/pugixml.cpp
OBJS=${SRCS:.cc=.o}

all: ${TARGET}

ForceFit: ${OBJS}
	${CXX} ${LDFLAGS} -o $@ ${OBJS}

clean:
	-rm *.o */*.o ${TARGET}

ForceFit.o: classes.h

classes.h: scanReader/crystalScan.h
classes.h: scanReader/gaussianScan.h
classes.h: scanReader/nwchemScan.h
classes.h: gradientCreator/crystalGradient.h
classes.h: molecularDynamics/dl_poly.h
classes.h: molecularDynamics/lammps.h
classes.h: molecularDynamics/tinker.h
classes.h: minimizer/powell1965.h

#gui/setWindow.o: gui/setWindow.h gui/mdQuestionWindow.h gui/geometryWindow.h classes.h
#gui/addSet.o: gui/addSet.h gui/setWindow.h
#gui/mdQuestionWindow.o: gui/mdQuestionWindow.h
#gui/minQuestionWindow.o: gui/minQuestionWindow.h
#gui/geometryWindow.o: gui/geometryWindow.h

scanReader/crystalScan.h: scanReader/scanReader.h

scanReader/gaussianScan.h: scanReader/scanReader.h

scanReader/nwchemScan.h: scanReader/scanReader.h

gradientCreator/crystalGradient.h: gradientCreator/gradientCreator.h

molecularDynamics/dl_poly.h: molecularDynamics/molecularDynamics.h

molecularDynamics/lammps.h: molecularDynamics/molecularDynamics.h

molecularDynamics/tinker.h: molecularDynamics/molecularDynamics.h

minimizer/powell1965.h: minimizer/minimizer.h
