ran2.o:ran2.h ran2.cpp
	g++  -O3  -c ran2.cpp
CCube.o:CCube.h CCube.cpp
	g++  -O3  -c CCube.cpp
CEos.o:CEos.h CEos.cpp
	g++  -O3  -c CEos.cpp
SolveU_mu.o:SolveU_mu.h SolveU_mu.cpp
	g++  -O3  -c SolveU_mu.cpp
Algorithm.o:Algorithm.h Algorithm.cpp
	g++  -O3  -c Algorithm.cpp
CHydro.o:CHydro.h CHydro.cpp 
	g++  -O3  -c CHydro.cpp
CAAcollision.o:CAAcollision.h CAAcollision.cpp
	g++  -O3  -c CAAcollision.cpp
main.o:main.cpp
	g++  -O3  -c main.cpp
hydro:main.o CHydro.o CAAcollision.o CEos.o Algorithm.o CCube.o SolveU_mu.o ran2.o 
	g++  -O3  CCube.o Algorithm.o main.o CHydro.o CAAcollision.o CEos.o SolveU_mu.o ran2.o -o hydro 

clean:
	rm -f *.o matrix hydro myjob.sh.* core* script/*
