#ifndef ALGORITHMH
#define ALGORITHMH
#include"CMatrix.h"
#include"PhyConst.h"
#include<cmath>
#include<iostream>
using namespace std;
typedef const double CD;
typedef const int CI;
typedef CMatrix_1D <double> M1D;
typedef CMatrix_2D <double> M2D;
typedef CMatrix_3D <double> M3D;

enum Rotation{ xyz, xzy, yxz, yzx, zxy, zyx };
// The value of IM in function SymmSplit 
// xyz means first along x, then along y, then along z
void ATsplit( M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D & JT, M3D & JA,\
	M3D &ED, M3D &BD, M3D &VX, M3D &VY, M3D & VZ, \
	CD & DT, CD & DX, CD & DY, CD & DZ, CD &Diff);

void AEvolve(CI & IW, M3D &PTT1, M3D &PTX1, M3D &PTY1, M3D &PTZ1, M3D &PJT1, M3D & PJA1,\
	M3D &VX, M3D &VY, M3D &VZ,CD &DT, CD &DX, CD &DY, CD &DZ, CD &Diff);


void Tsplit( M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D & JT, \
	M3D &ED, M3D &BD, M3D &VX, M3D &VY, M3D & VZ, \
	CD & DT, CD & DX, CD & DY, CD & DZ, CD &Diff);

void QuickTsplit( M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D & JT, \
	M3D &ED, M3D &BD, M3D &VX, M3D &VY, M3D & VZ, \
	CD & DT, CD & DX, CD & DY, CD & DZ, CD &Diff);

void Tsplit2( M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D & JT, \
	M3D &ED, M3D &BD, M3D &VX, M3D &VY, M3D & VZ, \
	CD & DT, CD & DX, CD & DY, CD & DZ, CD &DIFF);

void SymmSplit(CI & IM, M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D &JT, \
	M3D &PTT1, M3D &PTX1, M3D &PTY1, M3D &PTZ1, M3D &PJT1,\
	M3D &VX, M3D &VY, M3D &VZ,CD &DT, CD &DX, CD &DY, CD &DZ, CD &Diff);

void Evolve(CI & IW, M3D &PTT1, M3D &PTX1, M3D &PTY1, M3D &PTZ1, M3D &PJT1,\
	M3D &VX, M3D &VY, M3D &VZ,CD &DT, CD &DX, CD &DY, CD &DZ, CD &Diff);

void Shasta1D( M1D & A, M1D & U, M1D & B,\
	CI &ISL, CD &Tdw, CD &Diff);

#endif 
