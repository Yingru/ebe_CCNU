#include"Algorithm.h"
#include<iomanip>
#include<iostream>
using namespace std;
static long int STEP_XYZ_OR_YXZ=0;

void Shasta1D( M1D & A, M1D & U, M1D & B, CI &ISL, CD &TDW, CD &DIFF);
void AEvolve(CI & IW, M3D &PTT1, M3D &PTX1, M3D &PTY1, M3D &PTZ1, M3D &PJT1, M3D &PJA1,\
	M3D &VX, M3D &VY, M3D &VZ,CD &DT, CD &DX, CD &DY, CD &DZ, CD &DIFF)
{
    double TDW;
    int NW0,NW,N1L,N1H,N2L,N2H,N3L,N3H;
    switch(IW){
	case 1:
	    TDW=DT/DX;
	    NW0=PTT1.pag_l;
	    NW =PTT1.pag_h;
	    N1L=PTT1.row_l;
	    N1H=PTT1.row_h;
	    N2L=PTT1.col_l;
	    N2H=PTT1.col_h;
	    N3L=PTT1.pag_l;
	    N3H=PTT1.pag_h;
	    break;
	case 2:
	    TDW=DT/DY;
	    NW0=PTT1.row_l;
	    NW =PTT1.row_h;
	    N1L=PTT1.pag_l;
	    N1H=PTT1.pag_h;
	    N2L=PTT1.col_l;
	    N2H=PTT1.col_h;
	    N3L=PTT1.row_l;
	    N3H=PTT1.row_h;
	    break;

	case 3:
	    TDW=DT/DZ;
	    NW0=PTT1.col_l;
	    NW =PTT1.col_h;
	    N1L=PTT1.pag_l;
	    N1H=PTT1.pag_h;
	    N2L=PTT1.row_l;
	    N2H=PTT1.row_h;
	    N3L=PTT1.col_l;
	    N3H=PTT1.col_h;
	    break;
	default:
	    break;
    }
    M1D TT0(NW0,NW);
    M1D TX0(NW0,NW);
    M1D TY0(NW0,NW);
    M1D TE0(NW0,NW);
    M1D RJ0(NW0,NW);
    M1D JA0(NW0,NW);
    M1D VWN(NW0,NW);

    M1D ScT0(NW0,NW);
    M1D ScT1(NW0,NW);
    M1D ScT2(NW0,NW);
    M1D ScT3(NW0,NW);
    M1D ScT4(NW0,NW);
    M1D ScT5(NW0,NW);


    for(int I1=N1L; I1<=N1H; I1++)
	for(int I2=N2L; I2<=N2H; I2++){
	    for(int I3=N3L; I3<=N3H; I3++){
		switch(IW){
		    case 1:   //Along x
			TT0(I3) = PTT1(I3,I1,I2);
			TX0(I3) = PTX1(I3,I1,I2);
			TY0(I3) = PTY1(I3,I1,I2);
			TE0(I3) = PTZ1(I3,I1,I2);
			RJ0(I3) = PJT1(I3,I1,I2);
			JA0(I3) = PJA1(I3,I1,I2);
			VWN(I3) =   VX(I3,I1,I2);
			break;
		    case 2:   //Along y
			TT0(I3) = PTT1(I1,I3,I2);
			TX0(I3) = PTX1(I1,I3,I2);
			TY0(I3) = PTY1(I1,I3,I2);
			TE0(I3) = PTZ1(I1,I3,I2);
			RJ0(I3) = PJT1(I1,I3,I2);
			JA0(I3) = PJA1(I1,I3,I2);
			VWN(I3) =   VY(I1,I3,I2);
			break;
		    case 3:   //Along eta
			TT0(I3) = PTT1(I1,I2,I3);
			TX0(I3) = PTX1(I1,I2,I3);
			TY0(I3) = PTY1(I1,I2,I3);
			TE0(I3) = PTZ1(I1,I2,I3);
			RJ0(I3) = PJT1(I1,I2,I3);
			JA0(I3) = PJA1(I1,I2,I3);
			VWN(I3) =   VZ(I1,I2,I3);
			break;
		    default:
			break;
		}//End switch
	    }//End for I3

	    Shasta1D( TT0, VWN, ScT0, 0, TDW, DIFF);
	    Shasta1D( TX0, VWN, ScT1, 0, TDW, DIFF);
	    Shasta1D( TY0, VWN, ScT2, 0, TDW, DIFF);
	    Shasta1D( TE0, VWN, ScT3, 0, TDW, DIFF);
	    Shasta1D( RJ0, VWN, ScT4, 0, TDW, DIFF);
	    Shasta1D( JA0, VWN, ScT5, 0, TDW, DIFF);

	    for(int I3=N3L; I3<=N3H; I3++){
		switch(IW){
		    case 1:
			PTT1(I3,I1,I2)=TT0(I3); 
			PTX1(I3,I1,I2)=TX0(I3); 
			PTY1(I3,I1,I2)=TY0(I3); 
			PTZ1(I3,I1,I2)=TE0(I3); 
			PJT1(I3,I1,I2)=RJ0(I3); 
			PJA1(I3,I1,I2)=JA0(I3); 
			break;
		    case 2:
			PTT1(I1,I3,I2)= TT0(I3); 
			PTX1(I1,I3,I2)= TX0(I3); 
			PTY1(I1,I3,I2)= TY0(I3); 
			PTZ1(I1,I3,I2)= TE0(I3); 
			PJT1(I1,I3,I2)= RJ0(I3); 
			PJA1(I1,I3,I2)= JA0(I3); 
			break;
		    case 3:
			PTT1(I1,I2,I3)=TT0(I3); 
			PTX1(I1,I2,I3)=TX0(I3); 
			PTY1(I1,I2,I3)=TY0(I3); 
			PTZ1(I1,I2,I3)=TE0(I3); 
			PJT1(I1,I2,I3)=RJ0(I3); 
			PJA1(I1,I2,I3)=JA0(I3); 
			break;
		    default:
			break;
		}
	    }
	} 
}


//Time splitting method to solve 3D hydro evolution
void ATsplit( M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D & JT, M3D & JA, \
	M3D &ED, M3D &BD, M3D &VX, M3D &VY, M3D & VZ, \
	CD & DT, CD & DX, CD & DY, CD & DZ, CD &DIFF)
{
    const int NX0=TT.pag_l;
    const int NX =TT.pag_h;
    const int NY0=TT.row_l;
    const int NY =TT.row_h;
    const int NZ0=TT.col_l;
    const int NZ =TT.col_h;

    M3D PTT1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTX1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTY1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTZ1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJT1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJA1(NX0,NX, NY0,NY, NZ0,NZ);

    M3D PTT2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTX2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTY2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTZ2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJT2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJA2(NX0,NX, NY0,NY, NZ0,NZ);

    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		PTT1(I,J,K)  = TT(I,J,K);
		PTX1(I,J,K)  = TX(I,J,K);
		PTY1(I,J,K)  = TY(I,J,K);
		PTZ1(I,J,K)  = TZ(I,J,K);
		PJT1(I,J,K)  = JT(I,J,K);
		PJA1(I,J,K)  = JA(I,J,K);

		PTT2(I,J,K)  = TT(I,J,K);
		PTX2(I,J,K)  = TX(I,J,K);
		PTY2(I,J,K)  = TY(I,J,K);
		PTZ2(I,J,K)  = TZ(I,J,K);
		PJT2(I,J,K)  = JT(I,J,K);
		PJA2(I,J,K)  = JA(I,J,K);
	    }
    AEvolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, PJA1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    AEvolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, PJA1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    AEvolve(2, PTT2, PTX2, PTY2, PTZ2, PJT2, PJA2, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    AEvolve(1, PTT2, PTX2, PTY2, PTZ2, PJT2, PJA2, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		PTT2(I,J,K)  =(PTT1(I,J,K)+PTT2(I,J,K))/2.0;   
		PTX2(I,J,K)  =(PTX1(I,J,K)+PTX2(I,J,K))/2.0;
		PTY2(I,J,K)  =(PTY1(I,J,K)+PTY2(I,J,K))/2.0;
		PTZ2(I,J,K)  =(PTZ1(I,J,K)+PTZ2(I,J,K))/2.0;
		PJT2(I,J,K)  =(PJT1(I,J,K)+PJT2(I,J,K))/2.0;
		PJA2(I,J,K)  =(PJA1(I,J,K)+PJA2(I,J,K))/2.0;
	    }
    AEvolve(3, PTT2, PTX2, PTY2, PTZ2, PJT2, PJA2, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		TT(I,J,K)  =PTT2(I,J,K);   
		TX(I,J,K)  =PTX2(I,J,K);
		TY(I,J,K)  =PTY2(I,J,K);
		TZ(I,J,K)  =PTZ2(I,J,K);
		JT(I,J,K)  =PJT2(I,J,K);
		JA(I,J,K)  =PJA2(I,J,K);
	    }

}//End ATsplit();

void Evolve(CI & IW, M3D &PTT1, M3D &PTX1, M3D &PTY1, M3D &PTZ1, M3D &PJT1,\
	M3D &VX, M3D &VY, M3D &VZ,CD &DT, CD &DX, CD &DY, CD &DZ, CD &DIFF)
{
/* Just need to paralize this subroutine     */
    double TDW;
    int NW0,NW,N1L,N1H,N2L,N2H,N3L,N3H;
    switch(IW){
	case 1:
	    TDW=DT/DX;
	    NW0=PTT1.pag_l;
	    NW =PTT1.pag_h;
	    N1L=PTT1.row_l;
	    N1H=PTT1.row_h;
	    N2L=PTT1.col_l;
	    N2H=PTT1.col_h;
	    N3L=PTT1.pag_l;
	    N3H=PTT1.pag_h;
	    break;
	case 2:
	    TDW=DT/DY;
	    NW0=PTT1.row_l;
	    NW =PTT1.row_h;
	    N1L=PTT1.pag_l;
	    N1H=PTT1.pag_h;
	    N2L=PTT1.col_l;
	    N2H=PTT1.col_h;
	    N3L=PTT1.row_l;
	    N3H=PTT1.row_h;
	    break;

	case 3:
	    TDW=DT/DZ;
	    NW0=PTT1.col_l;
	    NW =PTT1.col_h;
	    N1L=PTT1.pag_l;
	    N1H=PTT1.pag_h;
	    N2L=PTT1.row_l;
	    N2H=PTT1.row_h;
	    N3L=PTT1.col_l;
	    N3H=PTT1.col_h;
	    break;
	default:
	    break;
    }
    M1D TT0(NW0,NW);
    M1D TX0(NW0,NW);
    M1D TY0(NW0,NW);
    M1D TE0(NW0,NW);
    M1D RJ0(NW0,NW);
    M1D VWN(NW0,NW);

    M1D ScT0(NW0,NW);
    M1D ScT1(NW0,NW);
    M1D ScT2(NW0,NW);
    M1D ScT3(NW0,NW);
    M1D ScT4(NW0,NW);


    for(int I1=N1L; I1<=N1H; I1++)
	for(int I2=N2L; I2<=N2H; I2++){
	    for(int I3=N3L; I3<=N3H; I3++){
		switch(IW){
		    case 1:   //Along x
			TT0(I3) = PTT1(I3,I1,I2);
			TX0(I3) = PTX1(I3,I1,I2);
			TY0(I3) = PTY1(I3,I1,I2);
			TE0(I3) = PTZ1(I3,I1,I2);
			RJ0(I3) = PJT1(I3,I1,I2);
			VWN(I3) =   VX(I3,I1,I2);
			break;
		    case 2:   //Along y
			TT0(I3) = PTT1(I1,I3,I2);
			TX0(I3) = PTX1(I1,I3,I2);
			TY0(I3) = PTY1(I1,I3,I2);
			TE0(I3) = PTZ1(I1,I3,I2);
			RJ0(I3) = PJT1(I1,I3,I2);
			VWN(I3) =   VY(I1,I3,I2);
			break;
		    case 3:   //Along eta
			TT0(I3) = PTT1(I1,I2,I3);
			TX0(I3) = PTX1(I1,I2,I3);
			TY0(I3) = PTY1(I1,I2,I3);
			TE0(I3) = PTZ1(I1,I2,I3);
			RJ0(I3) = PJT1(I1,I2,I3);
			VWN(I3) =   VZ(I1,I2,I3);
			break;
		    default:
			break;
		}//End switch
	    }//End for I3

	    Shasta1D( TT0, VWN, ScT0, 0, TDW, DIFF);
	    Shasta1D( TX0, VWN, ScT1, 0, TDW, DIFF);
	    Shasta1D( TY0, VWN, ScT2, 0, TDW, DIFF);
	    Shasta1D( TE0, VWN, ScT3, 0, TDW, DIFF);
	    Shasta1D( RJ0, VWN, ScT4, 0, TDW, DIFF);

	    for(int I3=N3L; I3<=N3H; I3++){
		switch(IW){
		    case 1:
			PTT1(I3,I1,I2)=TT0(I3); 
			PTX1(I3,I1,I2)=TX0(I3); 
			PTY1(I3,I1,I2)=TY0(I3); 
			PTZ1(I3,I1,I2)=TE0(I3); 
			PJT1(I3,I1,I2)=RJ0(I3); 
			break;
		    case 2:
			PTT1(I1,I3,I2)= TT0(I3); 
			PTX1(I1,I3,I2)= TX0(I3); 
			PTY1(I1,I3,I2)= TY0(I3); 
			PTZ1(I1,I3,I2)= TE0(I3); 
			PJT1(I1,I3,I2)= RJ0(I3); 
			break;
		    case 3:
			PTT1(I1,I2,I3)=TT0(I3); 
			PTX1(I1,I2,I3)=TX0(I3); 
			PTY1(I1,I2,I3)=TY0(I3); 
			PTZ1(I1,I2,I3)=TE0(I3); 
			PJT1(I1,I2,I3)=RJ0(I3); 
			break;
		    default:
			break;
		}
	    }
	} 
}
//Time splitting method to solve 3D hydro evolution using xyz yxz iteration
void QuickTsplit( M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D & JT, \
	M3D &ED, M3D &BD, M3D &VX, M3D &VY, M3D & VZ, \
	CD & DT, CD & DX, CD & DY, CD & DZ, CD &DIFF)
{
    const int NX0=TT.pag_l;
    const int NX =TT.pag_h;
    const int NY0=TT.row_l;
    const int NY =TT.row_h;
    const int NZ0=TT.col_l;
    const int NZ =TT.col_h;

    M3D PTT1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTX1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTY1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTZ1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJT1(NX0,NX, NY0,NY, NZ0,NZ);

    M3D PTT2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTX2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTY2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTZ2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJT2(NX0,NX, NY0,NY, NZ0,NZ);

    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		PTT1(I,J,K)  = TT(I,J,K);
		PTX1(I,J,K)  = TX(I,J,K);
		PTY1(I,J,K)  = TY(I,J,K);
		PTZ1(I,J,K)  = TZ(I,J,K);
		PJT1(I,J,K)  = JT(I,J,K);
	    }
        if((STEP_XYZ_OR_YXZ/2)%2){
    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    }
    else{
    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    }
    STEP_XYZ_OR_YXZ ++;

    Evolve(3, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);

    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		TT(I,J,K)  =PTT1(I,J,K);   
		TX(I,J,K)  =PTX1(I,J,K);
		TY(I,J,K)  =PTY1(I,J,K);
		TZ(I,J,K)  =PTZ1(I,J,K);
		JT(I,J,K)  =PJT1(I,J,K);
	    }

}

//Time splitting method to solve 3D hydro evolution
void Tsplit( M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D & JT, \
	M3D &ED, M3D &BD, M3D &VX, M3D &VY, M3D & VZ, \
	CD & DT, CD & DX, CD & DY, CD & DZ, CD &DIFF)
{
    const int NX0=TT.pag_l;
    const int NX =TT.pag_h;
    const int NY0=TT.row_l;
    const int NY =TT.row_h;
    const int NZ0=TT.col_l;
    const int NZ =TT.col_h;

    M3D PTT1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTX1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTY1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTZ1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJT1(NX0,NX, NY0,NY, NZ0,NZ);

    M3D PTT2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTX2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTY2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTZ2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJT2(NX0,NX, NY0,NY, NZ0,NZ);

    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		PTT1(I,J,K)  = TT(I,J,K);
		PTX1(I,J,K)  = TX(I,J,K);
		PTY1(I,J,K)  = TY(I,J,K);
		PTZ1(I,J,K)  = TZ(I,J,K);
		PJT1(I,J,K)  = JT(I,J,K);

		PTT2(I,J,K)  = TT(I,J,K);
		PTX2(I,J,K)  = TX(I,J,K);
		PTY2(I,J,K)  = TY(I,J,K);
		PTZ2(I,J,K)  = TZ(I,J,K);
		PJT2(I,J,K)  = JT(I,J,K);
	    }
    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    Evolve(2, PTT2, PTX2, PTY2, PTZ2, PJT2, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    Evolve(1, PTT2, PTX2, PTY2, PTZ2, PJT2, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		PTT2(I,J,K)  =(PTT1(I,J,K)+PTT2(I,J,K))/2.0;   
		PTX2(I,J,K)  =(PTX1(I,J,K)+PTX2(I,J,K))/2.0;
		PTY2(I,J,K)  =(PTY1(I,J,K)+PTY2(I,J,K))/2.0;
		PTZ2(I,J,K)  =(PTZ1(I,J,K)+PTZ2(I,J,K))/2.0;
		PJT2(I,J,K)  =(PJT1(I,J,K)+PJT2(I,J,K))/2.0;
	    }
    Evolve(3, PTT2, PTX2, PTY2, PTZ2, PJT2, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		TT(I,J,K)  =PTT2(I,J,K);   
		TX(I,J,K)  =PTX2(I,J,K);
		TY(I,J,K)  =PTY2(I,J,K);
		TZ(I,J,K)  =PTZ2(I,J,K);
		JT(I,J,K)  =PJT2(I,J,K);
	    }

}

void Tsplit2( M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D & JT, \
	M3D &ED, M3D &BD, M3D &VX, M3D &VY, M3D & VZ, \
	CD & DT, CD & DX, CD & DY, CD & DZ, CD &DIFF)
{	const int NX0=TT.pag_l;
    const int NX =TT.pag_h;
    const int NY0=TT.row_l;
    const int NY =TT.row_h;
    const int NZ0=TT.col_l;
    const int NZ =TT.col_h;


    M3D PTT1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTX1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTY1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTZ1(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJT1(NX0,NX, NY0,NY, NZ0,NZ);

    M3D PTT2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTX2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTY2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PTZ2(NX0,NX, NY0,NY, NZ0,NZ);
    M3D PJT2(NX0,NX, NY0,NY, NZ0,NZ);


    for(int IM=xyz; IM<=zyx; IM++){
	SymmSplit( IM, TT, TX, TY, TZ, JT, \
		PTT1, PTX1, PTY1, PTZ1, PJT1,\
		VX, VY, VZ,DT, DX, DY, DZ, DIFF);
	for(int I=NX0; I<=NX; I++)
	    for(int J=NY0; J<=NY; J++)
		for(int K=NZ0; K<=NZ; K++){
		    PTT2(I,J,K)  =PTT2(I,J,K)+PTT1(I,J,K)/6.0;   
		    PTX2(I,J,K)  =PTX2(I,J,K)+PTX1(I,J,K)/6.0;
		    PTY2(I,J,K)  =PTY2(I,J,K)+PTY1(I,J,K)/6.0;
		    PTZ2(I,J,K)  =PTZ2(I,J,K)+PTZ1(I,J,K)/6.0;
		    PJT2(I,J,K)  =PJT2(I,J,K)+PJT1(I,J,K)/6.0;
		}
    }
}



void SymmSplit(CI & IM, M3D &TT, M3D &TX, M3D &TY, M3D &TZ, M3D &JT, \
	M3D &PTT1, M3D &PTX1, M3D &PTY1, M3D &PTZ1, M3D &PJT1,\
	M3D &VX, M3D &VY, M3D &VZ,CD &DT, CD &DX, CD &DY, CD &DZ, CD &DIFF)
{
    const int NX0=TT.pag_l;
    const int NX =TT.pag_h;
    const int NY0=TT.row_l;
    const int NY =TT.row_h;
    const int NZ0=TT.col_l;
    const int NZ =TT.col_h;
    for(int I=NX0; I<=NX; I++)
	for(int J=NY0; J<=NY; J++)
	    for(int K=NZ0; K<=NZ; K++){
		PTT1(I,J,K)  = TT(I,J,K);
		PTX1(I,J,K)  = TX(I,J,K);
		PTY1(I,J,K)  = TY(I,J,K);
		PTZ1(I,J,K)  = TZ(I,J,K);
		PJT1(I,J,K)  = JT(I,J,K);
	    }

    switch(IM){/*Take an average on different evolution orders */
	case xyz:
	    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(3, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    break;
	case xzy:
	    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(3, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    break;
	case yxz:
	    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(3, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    break;
	case yzx:
	    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(3, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    break;
	case zxy:
	    Evolve(3, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    break;
	case zyx:
	    Evolve(3, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(2, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    Evolve(1, PTT1, PTX1, PTY1, PTZ1, PJT1, VX, VY, VZ, DT, DX, DY, DZ, DIFF);
	    break;
	default:
	    break;
    }

}


void Shasta1D( M1D & A, M1D & U, M1D & B,\
	CI &ISL, CD &TDW, CD &DIFF)
{/*Zalesak's paper */


    double EIM1, EII, EIP1, QUP,QUM,AIM1,AI,BIM1,BI;
    double AIMM,AIM,AIP,AIPP,UIMM,UIM,UI,UIP,UIPP,VIMM,VIM,VI,VIP,VIPP,\
	WMAXIM,WMAXI,WMAXIP,WMINIM,WMINI,WMINIP,\
	PPIM,QPIM,RPIM,RPI,RPIP,RMIM,RMI,RMIP,PPI,QPI,PPIP,QPIP,PMIM,QMIM,PMI,QMI,\
	CIP,CIM,PMIP,QMIP,ACIM,ACIP;
    const int NW0=A.col_l;
    const int NW=A.col_h;
    int IM,IMM, IP,IPP;

    M1D AB(NW0, NW);
    M1D AC(NW0, NW);

    for(int I=NW0+1; I<=NW-1; I++){
	EIM1=U(I-1)*TDW;
	EII =U(I)*TDW;
	EIP1=U(I+1)*TDW;

   
    if(EIM1>0.5){
            cout<<"TDW="<<TDW<<endl;
            cout<<"Pos="<<I<<endl;
    }
	assert(EIM1<=0.5);

	QUP=(0.5-EII)/(1.0+(EIP1-EII));
	QUM=(0.5+EII)/(1.0-(EIM1-EII));
	AIM1 = A(I) - A(I-1);
	AI   = A(I+1) - A(I);
	BIM1 = ISL*TDW*(B(I)-B(I-1));
	BI   = ISL*TDW*(B(I+1)-B(I));
	AB(I)=0.5*(QUP*QUP*AI-QUM*QUM*AIM1) +(QUP+QUM)*A(I) -0.5*ISL*TDW*(B(I+1)-B(I-1));
    }

    //AB(NW0)=2.0*AB(NW0+1)-AB(NW0+2);
    //AB(NW) =2.0*AB(NW-1)-AB(NW-2);

    //AB(NW0)=0.0;
    //AB(NW) =0.0;
    AB(NW0)=AB(NW0+1);
    AB(NW)=AB(NW-1);
    //-------- End of transport and diffuse stage

    for(int I=NW0; I<=NW; I++){
	AC(I)=A(I);
    }

    for(int I=NW0+2; I<=NW-2; I++){
	IMM=I-2;
	IM=I-1;
	IP=I+1;
	IPP=I+2;

	AIMM= DIFF*(AB(I-1)-AB(I-2));
	AIM = DIFF*(AB(I)-AB(I-1));
	AIP = DIFF*(AB(I+1)-AB(I));
	AIPP= DIFF*(AB(I+2)-AB(I+1));


	UIMM = max(AC(IMM), AB(IMM));
	UIM = max(AC(IM), AB(IM));
	UI = max(AC(I), AB(I));
	UIP = max(AC(IP), AB(IP));
	UIPP = max(AC(IPP), AB(IPP));

	VIMM = min(AC(IMM), AB(IMM));
	VIM = min(AC(IM), AB(IM));
	VI = min(AC(I), AB(I));
	VIP = min(AC(IP), AB(IP));
	VIPP = min(AC(IPP), AB(IPP));

	WMAXIM = max(max(UIMM,UIM),UI);
	WMAXI=max(max(UIM, UI), UIP);
	WMAXIP = max(max(UI, UIP),UIPP);

	WMINIM=min(min(VIMM, VIM), VI);
	WMINI=min(min(VIM, VI), VIP);
	WMINIP=min(min(VI, VIP), VIPP);

	PPIM = max(AIMM,0.0)-min(AIM,0.0);
	QPIM = WMAXIM - AB(IM);

	if(PPIM>0.0) 
	    RPIM = min( QPIM/PPIM, 1.0 );
	else RPIM = 0.0;


	PPI = max(AIM,0.0)-min(AIP,0.0);
	QPI = WMAXI - AB(I);
	if(PPI>0.0) RPI = min( QPI/PPI, 1.0 );
	else RPI = 0.0;

	PPIP = max(AIP,0.0)-min(AIPP,0.0);
	QPIP = WMAXIP - AB(IP);
	if(PPIP>0.0)
	    RPIP = min( QPIP/PPIP, 1.0 );
	else
	    RPIP = 0.0;

	PMIM = max(AIM,0.0)-min(AIMM,0.0);
	QMIM = AB(IM) - WMINIM;
	if(PMIM>0.0)
	    RMIM = min( QMIM/PMIM, 1.0);
	else 
	    RMIM = 0.0;


	PMI = max(AIP,0.0)-min(AIM,0.0);
	QMI = AB(I) - WMINI;
	if(PMI>0.0) 
	    RMI = min( QMI/PMI, 1.0);
	else RMI = 0.0;

	PMIP = max(AIPP,0.0)-min(AIP,0.0);
	QMIP = AB(IP) - WMINIP;
	if(PMIP>0.0)  
	    RMIP = min( QMIP/PMIP, 1.0);
	else RMIP = 0.0;

	if (AIP>=0.0) 
	    CIP = min(RPIP, RMI);
	else 
	    CIP = min(RPI, RMIP);

	if (AIM>=0.0) 
	    CIM = min(RPI, RMIM);
	else
	    CIM = min(RPIM, RMI);

	ACIM = CIM*AIM;
	ACIP = CIP*AIP;

	A(I) = AB(I) - (ACIP - ACIM);
    }

    //-------- Flux Corrected St>=------------
    //A(NW0+1)=2.0*A(NW0+2)-A(NW0+3);
    //A(NW-1) =2.0*A(NW-2)-A(NW-3);

    //A(NW0)=2.0*A(NW0+2)-A(NW0+4);
    //A(NW) =2.0*A(NW-2)-A(NW-4);

    //A(NW0+1)=0.0;
    //A(NW-1) =0.0;

    //A(NW0)=0.0;
    //A(NW) =0.0;

    A(NW0+1)=A(NW0+2);
    A(NW-1) =A(NW-2);

    A(NW0)=A(NW0+2);
    A(NW) =A(NW-2);
}



//int main()
//{
//	CMatrix_1D <double> A(-100,100),U(-100,100),B(-100,100);
//	for(int i=-100; i<=100; i++)
//	{
//		double x=i*0.2;
//		double wd=0.8;
//		if(i>=-50 && i<=50) A(i)=2.0;
//		else A(i)=1.0;
//
//		if(i>0) U(i)=0.5;
//		else if(i<0)U(i)=-0.5;
//		else U(i)=0.0;
//	}
//
//	double t=0.6;
//	while(t<=8.48){
//		for(int i=-100; i<=100; i++){
//			double x=double(i)*0.2;
//			cout<<setprecision(8)<<setw(25)<<A(i)<<endl;
//		}
//		Shasta1D(A,U,B,0,0.2,0.125);
//		t=t+0.04;
//	}
//	return 0;
//}

