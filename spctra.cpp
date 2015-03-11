#include"Int.h"
#include<cmath>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
const int NRY=41;
const int NPT=15;
const int NPHI=48;
using namespace std;

struct CParticle{
    int monval;
    string  name;
	double	mass;
	double	width;
	int	gspin;      /* spin degeneracy */
	int	baryon;
	int	strange;
	int	charm;
	int	bottom;
	int	gisospin;  /* isospin degeneracy */
	int	charge;
	int	decays;    /* amount of decays listed for this resonance */
	int	stable;     /* defines whether this particle is considered as stable */
};


int main(int argc, char** argv)
{
assert(argc==2);
int MC=atoi(argv[1]);
char fname[64];
if(MC>=0) sprintf(fname,"results/dNdYPtdPtdPhi_%d.dat",abs(MC));
else if(MC<0) sprintf(fname,"results/dNdYPtdPtdPhi_A%d.dat",abs(MC));

double dNdYPtdPtdPhi[NRY][NPT][NPHI];
double *pg, *wg; //gauss position and weight
double *pl, *wl; //laguree position and weight

pg = gaulep48; wg = gaulew48;
pl = gala15x;  wl = gala15w;

for(int k=0; k!=NRY; k++)
for(int i=0; i!=NPT; i++)
for(int j=0; j!=NPHI; j++){



}


return 0;
}


double	Edndp3(double eta, double pt, double phi, CParticle par)
{
    double yr, ptr, phir;
    double	phir, val;
    double        f1, f2;
    int     	pn, nry, npt, nphi;

    if(phirin < 0.0){
        printf("ERROR: phir %15.8le < 0 !!! \n", phirin);exit(0);}
    if(phirin > 2.0*PI){
        printf("ERROR: phir %15.8le > 2PI !!! \n", phirin);exit(0);}
    if(phirin < 0.5*PI) phir = phirin;
    else{
        if(phirin < PI) phir = PI - phirin;
        else{
            if(phirin < 1.5*PI) phir = phirin - PI;
            else phir = 2.0*PI - phirin;
        }
    }

    pn = partid[MHALF + res_num];

    nry = 1; 
    while((yr > RY[nry])&&(nry<(particle[pn].nry-1))) nry++; 
    nphi = 1; 
    while((phir > PHI[nphi])&&(nphi<(particle[pn].nphi-1))) nphi++; 
    npt = 1; 
    while((ptr > particle[pn].pt[npt]) &&
            (npt<(particle[pn].npt - 1))) npt++; 

    /* phi and rapidity interpolation */

    double fy1, fy2;

    fy1 = lin_int(PHI[nphi-1], PHI[nphi], 
            particle[pn].dNdyptdptdphi[nry-1][npt-1][nphi-1], 
            particle[pn].dNdyptdptdphi[nry-1][npt-1][nphi], phir);
    fy2 = lin_int(PHI[nphi-1], PHI[nphi], 
            particle[pn].dNdyptdptdphi[nry][npt-1][nphi-1], 
            particle[pn].dNdyptdptdphi[nry][npt-1][nphi], phir);

    f1 = lin_int(RY[nry-1], RY[nry], fy1, fy2, yr);

    fy1 = lin_int(PHI[nphi-1], PHI[nphi], 
            particle[pn].dNdyptdptdphi[nry-1][npt][nphi-1], 
            particle[pn].dNdyptdptdphi[nry-1][npt][nphi], phir);
    fy2 = lin_int(PHI[nphi-1], PHI[nphi], 
            particle[pn].dNdyptdptdphi[nry][npt][nphi-1], 
            particle[pn].dNdyptdptdphi[nry][npt][nphi], phir);

    f2 = lin_int(RY[nry-1], RY[nry], fy1, fy2, yr);


    if(ptr > PTCHANGE){
        f1 = log(f1); 
        f2 = log(f2);
    }
    /* pt interpolation */
    val = lin_int(particle[pn].pt[npt-1],particle[pn].pt[npt], f1, f2, ptr);
    if(ptr > PTCHANGE)
        val = exp(val);

    /*
       printf(" nphi  %i npt %i \n", nphi,npt);
       printf(" f1  %15.8le %15.8le  \n", f1, f2);
       printf(" phi  %15.8lf %15.8lf  \n", PHI[nphi-1], PHI[nphi]); 
       printf(" pt   %15.8lf %15.8lf  \n", particle[pn].pt[npt-1],particle[pn].pt[npt]);
       printf(" phi  %15.8lf pt %15.8lf    val %15.8lf \n", phir, ptr,val); 
       printf(" phi %15.8le %15.8le \n",particle[pn].dNdptdphi[npt][nphi-1],
       particle[pn].dNdptdphi[npt][nphi]);
       printf(" pt  %15.8le %15.8le \n",particle[pn].dNdptdphi[npt-1][nphi-1],
       particle[pn].dNdptdphi[npt-1][nphi]);

       exit(0);
     */
    return val;
}



