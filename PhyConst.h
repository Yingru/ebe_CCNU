#ifndef __PHYCONST
#define __PHYCONST
static long NSEED=234567;

namespace PhyConst{
	const double pi=3.1415926535;
	const double hbarc=0.19733;
	const double acu=1.0e-14;   //numerical precision
	enum {EOSI, EOSQ, SMEOSQ, EOSL, EOSL2};
	//note EOSI=0, EOSQ=1, SMEOSQ=2, EOSL=3
}

struct HydroConfig{               //Input parameters
       int Ieos;                  //Equation of State
       double b;                  //ImpactParameter
       double t0;                 //Initial time (unit: fm)
       double e0;                 //Initial energy density e0(b=0,r=0) (unit: Gev/fm^3)
       double edec;               //Decoupling energy density (unit: GeV/fm^3)

       int Ieos2dec;              //only for IEOS=2 (EOSI)  0: decouple by gluon 1: decouple by pions   
       int Iglaubercgc;           //initializion by Glauber or CGC  (0: Glauber, 1:CGC)      
       int Iein;                  //type of initialization  0:initialize by energy density 1: initialize by entropy density  
       double Hwn;                //% ratio of WN   (1-HWN) % of BC        

       double A;                   //mass of atom
       double Si0;                 // cross section for NN (unit: fm^2) 40mb, sqrt(SNN)=130GeV       
       double Ro0;                 // nuclear density  (unit: fm^-3)
       double Eta;                 // diffussion parameter of nuclear radii  (unit: fm)
       double Ra;                  // nuclear radii (unit: fm)
       double Eta_flat;            // eta_flat     half width of the rapidity platform in the central rapidity
       double Eta_gw;              // Gaussion width at large rapidity exp(-(abs(eta)-eta_flat)**2.0/(2.0*Gw**2.0))

       double ViscousC;            // eta/s (const.) (0: no shear viscosity)
       double VisBeta;             // \tau_Pi=VisBeta*6.0\eta /(ST)
       double VisBulk;             // VisBulk=C;  Xi/s= C* (Xi/S)_min  (C=0, no bulk vis; or C>1 ) 
       double IRelaxBulk;          // type of bulk relaxation time (0: critical slowing down; 1: contant relaxation time; 2: 1.5/(2 Pi T))
       double BulkTau;             // constant relaxation time (fm/c) (require input IRelaxBulk=1)  
       int NXD;
       int NYD;                   // freeze-out step in x and y dirction 
       int NZD;                   // freeze-out step in x and y dirction 
       int NTD;                   // freeze-out step in \tau dirction 
       double dt;                 //time step
       double dx;                 //
       double dy;
       double deta;              //rapidity direction grid size
};
#endif
