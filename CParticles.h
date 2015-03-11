#ifndef __CPARTICLE
#define __CPARTICLE
#include<string>
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

void readspec();
#endif
