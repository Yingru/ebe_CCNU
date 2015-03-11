/*    reso.c   */
/* 
 ********************************************************************************
 This is the public version of the resonance decay calculation using the output 
 generated by the hydrodynamical code azhydro0p2.  The majority of the code was 
 developed by Josef Sollfrank and Peter Kolb.  Additional work and comments 
 were added by Evan Frodermann, September 2005.
 Please refer to the papers 
 J. Sollfrank, P. Koch, and U. Heinz, Phys. Lett B 252 (1990) and 
 J. Sollfrank, P. Koch, and U. Heinz, Z. Phys. C 52 (1991) 
 for a description of the formalism utilized in this program.
 ********************************************************************************
 */


#include	<string.h>
#include	<stdio.h>
#include	<math.h>

#include	"reso.h"
#include	"functions.h"

//*****************************************************
//This is the main program that calls all the needed subroutines.
main(int argc, char**argv) {
    //    assert(argc==2);

    IEVENT = atoi(argv[1]);

    FILE	*datin;

    char	infile[FILEDIM];
    char	outdir[FILEDIM];
    char    specfile[FILEDIM];
    char	dummy[200];

    int max, maxdecay, bound;

    printf("START of reso decays !\n");
    //Read in the data from "reso.inp, including the results folder and the spectra data      
    datin = fopen("reso.inp","r");
    fscanf(datin,"%s%s",infile,dummy);
    fscanf(datin,"%s%s",specfile,dummy);
    fscanf(datin,"%s%s",outdir,dummy);
    fscanf(datin,"%i%s",&bound,dummy);
    fclose(datin);
    //Read in the spectra and decays using "resoweak.dat" as a database of particles
    readspec(infile,specfile, &max, &maxdecay);
    //The main module that calculates the resonance decay feeddown
    cal_reso_decays(max, maxdecay, bound); 
    //Writes the spectra to specified data files.
    writespec(max,outdir);


    /* Test dNdp3(y,pt,phi,reso_num)   */
/*    int pn=max-1;
    int part = particle[pn].monval;
    int npt, nphi, nry;

    npt = particle[pn].npt;
    nphi = particle[pn].nphi;
    nry = particle[pn].nry;
    double Y[nry];

    int h, l, i;
    double ry_max=5.0;
    double ry_min=-5.0;
    double dry=(ry_max-ry_min)/(nry-1.0);
    for(i=0; i<nry; i++){
        Y[i] = ry_min+i*dry;
    }

    for (h = 0; h < nry; h++){
        double y = Y[h];                         
        for (l = 0; l < npt; l++)
        {
            for (i = 0; i < nphi; i++)
            {

            printf("%11.4le ",Edndp3 (y, particle[pn].pt[l], PHI[i], part));

            }
            printf("\n");
        }
            printf("\n\n");
    }

*/
//    int part = 211;
//    double dn=Edndp3(-6.0, 4.0, 1.57, 211);
//    printf("%11.4le ",dn);
    /////////////////////////////////////
    printf("END of reso decays ! \n");

    return 0;	/* ok */
}





