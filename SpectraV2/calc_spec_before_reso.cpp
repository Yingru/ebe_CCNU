#include"CSpectra.h"

int main(int argc, char** argv)
{

    int IEVENT=0;
    if(argc != 2){
        cout<<"useage: ./calc_spec_before_reso IEVENT"<<endl;
        exit(0);
    }
    else{
        IEVENT = atoi(argv[1]);
    }

    CSpectra spec(1);
    //spec.ReadParticles("../Reso3D/resoweak.dat");
    char particles_data_table[]="../s95p-v1/pdg05.dat";
    spec.ReadParticles(particles_data_table);

    for(deque<CParticle>::size_type i=0; i!=spec.particles.size();++i){
        int m=spec.particles[i].monval;
        cout<<spec.particles[i].monval<<' '<<endl;
    }

    spec.Set_FrzT(0.136);
    spec.ReadHyperSF(IEVENT);
    spec.CalcSpectra_BeforeReso(IEVENT);
    return 0;
}

