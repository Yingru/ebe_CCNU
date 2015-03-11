#include"CParticles.h"
#include<iostream>
#include<fstream>
#include<vector>

using namespace std;

vector<CParticle> particles;

void readspec(){
    double buff;
    char oneline[256];
    ifstream fin("Reso3D/resoweak.dat");
    CParticle p;
    while(!fin.eof()){
        fin>>p.monval>>p.name>>p.mass>>p.width         \
            >>p.gspin>>p.baryon>>p.strange>>p.charm     \
            >>p.bottom>>p.gisospin>>p.charge>>p.decays;

        fin>>buff>>buff>>buff; //Not used data
        fin.getline(oneline,256);//Skip the current line
        for(int k=0; k<p.decays; k++)
            fin.getline(oneline,256);//Skip the daughters' data
        particles.push_back(p);
    }
    particles.pop_back();  //The last line of the file is empty
    fin.close();
}

//#define __TESTPARTICLE
#ifdef __TESTPARTICLE

int main()
{
    readspec();
    //  for(vector<CParticle>::iterator it=particles.begin(); it !=particles.end(); ++it)cout<< (*it).name<<endl;
    vector<CParticle>::size_type i;
    for(i=0; i!=particles.size(); i++)
        cout<<particles[i].name<<endl;
}
#endif
