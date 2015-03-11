#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cstdio>

using namespace std;

int main()
{
    const int NX=141;
    const int NY=141;
    const int NEVENT=50;
    float ED_average[NX][NY]; //={0.0};
    for(int i=0; i<NX; i++)
	for(int j=0; j<NY; j++){
	    ED_average[i][j] = 0.0;
	}

    char File_ED[50];
    for(int NT=50; NT<=NEVENT+50; NT++){
	sprintf(File_ED,"results/ED_XY%d.dat",NT); 
	ifstream fin(File_ED);
	float ed;
	for(int i=0; i<NX; i++)
	    for(int j=0; j<NY; j++){
		fin>>ed;
		ED_average[i][j] += ed/NEVENT;
	    }
	fin.close();
    }
    ofstream fout("ED_Average.dat");
    for(int i=0; i<NX; i++){
	for(int j=0; j<NY; j++){
	    fout<<ED_average[i][j]<<" ";
	}
	cout<<ED_average[i][NY/2]<<endl;
    }

    fout.close();

}
