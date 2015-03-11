#include"CHydro.h"

using namespace std;

int main(int argc, char **argv)
{


    int IEVENT=0;
    if(argc != 2){
        cout<<"useage: ./hydro IEVENT"<<endl;
    }
    else{
    IEVENT=atoi(argv[1]);
    }

    //CHydro3D hydro(-152,152, -152,152, -42,42);
    CHydro3D hydro(-81,81, -81,81, -51,51);
    //CHydro3D hydro(-71,71, -71,71, -41,41);

    char HydroSetting[50];
    char DirEvent[50];

    int I=1;
    sprintf(HydroSetting,"Vishydro%d.inp",I);
    sprintf(DirEvent,"results/event%d",IEVENT);

    hydro.IEVENT = IEVENT;

    struct stat st;
    if(!(stat(DirEvent,&st)==0)){
        //Judge if directory exist, if not, create one
        if(mkdir(DirEvent,0777)==-1){
            cout<<"mkdir error!";
            exit(0);
        }
    }//create directory results/eventi


    hydro.Initialization(HydroSetting);

    hydro.Evolution();  
    /*Do hydro evolution and frz out hypersurface calculation */
    /*Frz out hypersf is stored in results/HyperSfI.dat       */

    return 0;

}
