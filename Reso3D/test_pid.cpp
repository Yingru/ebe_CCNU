#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cassert>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>

using namespace std;

int main(int argc, char** argv)
{


  if(mkdir("mkdir_test",0777) ==-1)cout<<"mkdir erro!"<<endl; 
  else cout<<"mkdir succeed!"<<endl;

return 0;
}
