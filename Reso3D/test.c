#include<stdio.h>
#include<stdlib.h>

int main()
{
    bool bsmall=(1.0E-3380)>(1.0E-3300) ? true : false;
    printf("bsmall=%d",bsmall);
    return 0;
}

