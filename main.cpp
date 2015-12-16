#include "DefClass.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>
#include "mkl.h"

int main (int argc, char *argv[])
{
    int Idx,Nem,norm,diag;
    double t,Dis,alpha,ep0;

    if (argc==0)
    {
        Idx=(int)atof(argv[1]);
        Nem=(int)atof(argv[2]);
        norm=(int)atof(argv[3]);
        diag=(int)atof(argv[4]);
        t=(double)atof(argv[5]);
        Dis=(double)atof(argv[6]);
        alpha=(double)atof(argv[7]);
        ep0=(double)atof(argv[8]);

        EM one(Idx,Nem,norm,diag,t,Dis,alpha,ep0);
        one.go();
    }

    return 0;
}
