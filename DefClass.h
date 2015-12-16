#ifndef DEFCLASS_H
#define DEFCLASS_H
#include <fstream>
#include <complex>
#include "mkl.h"
class WG
{
    public:
        WG();
    protected:
        std::complex<double> *psi_;
        double *H_,*vec_,*eig_;
        double Dis_,alpha_,t_,ep0_,*Drand_;
        int N_,norm_,diag_,*IA_,*JA_;
        VSLStreamStatePtr stream_;

        void getPara(int,int,double,double,double,double);
        void geneDis();
        void eigensolve();
        void sim();
        void disp();
        void getI(double *);
};

class EM : protected WG
{
    public:
        EM(int,int,int,int,double,double,double,double);
        void go();
    private:
        double *En_;
        int I_,Nem_;
        std::ofstream out_;
        void record();
};
#endif
