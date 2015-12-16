#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <complex>
#include "DefClass.h"
#include "feast.h"
#include "feast_sparse.h"
#include "mkl.h"

#define METHOD VSL_RNG_METHOD_UNIFORM_STD
#define I std::complex<double> (0,1)

WG::WG()
{
}

void WG::getPara(int norm,int diag,double t,double Dis,double alpha,double ep0)
{
    norm_=norm;
    diag_=diag;
    t_=t;
    Dis_=Dis;
    alpha_=alpha;
    ep0_=ep0;
    vslNewStream(&stream_,VSL_BRNG_MT19937,777);

    psi_=new std::complex<double> [N_];
//    vec_=new std::complex<double> [N_*N_];
//    eig_=new std::complex<double> [N_];
//    psi_=new double [2*N_];
    vec_=new double [2*N_*N_];
    eig_=new double [2*N_];
    IA_=new int [N_+1];

    if (diag==0)
    {
        Drand_=new double [N_-1];
//        H_=new std::complex<double> [N_-1];
        H_=new double [2*(N_-1)];
        JA_=new int [N_-1];
    }
    else
    {
        Drand_=new double [N_];
//        H_=new std::complex<double> [2*N_-1];
        H_=new double [2*(2*N_-1)];
        JA_=new int [2*N_-1];
    }
}

void WG::geneDis()
{
    int i;
    if (diag_==0)
    {
        vdRngUniform(METHOD,stream_,N_-1,Drand_,-Dis_,Dis_);
        for (i=0;i<N_-1;i++)
        {
//            *(H_+2*i)=exp(-2.0**(Drand_+i))+I*alpha_*(1.0+*(Drand_+i))*exp(-2.0**(Drand_+i));
            *(H_+2*i)=exp(-2.0**(Drand_+i));
            *(H_+2*i+1)=alpha_*(1.0+*(Drand_+i))*exp(-2.0**(Drand_+i));
            *(IA_+i)=i+1;
            *(JA_+i)=i+2;
        }
        *(IA_+N_-1)=N_-1;
        *(IA_+N_)=N_;
    }
    else
    {
        vdRngUniform(METHOD,stream_,N_,Drand_,-Dis_,Dis_);
        double tep=0.0;
        for (i=0;i<N_;i++)
        {
            tep+=*(Drand_+i);
        }
        tep=tep/N_;
        for (i=0;i<N_-1;i++)
        {
//            *(H_+2*i)=6.64*(*(Drand_+i)/tep-1.0)/(1.0-ep0_);
            *(H_+4*i)=6.64*(*(Drand_+i)/tep-1.0)/(1.0-ep0_);
            *(H_+4*i+1)=0.0;
//            *(H_+2*i+1)=((*(Drand_+i)+*(Drand_+i+1))/tep-1.0-ep0_)/(1.0-ep0_)+I*alpha_;
            *(H_+4*i+2)=((*(Drand_+i)+*(Drand_+i+1))/tep-1.0-ep0_)/(1.0-ep0_);
            *(H_+4*i+3)=alpha_;
        }
        *(H_+4*N_-2)=6.64*(*(Drand_+N_-1)-1.0)/(1.0-ep0_);
        *(H_+4*N_-3)=0.0;
    }
}

void WG::eigensolve()
{
    char UPLO='U';
    int feastparam[64],loop,M0,M,info;
    double epsout,r,*res;
//    std::complex<double> Emid=(0,0);
    double Emid[2];

    Emid[0]=0;
    Emid[1]=0;
    res=new double [N_];

    feastinit(feastparam);

    r=4.0;
    feastparam[17]=50;
    zfeast_scsrev(&UPLO,&N_,H_,IA_,JA_,feastparam,&epsout,&loop,Emid,&r,&M0,eig_,vec_,&M,res,&info);
    if (info!=0)
    {
        std::cout << "info = " << info << std::endl;
    }
}

void WG::sim()
{
    std::complex<double> psiinit[N_];
    int i,j;

    for (i=0;i<N_;i++)
    {
//        *(psiinit+i)=std::conj(*(vec_+N_*i+(N_-1)/2));
        *(psiinit+i)=std::conj(*(vec_+2*N_*i+N_-1)+I**(vec_+2*N_*i+N_));
        *(psiinit+i)=*(psiinit+i)*exp(I*(*(eig_+2*i)+I**(eig_+2*i+1))*t_);
    }
    for (i=0;i<N_;i++)
    {
        for (j=0;j<N_;j++)
        {
            *(psi_+i)+=(*(vec_+2*i*N_+2*j)+I**(vec_+2*i*N_+2*j+1))**(psiinit+j);
        }
    }
}

void WG::getI(double *In)
{
    int i;
    for (i=0;i<N_;i++)
    {
        *(In+i)=std::norm(*(psi_+i));
    }
}

void WG::disp()
{

}

EM::EM(int Idx,int Nem,int norm,int diag,double t,double Dis,double alpha,double ep0)
{
    N_=51; // 51 waveguides
    I_=Idx;
    Nem_=Nem;
    WG::getPara(norm,diag,t,Dis,alpha,ep0);
    En_=new double [Nem_*N_];

    char filename[20];
    sprintf(filename,"D%d.txt",I_);
    out_.open(filename,std::ostream::out);
}

void EM::go()
{
    int i;
    for (i=0;i<Nem_;i++)
    {
        WG::geneDis();
        WG::eigensolve();
        WG::sim();
        WG::getI(En_+i*N_);
    }
    void record();
}

void EM::record()
{
    double tE1[N_],tE2[N_];
    int i,j;

    for (i=0;i<N_;i++)
    {
        for (j=0;j<Nem_;j++)
        {
            *(tE1+i)+=*(En_+j*N_+i)**(En_+j*N_+i);
            *(tE2+i)+=*(En_+j*N_+i);
        }
        *(tE1+i)=*(tE1+i)/Nem_;
        *(tE2+i)=*(tE2+i)**(tE2+i)/(Nem_*Nem_);
    }
    for (i=0;i<N_;i++)
    {
        out_ << tE1[i]/tE2[i] << std::endl;
    }
}
