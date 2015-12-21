#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include "feast.h"
#include "feast_sparse.h"
#include "feast_banded.h"
#include "mkl.h"
#include "mkl_vsl.h"

#define N_ 51
#define LDA_ 2
#define METHOD VSL_RNG_METHOD_UNIFORM_STD

double _Complex *psi_;
double *H_,*vec_,*eig_,*Drand_,*En_;
double Dis_,alpha_,t_,ep0_;
int norm_,diag_,*IA_,*JA_,I_,Nem_,iNem_;
VSLStreamStatePtr stream_;
FILE *out_;

int main (int argc, char *argv[])
{
    int Idx,Nem,norm,diag,i;
    double t,Dis,alpha,ep0;

    if (argc!=0)
    {
        Idx=(int)atof(argv[1]);
        Nem=(int)atof(argv[2]);
        norm=(int)atof(argv[3]);
        diag=(int)atof(argv[4]);
        t=(double)atof(argv[5]);
        Dis=(double)atof(argv[6]);
        alpha=(double)atof(argv[7]);
        ep0=(double)atof(argv[8]);
    }

    void WG_getPara(int,int,int,int,double,double,double,double);
    void WG_getDis();
    void WG_eigensolve();
    void WG_sim();
    void WG_getI(double *);
    void EM_record();
    WG_getPara(Idx,Nem,norm,diag,t,Dis,alpha,ep0);

    for (i=0;i<Nem_;i++)
    {
        iNem_=i;
        WG_getDis();
        WG_eigensolve();
        WG_sim();
        WG_getI(En_+i*N_);
    }

    EM_record();

    return 0;
}

void WG_getPara(int Idx,int Nem,int norm,int diag,double t,double Dis,double alpha,double ep0)
{
    I_=Idx;
    Nem_=Nem;
    norm_=norm;
    diag_=diag;
    t_=t;
    Dis_=Dis;
    if (Dis<1e-6)
    {
        Dis_=1e-8;
    }
    alpha_=alpha;
    ep0_=ep0;
    vslNewStream(&stream_,VSL_BRNG_MT19937,time(NULL));

    psi_=(double _Complex *)malloc(N_*sizeof(double _Complex));
    vec_=(double *)malloc(2*N_*N_*sizeof(double));
    eig_=(double *)malloc(2*N_*sizeof(double));
    IA_=(int *)malloc((N_+1)*sizeof(int));
    En_=(double *)malloc(Nem_*N_*sizeof(double));

    if (diag==0)
    {
        Drand_=(double *)malloc(Nem_*(N_-1)*sizeof(double));
        H_=(double *)malloc(2*LDA_*N_*sizeof(double));
        JA_=(int *)malloc((N_-1)*sizeof(int));
        vdRngUniform(METHOD,stream_,Nem_*(N_-1),Drand_,-Dis_,Dis_);
    }
    else
    {
        Drand_=(double *)malloc(Nem_*N_*sizeof(double));
        H_=(double *)malloc(2*(2*N_-1)*sizeof(double));
        JA_=(int *)malloc((2*N_-1)*sizeof(int));
        vdRngUniform(METHOD,stream_,Nem_*N_,Drand_,-Dis_+1.0,Dis_+1.0);
    }

    char filename[20];
    sprintf(filename,"D%d.txt",I_);
    if ((out_=fopen(filename,"w"))==NULL)
    {
        printf("file open error\n");
    }
    vslDeleteStream(&stream_);
    return;
}

void WG_getDis()
{
    int i;
    if (diag_==0)
    {
        *H_=0.0;
        *(H_+1)=0.0;
        *(H_+2)=0.0;
        *(H_+3)=0.0;
        for (i=0;i<N_-1;i++)
        {
            *(H_+4*i+3+1)=exp(-2.0**(Drand_+(N_-1)*iNem_+i));
            *(H_+4*i+3+2)=alpha_*(1.0+*(Drand_+(N_-1)*iNem_+i))*exp(-2.0**(Drand_+(N_-1)*iNem_+i));
            *(H_+4*i+3+3)=0.0;
            *(H_+4*i+3+4)=0.0;
        }
    }
    else
    {
        double tep=0.0;
        for (i=0;i<N_;i++)
        {
            tep+=*(Drand_+N_*iNem_+i);
        }
        tep=tep/N_;
        for (i=0;i<N_-1;i++)
        {
            *(H_+4*i)=6.64*(*(Drand_+N_*iNem_+i)/tep-1.0)/(1.0-ep0_);
            *(H_+4*i+1)=0.0;
            *(H_+4*i+2)=((*(Drand_+N_*iNem_+i)+*(Drand_+N_*iNem_+i+1))/tep-1.0-ep0_)/(1.0-ep0_);
            *(H_+4*i+3)=alpha_;
            *(IA_+i)=2*i+1;
            *(JA_+2*i)=i+1;
            *(JA_+2*i+1)=i+2;
        }
        *(H_+4*N_-4)=6.64*(*(Drand_+N_*iNem_+N_-1)/tep-1.0)/(1.0-ep0_);
        *(H_+4*N_-3)=0.0;
        *(IA_+N_-1)=2*N_-1;
        *(IA_+N_)=2*N_;
        *(JA_+2*N_-2)=N_;
    }
}

void WG_eigensolve()
{
    char UPLO='U';
    int feastparam[64],loop,M0,M,info,N,klu,LDA;
    double epsout,r,*res;
    double Emid[2];

    N=N_;
    M0=N_;
    Emid[0]=0;
    Emid[1]=0;
    res=(double *)malloc(N_*sizeof(double));
    klu=1;
    LDA=2;

    feastinit(feastparam);

    r=4.0;
    feastparam[17]=50;
    feastparam[9]=1;
    zfeast_sbev(&UPLO,&N,&klu,H_,&LDA,feastparam,&epsout,&loop,Emid,&r,&M0,eig_,vec_,&M,res,&info);
    if (info!=0)
    {
        printf("info = %d \n",info);
    }
    free(res);
}

void WG_sim() {
    double _Complex psiinit[N_];
    int i,j;

    for (i=0;i<N_;i++)
    {
        *(psiinit+i)=conj(*(vec_+2*N_*i+N_-1)+_Complex_I**(vec_+2*N_*i+N_));
        *(psiinit+i)=*(psiinit+i)*cexp(I*(*(eig_+2*i)+I**(eig_+2*i+1))*t_);
    }
    for (i=0;i<N_;i++)
    {
        *(psi_+i)=0.0;
        for (j=0;j<N_;j++)
        {
            *(psi_+i)=*(psi_+i)+(*(vec_+2*j*N_+2*i)+I**(vec_+2*j*N_+2*i+1))**(psiinit+j);
        }
    }
}

void WG_getI(double *In)
{
    int i;
    for (i=0;i<N_;i++)
    {
        *(In+i)=pow(creal(*(psi_+i)),2.0)+pow(cimag(*(psi_+i)),2.0);
    }
}

void EM_record()
{
    double tE1[N_],tE2[N_];
    int i,j;

    for (i=0;i<N_;i++)
    {
        *(tE1+i)=0.0;
        *(tE2+i)=0.0;
        for (j=0;j<Nem_;j++)
        {
            *(tE1+i)+=pow(*(En_+j*N_+i),2.0);
            *(tE2+i)+=*(En_+j*N_+i);
        }
        *(tE1+i)=*(tE1+i)/Nem_;
        *(tE2+i)=*(tE2+i)**(tE2+i)/(Nem_*Nem_);
    }
    for (i=0;i<N_;i++)
    {
//        fprintf(out_,"%f\t%f\t%f\n",tE1[i],tE2[i],En_[i]);
        fprintf(out_,"%f\n",tE1[i]/tE2[i]);
//        fprintf(out_,"%f\t%f\t%f\n",*(eig_+2*i),*(eig_+2*i+1),En_[i]);
    }
//    for (i=0;i<Nem_*N_;i++)
//    {
//        fprintf(out_,"%f\n",*(En_+i));
//    }
}
