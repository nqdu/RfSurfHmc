#include<Eigen/Dense>
#include<iostream>
#include<string>
#include<surfdisp96.h>
using namespace Eigen;
/*
    Compute dispersion curves with fortran subroutines.
    here the subroutine is wrapped by a c++ interface

    Parameters:
        thkm,vpm,vsm,rhom   : parameters for each layer
        iflsph              : 0,flat earth,1 spherical symmetrical earth
        type                : wavetype(R,L) + phase type(c,g),string with 2 characters 
        mode                : 1 fundamental mode, 2 and higher for higher modes
        t                   : period vector
    
    Returns: 
        cg                  : dispersion velocity 

*/
VectorXd disp(VectorXf thkm,VectorXf vpm,VectorXf vsm,
        VectorXf rhom,VectorXd t,std ::string type="Rc",
        int iflsph=0,int mode=1)
{   
    int nz = thkm.size(),nt = t.size();
    VectorXd cg(nt);
    int igr,iwave;
    if(type=="Rc"){
        iwave = 2;
        igr = 0;
    }
    else if(type=="Rg"){
        iwave = 2;
        igr = 1;
    }
    else if(type=="Lc"){
        iwave = 1;
        igr = 0;
    }
    else if (type=="Lg"){
        iwave = 1;
        igr = 1;
    }
    else{
        std ::cout << "wavetp should be one of {Rc,Rg,Lc,Lg} " << std::endl;
        exit(0);
    }

    surfdisp96_(&thkm(0),&vpm(0),&vsm(0),&rhom(0),&nz,
            &iflsph,&iwave,&mode,&igr,&nt,&t(0),&cg(0));
    return cg;
}

/*
    compute vp and rho by using empirical relations
*/
void empirical_relation(float *vs,float *vp,float *rho,int nz)
{
    for(int i=0;i<nz;i++){
        vp[i] = 0.9409 + 2.0947*vs[i] - 0.8206*pow(vs[i],2) +
                0.2683*pow(vs[i],3) - 0.0251*pow(vs[i],4);
        rho[i] = 1.6612*vp[i] - 0.4721*pow(vp[i],2) + 
                0.0671*pow(vp[i],3) - 0.0043*pow(vp[i],4) 
                + 0.000106*pow(vp[i],5); 
    }
}