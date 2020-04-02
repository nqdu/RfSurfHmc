#include<celldisp.h>
#include<fstream>
#include<math.h>
#include<disp.h>
#include<gaussian.h>

/*
    Forward computation : G(m)
*/
VectorXd LayerModel :: forward()
{
    VectorXd data(nt);
    VectorXf vp(nz),rho(nz);
    int k = 0;

    // compute vp and rho according to empirical relations
    empirical_relation(&v(0),&vp(0),&rho(0),nz);

    // forward computation
    if(ntRc > 0){
        data.segment(k,ntRc)= disp(thk,vp,v,rho,tRc,"Rc",iflsph);
        k += ntRc;
    }
    if(ntRg > 0){
        data.segment(k,ntRg)= disp(thk,vp,v,rho,tRg,"Rg",iflsph);
        k += ntRg;        
    }
    if(ntLc > 0){
        data.segment(k,ntLc)= disp(thk,vp,v,rho,tLc,"Lc",iflsph);
        k += ntLc;        
    }   
    if(ntLg > 0){
        data.segment(k,ntLg)= disp(thk,vp,v,rho,tLg,"Lg",iflsph);      
    } 

    return data; 
}
/*
    compute depthkernel of this model 
*/
MatrixXd LayerModel:: gradient()
{
    MatrixXd kernel(nz,nt);
    //#pragma omp parallel num_threads(2)
    //#pragma omp parallel for 
    for(int k=0;k<nz;k++){
        VectorXd c1(nt),c2(nt);
        float dv,v0;
        dv = 0.01;

        // compute kernel
        v0 = v(k);
        v(k) = v0 * (1 - 0.5 * dv );
        c1 = forward();
        v(k) = v0 * (1 + 0.5 * dv );
        c2 = forward();
        kernel.row(k) = (c2 - c1) / (dv * v0);
        v(k) = v0;
    }

    return kernel;
}