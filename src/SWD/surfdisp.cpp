#include"surfdisp.hpp"
#include<math.h>

/**
 * convert phase/group velocity of flat earth to that of spherical earth. 
 * Refer to Schwab, F. A., and L. Knopoff (1972). Fast surface wave and free
 * mode computations, in  Methods in Computational Physics, Volume 11,
 * Love Wave Equations  44, 45 , 41 pp 112-113. 
 * Rayleigh Wave Equations 102, 108, 109 pp 142, 144.
 * @param t period
 * @param c phase/group velo
 * @param wavetp wavetype one of [Rc,Rg,Lc,Lg]
 * 
 * @return cnew phase/group velocity of spherical earth
 */
double _flat2sphere(double t,double c,std::string wavetp)
{
    double ar = 6371.0; // earth radius
    double tm;
    double omega = 2.0 * M_PI / t; // angular frequency

    // check input parameters
    bool flag = wavetp == "Rc" || wavetp == "Rg" || 
                wavetp == "Lc" || wavetp == "Lg";
    if(flag == false){
        std::cout << "cnew = _flat2sphere(double t,double c,std::string wavetp)\n";
        std::cout << "parameter wavetp should be in one of [Rc,Rg,Lc,Lg]\n ";
        exit(0);
    }

    if(wavetp[0] == 'L'){ //love wave
        tm = 1. + pow(1.5 * c / (ar * omega),2); 
    }
    else{ // Rayleigh wave
        tm = 1. + pow(0.5 * c / (ar * omega),2);
    }
    tm = sqrt(tm);

    // convert to spherical velocity
    double cnew;
    if(wavetp[1] == 'c'){ // phase velocity
        cnew = c / tm;
    }
    else{ // group velocity
        cnew = c * tm;
    }

    return cnew;
}

/**
 * calculates the dispersion values for any layered model, any frequency, and any mode.
 * @param nlayer no. of layers
 * @param thkm,vpm,vsm,rhom model, shape(nlayer)
 * @param kmax no. of periods used
 * @param t,cp period and phase velocity, shape(kmax)
 * @param sphere true for spherical earth, false for flat earth 
 * @param wavetype one of [Rc,Rg,Lc,Lg]
 * @param mode i-th mode of surface wave, 0 fundamental, 1 first higher, ....
 * @param keep_flat keep flat earth phase/group velocity or convert it to spherical
 */
int _surfdisp(float *thk,float *vp,float *vs,float *rho,
            int nlayer,double *t,double *cg,int kmax,std::string wavetype,
            int mode,bool sphere,bool keep_flat)
{
    int iwave,igr,ifsph=0;
    if(wavetype=="Rc"){
        iwave = 2;
        igr = 0;
    }
    else if(wavetype == "Rg"){
        iwave = 2;
        igr = 1;
    }
    else if(wavetype=="Lc"){
        iwave = 1;
        igr = 0;
    }
    else if(wavetype=="Lg"){
        iwave = 1;
        igr = 1;
    }
    else{
        std::cout <<"wavetype should be one of [Rc,Rg,Lc,Lg]"<<std::endl;
        exit(0);
    }

    if(sphere == true) ifsph = 1;
    int ierr;
    surfdisp96_(thk,vp,vs,rho,nlayer,ifsph,iwave,mode+1,igr,kmax,t,cg,&ierr);

    if(ierr == 1){
        return ierr;
    } 
    if(sphere == true && keep_flat == false){
        for(int i=0;i<kmax;i++){
            cg[i] = _flat2sphere(t[i],cg[i],wavetype);
        }
    }

    return ierr;
}

/** Love group velocity with analytical method  
 * @param nlayer no. of layers
 * @param thk,vs,rho 1-D model
 * @param kmax no. of periods
 * @param t cg periods and group velocity
 * @param mode i-th mode of surface wave, 1 fundamental, 2 first higher, ....
 * @param sphere true for spherical earth, false for flat earth 
 */
int _LoveGroup(float *thk,float *vs,float *rho,
              int nlayer,double *t,double *cg,int kmax,
              int mode,bool sphere)
{
    // check if earth model is spherical
    int iflsph = 0; if(sphere == true) iflsph = 1;

    // temporay arrays
    float vp[kmax];
    double cp[kmax],uu[nlayer],tt[nlayer],dcdh[nlayer],
            dcdr[nlayer],dcdb[nlayer];

    // compute phase velocity
    for(int i=0;i<kmax;i++) vp[i] = 1.732 * vs[i];
    int ierr = _surfdisp(thk,vp,vs,rho,nlayer,t,cp,kmax,"Lc",mode,sphere,true);
    if(ierr == 1) return ierr;

    // convert phase to group velocity
    for(int i=0;i<kmax;i++){
        slegn96_(thk,vs,rho,nlayer,t+i,cp+i,cg+i,uu,tt,dcdb,dcdh,dcdr,iflsph);
    }
    return ierr;
}

/** Rayleigh group velocity with analytical method  
 * @param nlayer no. of layers
 * @param thk,vp,vs,rho 1-D model
 * @param kmax no. of periods
 * @param t cg periods and group velocity
 * @param mode i-th mode of surface wave, 1 fundamental, 2 first higher, ....
 * @param sphere true for spherical earth, false for flat earth 
 */
int _RayleighGroup(float *thk,float *vp,float *vs,float *rho,
              int nlayer,double *t,double *cg,int kmax,
              int mode,bool sphere)
{
    // check if earth model is spherical
    int iflsph = 0; if(sphere == true) iflsph = 1;

    // temporay arrays
    double cp[kmax],ur[nlayer],uz[nlayer],tr[nlayer],tz[nlayer],
            dcdh[nlayer],dcda[nlayer],dcdr[nlayer],dcdb[nlayer];

    // compute phase velocity
    int ierr = _surfdisp(thk,vp,vs,rho,nlayer,t,cp,kmax,"Rc",mode,sphere,true);
    if(ierr == 1) return ierr;

    // convert phase to group velocity
    for(int i=0;i<kmax;i++){
        sregn96_(thk,vp,vs,rho,nlayer,t+i,cp+i,cg+i,
                ur,uz,tr,tz,dcda,dcdb,dcdh,dcdr,iflsph);
    }

    return ierr;
}


/**
 * compute phase velocity sensitivity kernel, and group velocity 
 * for layer-based model. note that if in spherical model, 
 * input cp should be cp_flat
 * @param nlayer no. of layers
 * @param thk,vp,vs,rhom model, shape(nlayer)
 * @param nt no.of periods used
 * @param t,cp,cg period and phase/group velocity
 * @param dc2da(b,h,r) phase velocity sensitivity kernel for vp,vs,thick and rho shape(nlayer,nt)
 * @param du2da(b,h,r) phase velocity sensitivity kernel for vp,vs,thick and rho shape(nlayer,nt)
 * @param sphere true for spherical earth, false for flat earth 
 * @param wavetype one of [Rc,Rg,Lc,Lg]
 * @param mode i-th mode of surface wave, 0 fundamental, 1 first higher, ....
 */
int _SurfKernel(float *thk,float *vp,float *vs,float *rho,
                int nlayer,double *t,double *c,int nt,double *dcda,
                double *dcdb,double *dcdr,double *dcdh,
                std::string wavetp,int mode,bool sphere=true)
{
    // check input parameters
    bool flag = wavetp == "Rc" || wavetp == "Rg" || 
                wavetp == "Lc" || wavetp == "Lg";
    if(flag == false){
        std::cout << "cnew = _flat2sphere(double t,double c,std::string wavetp)\n";
        std::cout << "parameter wavetp should be in one of [Rc,Rg,Lc,Lg]\n ";
        exit(0);
    }

    // compute velocity and kernel
    bool keep_flat = true;
    int iflsph = 0; if(sphere) iflsph = 1;

    int ierr;
    if(wavetp == "Rc"){ // rayleigh phase
        // compute cf
        //printf("enterin");
       ierr =  _surfdisp(thk,vp,vs,rho,nlayer,t,c,nt,wavetp,mode,sphere,keep_flat);
       if(ierr == 1) {
           for(int i=0;i<nlayer;i++){
               std::cout << vs[i] << std::endl;
           }
           return ierr;
       }
       //printf("entermd");

        // compute kernel
        double cg,ur[nlayer],uz[nlayer],tr[nlayer],tz[nlayer];
        for(int i=0;i<nt;i++){
            int k = i * nlayer;
            sregn96_(thk,vp,vs,rho,nlayer,t+i,c+i,&cg,ur,uz,tr,tz,dcda+k,
                    dcdb+k,dcdh+k,dcdr+k,iflsph);
        }
        //printf("enterout\n");
        
    }
    else if(wavetp == "Rg"){
        
        // compute cf
        double cp[nt],cp1[nt],cp2[nt],t1[nt],t2[nt];
        for(int i=0;i<nt;i++){
            t1[i] = t[i] * (1.0 + 0.05);
            t2[i] = t[i] * (1.0 - 0.05);
        }
        ierr = _surfdisp(thk,vp,vs,rho,nlayer,t,cp,nt,"Rc",mode,sphere,keep_flat);
        int ierr1 = _surfdisp(thk,vp,vs,rho,nlayer,t1,cp1,nt,"Rc",mode,sphere,keep_flat);
        int ierr2 = _surfdisp(thk,vp,vs,rho,nlayer,t2,cp2,nt,"Rc",mode,sphere,keep_flat);
        ierr = ierr + ierr1 + ierr2 > 0;
        //printf("rayleigh group %d\n",ierr);
        if(ierr == 1) return ierr;
        

        // compute kernel
        double ur[nlayer],uz[nlayer],tr[nlayer],tz[nlayer];
        double dcdb1[nlayer*nt],dcda1[nlayer*nt],
                dcdr1[nlayer*nt],dcdh1[nlayer*nt];
        for(int i=0;i<nt;i++){
            int k = i * nlayer;
            sregnpu_(thk,vp,vs,rho,nlayer,t+i,cp+i,c+i,ur,uz,tr,tz,
                    t1+i,cp1+i,t2+i,cp2+i,dcda1+k,dcdb1+k,dcdh1+k,
                    dcdr1+k,dcda+k,dcdb+k,dcdh+k,dcdr+k,iflsph);
        }
    }
    else if(wavetp == "Lc"){ // Love phase
        // compute cf
        ierr = _surfdisp(thk,vp,vs,rho,nlayer,t,c,nt,wavetp,mode,sphere,keep_flat);
        
        if(ierr == 1) return ierr;

        // compute kernel
        double cg,uu[nlayer],tt[nlayer];
        for(int i=0;i<nt;i++){
            int k = i * nlayer;
            slegn96_(thk,vs,rho,nlayer,t+i,c+i,&cg,
                    uu,tt,dcdb+k,dcdh+k,dcdr+k,iflsph);
        }
    }
    else{
        // compute cf
        double cp[nt],cp1[nt],cp2[nt],t1[nt],t2[nt];
        for(int i=0;i<nt;i++){
            t1[i] = t[i] * (1.0 + 0.05);
            t2[i] = t[i] * (1.0 - 0.05);
        }
        ierr = _surfdisp(thk,vp,vs,rho,nlayer,t,cp,nt,"Lc",mode,sphere,keep_flat);
        int ierr1 = _surfdisp(thk,vp,vs,rho,nlayer,t1,cp1,nt,"Lc",mode,sphere,keep_flat);
        int ierr2 = _surfdisp(thk,vp,vs,rho,nlayer,t2,cp2,nt,"Lc",mode,sphere,keep_flat);
        ierr = ierr || ierr1 || ierr2;
        if(ierr == 1) return ierr;

        // compute kernel
        double uu[nlayer],tt[nlayer];
        double dcdb1[nlayer*nt],dcdr1[nlayer*nt],dcdh1[nlayer*nt];
        for(int i=0;i<nt;i++){
            int k = i * nlayer;
            slegnpu_(thk,vs,rho,nlayer,t+i,cp+i,c+i,uu,tt,t1+i,cp1+i,
                    t2+i,cp2+i,dcdb1+k,dcdh1+k,dcdr1+k,dcdb+k,dcdh+k,
                    dcdr+k,iflsph);
        }
    }
    
    return ierr;
}