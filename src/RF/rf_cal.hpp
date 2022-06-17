#pragma once
extern "C"{

// base on nanqiao's surfdisp.hpp
// That's the cpp warpper , packing the fortran code

/**
 * computer the rf 
**/

void receiverFunc_par_(double *thk,double *vp,double *vs,double *rho,
                      int nlayer,int nt,double dt,double ray_p,
                      double gauss,double time_shift,double water_level,
                      int rf_type,int par_type,double *rcvfun,double *rcv_fun_p);

void receiverFunc_(double *thk,double *vp,double *vs,double *rho,
                  int nlayer,int nt,double dt,double ray_p,
                  double gauss,double time_shift,double water_level,
                   int rf_type,double *rcv_fun);

}

