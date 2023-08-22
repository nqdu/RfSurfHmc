extern "C"{

// base on nanqiao's surfdisp.hpp
// That's the cpp warpper , packing the fortran code

/**
 * computer the rf 
**/

void cal_rf_par_time_(double *thk,double *vp,double *vs,double *rho,
                       double *qa,double *qb,
                      int nlayer,int nt,double dt,double ray_p,
                      double gauss,double time_shift,
                      int rf_type,int par_type,double *rcv_fun,
                      double *rcv_fun_p);

void cal_rf_time_(double *thk,double *vp,double *vs,double *rho,
                    double *qa,double *qb,
                  int nlayer,int nt,double dt,double ray_p,
                  double gauss,double time_shift,
                   int rf_type,double *rcv_fun);


void cal_rf_par_freq_(double *thk,double *vp,double *vs,double *rho,
                       double *qa,double *qb,
                      int nlayer,int nt,double dt,double ray_p,
                      double gauss,double time_shift,double water,
                      int rf_type,int par_type,double *rcv_fun,
                      double *rcv_fun_p);

void cal_rf_freq_(double *thk,double *vp,double *vs,double *rho,
                    double *qa,double *qb,
                  int nlayer,int nt,double dt,double ray_p,
                  double gauss,double time_shift,double water,
                   int rf_type,double *rcv_fun);

}

