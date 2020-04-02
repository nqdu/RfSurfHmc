#pragma once
#include<Eigen/Dense>
#include<iostream>
#include<string>
using Eigen :: VectorXf;
using Eigen :: VectorXd;

VectorXd disp(VectorXf thkm,VectorXf vpm,VectorXf vsm,
        VectorXf rhom,VectorXd t,std ::string wavetp="Rc",
        int iflsph=0,int mode=1);
void empirical_relation(float *vs,float *vp,float *rho,int nz);