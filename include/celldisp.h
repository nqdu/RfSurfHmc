#pragma once
#include<Eigen/Core>
#include<unsupported/Eigen/CXX11/Tensor>
#include<iostream>
#include<string>
using Eigen :: VectorXf;
using Eigen :: VectorXd;
using Eigen :: MatrixXd;

/*
    c++ class for 1-D layer model
*/
class LayerModel{
    public:
    int nz; // no. of layers for this model
    int iflsph;  // 0 for flat earth and 1 for spherical earth 
    int ntRc,ntRg,ntLc,ntLg; // no. of period points for Rayleigh/Love phase/group velocity
    int nt;     // no. of all period points

    VectorXf thk,v; // parameters for each layer
    VectorXd tRc,tRg,tLc,tLg; // period vector
    VectorXd d;               // data vector 

    int read_inputfile(std :: string filename);
    VectorXd forward();
    MatrixXd gradient();
};