#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include<disp.h>
using Eigen :: VectorXf;
using Eigen :: VectorXd;
using std:: cout;
using std:: endl;
using std:: ifstream;
using std:: ofstream;

int main(){
    // read velocity model
    int nz;
    ifstream infile;
    ofstream outfile;
    infile.open("input/true.dat");
    infile >> nz;
    VectorXf thk(nz),vp(nz),vs(nz),rho(nz);
    for(int i=0;i<nz;i++){
        infile >> thk(i) >> vs(i);
    }
    std::cout << thk<<std::endl;
    empirical_relation(&vs(0),&vp(0),&rho(0),nz);

    infile.close();

    // compute dispersion curve
    int nt = 20;
    VectorXd T(nt);
    for(int i=0;i<nt;i++){
        T(i) = 4.0 + i * 2.0;
    }
    VectorXd cg(nt);
    cg = disp(thk,vp,vs,rho,T,"Rc");
    outfile.open("Rc.true");
    for(int i=0;i<nt;i++){
        outfile << T(i) << " " << cg(i) << endl;
    }
    outfile.close();

    cg = disp(thk,vp,vs,rho,T,"Rg");
    outfile.open("Rg.true");
    for(int i=0;i<nt;i++){
        outfile << T(i) << " " << cg(i) << endl;
    }
    outfile.close();

    /*
    cg = disp(thk,vp,vs,rho,T,"Lc");
    outfile.open("Lc.true");
    for(int i=0;i<nt;i++){
        outfile << T(i) << " " << cg(i) << endl;
    }
    outfile.close();

    cg = disp(thk,vp,vs,rho,T,"Lg");
    outfile.open("Lg.true");
    for(int i=0;i<nt;i++){
        outfile << T(i) << " " << cg(i) << endl;
    }
    outfile.close();
    */
    return 0;
    
}