#pragma once 
#include<Eigen/Core>
#include<iostream>
#include<string>
#include<fstream> 
#include<math.h>
#include"gaussian.h" 
using namespace Eigen; 
template<class T> 
class HamiltonMC{
        private: 
    int n; // model dimension 
    int nd; // data dimension 
    Matrix<float,-1,2> boundary; // model boundaries
    Vector2d eps; // time interval 
    Vector2i iteration; // iteration
    MatrixXd S; // Tickonov regularization matrix  
    int niter; 
    double potential(T &model); 
    double kinetic(VectorXd p); 
    VectorXd gradient(T &model); 
    void save_disp(T &model,std::string outdir,VectorXd data,int istep=0); 

    public: 
    int read_input(T &model,std::string inputfile); 
    int run(T &model,std::string outdir); 
}; 
template<class T> double HamiltonMC<T>:: potential(T &model)
{ 
    VectorXd data(nd),tmp(nd),tmp1(n-2),x(n); 
    double s; 
    data = model.forward(); 
    tmp = data - model.d;
    x = model.v.template cast<double>();
    tmp1 = S * x;
    s = 0.5 * tmp.transpose() * tmp ; 
    s = s  + 0.01 * tmp1.transpose() * tmp1;
    return s; 
} 

template<class T> double HamiltonMC<T>:: kinetic(VectorXd p) 
{ 
    double s; 
    s = 0.5 * p.transpose() * p;
    return s; 
} 

template<class T> VectorXd HamiltonMC<T>:: gradient(T &model) 
{ 
    MatrixXd kernel(n,nd); 
    VectorXd data(nd),s(n),s1(n);

    data = model.forward(); 
    kernel = model.gradient();
    s1 = 0.01 * 2 * S.transpose() * (S * model.v.template cast<double>()); 
    // compute du_over_dm 
    s = kernel * (data - model.d) + s1; 
    return s; 
} 
// model storage information of your inversion problem 
template<class T> int HamiltonMC<T>:: read_input(T &model,std::string inputfile) 
{ 
    std :: ifstream infile; infile.open(inputfile); // read init model 
    infile >> model.nz >> model.iflsph; 
    n = model.nz; 
    if(model.iflsph){ 
        std::cout << "Spherical Symmetric Earth Model\nNo. of layers = " << n << std::endl; 
    }
    else{ 
        std::cout << "Flat Earth Model\nNo. of layers = " << n << std::endl; 
    } 
    model.thk.resize(n); 
    model.v.resize(n);
    boundary.resize(n,2); 
    for(int i = 0; i < n;i++){ 
        infile >> model.thk(i) >> boundary(i,0)>>boundary(i,1);
	//boundary(i,0) = 2.0; boundary(i,1) = 5.0; 
        float v0 = uniform(boundary(i,0),boundary(i,1)); 
        model.v(i) = v0; 
        // model.v(i) = 0.5 * ( boundary(i,0) + boundary(i,1)); 
    } 
    std::cout<<"thickness and S-wave velocity of each layer :" <<std::endl; 
    for(int i=0;i<n;i++){ 
        if(i!=n-1){
            std::cout << model.thk(i) <<" " << model.v(i)<<std::endl; 
        } 
        else{ 
            std::cout<<"half space "<< model.v(i) <<std::endl;
        }
    }

    // read dispersion curve
    int k = 0;
    infile >> model.ntRc >> model.ntRg >> model.ntLc >> model.ntLg;
    int ntRc = model.ntRc,ntRg= model.ntRg,ntLc=model.ntLc,ntLg = model.ntLg;
    model.nt = model.ntRc + model.ntRg + model.ntLc + model.ntLg;  
    nd = model.nt;
    model.d.resize(nd);
    if(model.ntRc > 0){
        model.tRc.resize(ntRc);
        for(int i=0;i<model.ntRc;i++) infile >> model.tRc(i) >> model.d(i+k);
        std::cout<<"Rayleigh wave phase period points(s):" << std::endl;
        for(int i=0;i<ntRc;i++) std::cout << model.tRc(i) << " ";
        std::cout<<std::endl;
    }
    k += model.ntRc;
    if(ntRg > 0){
        model.tRg.resize(ntRg);
        for(int i=0;i<ntRg;i++) infile >> model.tRg(i)>> model.d(i+k);
        std::cout<<"Rayleigh wave group period points(s):" << std::endl;
        for(int i=0;i<ntRg;i++) std::cout << model.tRg(i) << " ";
        std::cout<<std::endl;
    }
    k += ntRg;
    if(ntLc > 0){
        model.tLc.resize(ntLc);
        for(int i=0;i<ntLc;i++) infile >> model.tLc(i)>> model.d(i+k);
        std::cout<<"Love wave phase period points(s):" << std::endl;
        for(int i=0;i<ntLc;i++) std::cout << model.tLc(i) << " ";
        std::cout<<std::endl;
    }
    k += ntLc;
    if(ntLg > 0){
        model.tLg.resize(ntLg);
        for(int i=0;i<ntLg;i++) infile >> model.tLg(i)>> model.d(i+k);
        std::cout<<"Love wave group period points(s):" << std::endl;
        for(int i=0;i<ntLg;i++) std::cout << model.tLg(i) << " ";
        std::cout<<std::endl; 
    }

    // read regulariztion and damping factors
    //infile >> weight >> damp;  

    // read step length and no. of steps
    infile >> eps(0) >> eps(1) >> iteration(0) >> iteration(1);
    infile >> niter;
    
    infile.close();

    // compute Tickonov regularization
    S.resize(n-2,n);
    S.setZero();
    for(int i=0;i<n-2;i++){
        S(i,i) = 1.0;
        S(i,i+1) = -2.0;
        S(i,i+2) = 1.0;
    }

    return 1;
}

template<class T>
void HamiltonMC<T> :: save_disp(T &model,std::string outdir,VectorXd data,int istep)
{
    int k = 0;
    std ::ofstream file;
    std::string filename;
    filename= "mkdir -p " + outdir;
    int info = system(filename.data());
    if(model.ntRc > 0){
        filename = outdir+"/Rc";
        file.open(filename + std ::to_string(istep));
        for(int i=0;i<model.ntRc;i++){
            file << model.tRc(i) << " " << data(k+i) << std::endl;
        }
        file.close();
        k += model.ntRc;
    }
    if(model.ntRg > 0){
        filename = outdir+"/Rg";
        file.open(filename + std ::to_string(istep));
        for(int i=0;i<model.ntRg;i++){
            file << model.tRg(i) << " " << data(k+i) << std::endl;
        }
        file.close();
        k += model.ntRg;
    }
    if(model.ntLc > 0){
        filename = outdir+"/Lc";
        file.open(filename + std ::to_string(istep));
        for(int i=0;i<model.ntLc;i++){
            file << model.tLc(i) << " " << data(k+i) << std::endl;
        }
        file.close();
        k += model.ntLc;
    }
    if(model.ntLg > 0){
        filename = outdir+"/Lg";
        file.open(filename + std ::to_string(istep));
        for(int i=0;i<model.ntLg;i++){
            file << model.tLg(i) << " " << data(k+i) << std::endl;
        }
        file.close();
    }
}

template<class T>
int HamiltonMC<T>::run(T &model,std::string outdir)
{
    int istep=0,count=0,L;
    bool accept = true;
    MatrixXd kernel(n,nd),ModelHMC(niter,n);
    VectorXd p(n),x(n),data(nd),tmp(n);
    VectorXd xtmp(n),ptmp(n);
    VectorXf v0(n);
    double H0,H1,dt;
    std ::ofstream file;
    H1 = NAN;
    while(istep < niter){
        count += 1 ;
        dt = uniform(eps(0),eps(1));
        L = uniformi(iteration(0),iteration(1));
        std::cout<<istep+1<<"-th search step"<<std::endl;
        std::cout<<"L="<<L<<" dt="<<dt<<std::endl;

        // generate momentum with gaussian distribution
        for(int i=0;i<n;i++){
            p(i) = gaussian()*1.0e-2 ;
        }

        // renew x and p via Hamiltonian update
        // leap-frog scheme is used
        if(accept){ // if the new point is accepted
            x = model.v.template cast<double>();
            v0 = model.v * 1.0;

            // compute current Hamiltonian
            if(isnan(H1)){
                H0 = potential(model) + kinetic(p);
                tmp = gradient(model);
            }
            else{
                H0 = H1;
            }
        }
        else{ // if the new point is not accepted 
            model.v = v0;
            x = model.v.template cast<double>();
            tmp = gradient(model);  
        }
        
        // leap-frog scheme
        p -= 0.5 * dt* tmp;
        for(int i = 0; i <L;i++){
            std::cout<<"leap frog evolution: "<<i+1<<"-th step"<<std::endl;
            x += dt * p;
            model.v = x.cast<float>();
            xtmp = x;
            ptmp = p;

            // tackle boundaries
            for(int k=0;k<n;k++){
                double up = boundary(k,1),low = boundary(k,0);
                while(xtmp(k) > up || xtmp(k) <low){
                    if(xtmp(k) >up){
                        xtmp(k) = 2 * up - xtmp(k);
                        ptmp(k) = - ptmp(k);
                    }
                    if(xtmp(k) < low){
                        xtmp(k) = 2* low  - xtmp(k);
                        ptmp(k) = - ptmp(k);
                    }
                }
            }
            x = xtmp;
            p = ptmp;
            model.v = x.cast<float>();
            if(i != L-1){
                tmp = gradient(model); 
                p -= dt * tmp;
            }
        }
       //Make a half step for momentum at the end.
       tmp = gradient(model);
       p -= 0.5 * dt * tmp;

       // negate p to make proposal symmetric
       p = -p;

       // compute new Hamiltonian
       H1 = potential(model) + kinetic(p);

       // Accept or reject the state at end of trajectory, returning either
       // the position at the end of the trajectory or the initial position 
       if(uniform() < exp(-H1 + H0)){
           std::cout << H1 << " " << H0 << std::endl;
           ModelHMC.row(istep) = x;
           istep += 1;
           std::cout << "iteration " << istep << "-th model is accepted" << std::endl;
           accept = true;

           //save the accepted models
           std ::string filename = outdir+"/iter";
           file.open(filename + std ::to_string(istep));
           for(int j=0;j<n;j++){
                file << model.thk(j) << " " << x(j) << std::endl;
           }
           file.close();

           // save dispersion curves
           data = model.forward();
           save_disp(model,outdir,data,istep);
       }
       else{
           accept = false;
       }
        
    }
    std::cout << "the accepting rate is " << 100. * niter / count << "%"<< std::endl;

    // compute mean model and corresponding dispersion curve
    x = ModelHMC.colwise().sum() / niter;
    model.v=x.cast<float>();
    data = model.forward();
    save_disp(model,".",data);

    return 0;
}
