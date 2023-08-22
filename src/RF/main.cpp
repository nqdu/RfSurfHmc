#include"rf_cal.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <tuple>
namespace py = pybind11;
using py::arg;


typedef py::array_t<double> dtensor;
typedef py::array_t<float> ftensor;

dtensor rf_time(dtensor &thk_in, dtensor &rho_in, dtensor &vp_in, 
                dtensor &vs_in,dtensor &qa_in,dtensor &qb_in,
                double ray_p, int nt, double dt, double gauss,
                double time_shift, const std::string &rf_type)
{
    int nlayer = thk_in.size();
    int int_rf_type = 0;
    dtensor rcv_fun(nt);
    
    if( (rf_type == "P") || (rf_type == "p") )
    {
        int_rf_type = 1;
    }
    else if ( (rf_type == "S") || (rf_type == "s") )
    {
        int_rf_type = 2;
        time_shift = - time_shift;
    }
    else
    {
        std::cout << "rf_type should be one of [P,p,S,s]" << std::endl;
        exit(-1);
    }

    // local arrays
    double vp[nlayer],vs[nlayer],thk[nlayer],rho[nlayer];
    double qa[nlayer],qb[nlayer];

    // copy input arrays
    auto thk0 = thk_in.unchecked<1>(), vp0 = vp_in.unchecked<1>();
    auto vs0 = vs_in.unchecked<1>(), rho0 = rho_in.unchecked<1>();
    auto qa0 = qa_in.unchecked<1>(), qb0 = qb_in.unchecked<1>();
    for(int i = 0; i < nlayer; i ++){
        vs[i] = vs0(i); vp[i] = vp0(i);
        rho[i] = rho0(i); thk[i] = thk0(i);
        qa[i] = qa0(i); qb[i] = qb0(i);
    }
    
    cal_rf_time_(thk,vp,vs,rho,qa,qb,nlayer,nt,dt,ray_p,gauss,time_shift,
                 int_rf_type,rcv_fun.mutable_data());
    
    return rcv_fun;

}


dtensor rf_freq(dtensor &thk_in, dtensor &rho_in, dtensor &vp_in, 
                dtensor &vs_in,dtensor &qa_in,dtensor &qb_in,
                double ray_p, int nt, double dt, double gauss,
                double time_shift,double water, const std::string &rf_type)
{
    int nlayer = thk_in.size();
    int int_rf_type = 0;
    dtensor rcv_fun(nt);
    
    if( (rf_type == "P") || (rf_type == "p") )
    {
        int_rf_type = 1;
    }
    else if ( (rf_type == "S") || (rf_type == "s") )
    {
        int_rf_type = 2;
        time_shift = - time_shift;
    }
    else
    {
        std::cout << "rf_type should be one of [P,p,S,s]" << std::endl;
        exit(-1);
    }

    // local arrays
    double vp[nlayer],vs[nlayer],thk[nlayer],rho[nlayer];
    double qa[nlayer],qb[nlayer];

    // copy input arrays
    auto thk0 = thk_in.unchecked<1>(), vp0 = vp_in.unchecked<1>();
    auto vs0 = vs_in.unchecked<1>(), rho0 = rho_in.unchecked<1>();
    auto qa0 = qa_in.unchecked<1>(), qb0 = qb_in.unchecked<1>();
    for(int i = 0; i < nlayer; i ++){
        vs[i] = vs0(i); vp[i] = vp0(i);
        rho[i] = rho0(i); thk[i] = thk0(i);
        qa[i] = qa0(i); qb[i] = qb0(i);
    }
    
    cal_rf_freq_(thk,vp,vs,rho,qa,qb,nlayer,nt,dt,ray_p,gauss,time_shift,
                 water,int_rf_type,rcv_fun.mutable_data());
    
    return rcv_fun;

}


auto rf_par_time(dtensor &thk_in, dtensor &rho_in, dtensor &vp_in, 
                dtensor &vs_in,dtensor &qa_in,dtensor &qb_in,
            double ray_p, int nt, double dt, double gauss,
            double time_shift,std::string rf_type, 
            std::string par_type)
{
    int nlayer = thk_in.size();
    int int_rf_type = 0;
    int int_par_type = 0;

    if( (rf_type == "P") || (rf_type == "p") )
    {
        int_rf_type = 1;
    }
    else if ( (rf_type == "S") || (rf_type == "s") )
    {
        int_rf_type = 2;
       time_shift =  -time_shift;
    }
    else
    {
        std::cout << "rf_type should be one of [P,p,S,s]" << std::endl;
        exit(-1);
    }

    
    if( (par_type == "alpha") || (par_type == "vp") )
    {
        int_par_type = 2;
    } 
    else if( (par_type == "beta") || (par_type == "vs") )
    {
        int_par_type = 3;
    }
    else if(par_type == "rho")
    {
        int_par_type = 1;
    }
    else if( (par_type == "thick") || (par_type == "h") )
    {
        int_par_type = 4;
    }
    else
    {
        std::cout << "par_type should be one of [vp,vs,rho,thick]" << std::endl;
        exit(-1);       
    }


    // local arrays
    double vp[nlayer],vs[nlayer],thk[nlayer],rho[nlayer];
    double qa[nlayer],qb[nlayer];

    // copy input arrays
    auto thk0 = thk_in.unchecked<1>(), vp0 = vp_in.unchecked<1>();
    auto vs0 = vs_in.unchecked<1>(), rho0 = rho_in.unchecked<1>();
    auto qa0 = qa_in.unchecked<1>(), qb0 = qb_in.unchecked<1>();
    for(int i = 0; i < nlayer; i ++){
        vs[i] = vs0(i); vp[i] = vp0(i);
        rho[i] = rho0(i); thk[i] = thk0(i);
        qa[i] = qa0(i); qb[i] = qb0(i);
    }
    

    dtensor rcv_fun(nt),rcv_fun_p(nt*nlayer);

    cal_rf_par_time_(thk,vp,vs,rho,qa,qb,nlayer,nt,dt,ray_p,gauss,time_shift,
                    int_rf_type,int_par_type,rcv_fun.mutable_data(),
                    rcv_fun_p.mutable_data());

    rcv_fun_p.resize({nlayer,nt});

    return std::make_tuple(rcv_fun,rcv_fun_p);
   
}

auto rf_par_freq(dtensor &thk_in, dtensor &rho_in, dtensor &vp_in, 
                dtensor &vs_in,dtensor &qa_in,dtensor &qb_in,
                double ray_p, int nt, double dt, double gauss,
                double time_shift,double water,std::string rf_type, 
                std::string par_type)
{
    int nlayer = thk_in.size();
    int int_rf_type = 0;
    int int_par_type = 0;

    if( (rf_type == "P") || (rf_type == "p") )
    {
        int_rf_type = 1;
    }
    else if ( (rf_type == "S") || (rf_type == "s") )
    {
        int_rf_type = 2;
       time_shift =  -time_shift;
    }
    else
    {
        std::cout << "rf_type should be one of [P,p,S,s]" << std::endl;
        exit(-1);
    }

    
    if( (par_type == "alpha") || (par_type == "vp") )
    {
        int_par_type = 2;
    } 
    else if( (par_type == "beta") || (par_type == "vs") )
    {
        int_par_type = 3;
    }
    else if(par_type == "rho")
    {
        int_par_type = 1;
    }
    else if( (par_type == "thick") || (par_type == "h") )
    {
        int_par_type = 4;
    }
    else
    {
        std::cout << "par_type should be one of [vp,vs,rho,thick]" << std::endl;
        exit(-1);       
    }


    // local arrays
    double vp[nlayer],vs[nlayer],thk[nlayer],rho[nlayer];
    double qa[nlayer],qb[nlayer];

    // copy input arrays
    auto thk0 = thk_in.unchecked<1>(), vp0 = vp_in.unchecked<1>();
    auto vs0 = vs_in.unchecked<1>(), rho0 = rho_in.unchecked<1>();
    auto qa0 = qa_in.unchecked<1>(), qb0 = qb_in.unchecked<1>();
    for(int i = 0; i < nlayer; i ++){
        vs[i] = vs0(i); vp[i] = vp0(i);
        rho[i] = rho0(i); thk[i] = thk0(i);
        qa[i] = qa0(i); qb[i] = qb0(i);
    }
    

    dtensor rcv_fun(nt),rcv_fun_p(nt*nlayer);

    cal_rf_par_freq_(thk,vp,vs,rho,qa,qb,nlayer,nt,dt,ray_p,gauss,time_shift,
                    water,int_rf_type,int_par_type,rcv_fun.mutable_data(),
                    rcv_fun_p.mutable_data());

    rcv_fun_p.resize({nlayer,nt});

    return std::make_tuple(rcv_fun,rcv_fun_p);
   
}

PYBIND11_MODULE(librf,m){
    m.doc() = "Receiver function and partial derivative\n";
    m.def("rf_time",&rf_time,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"),arg("qa"),arg("qb"), arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"),
          arg("rf_type"),
          "Time domain receiver function calculating c++ wrapper");
    m.def("rf_freq",&rf_freq,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"),arg("qa"),arg("qb"), arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"), arg("water"),
          arg("rf_type"),
          "Frequency domain receiver function calculating c++ wrapper");

    m.def("rf_par_time",&rf_par_time,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"), arg("qa"),arg("qb"),arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"),
          arg("rf_type"), arg("par_type"),
          "Time domain receiver function kernel calculating c++ wrapper");


    m.def("rf_par_freq",&rf_par_freq,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"), arg("qa"),arg("qb"),arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"),arg("water"),
          arg("rf_type"), arg("par_type"),
          "Frequency domain receiver function kernel calculating c++ wrapper");
}
