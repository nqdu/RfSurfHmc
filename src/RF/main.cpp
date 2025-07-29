#include"rf_cal.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <tuple>


namespace py = pybind11;
using py::arg;

const auto FCST =  py::array::c_style | py::array::forcecast ;

typedef py::array_t<double,FCST> dmat;
typedef py::array_t<float> fmat;


dmat forward(dmat &thk, dmat &rho, dmat &vp, 
            dmat &vs,dmat &qa,dmat &qb,
            double ray_p, int nt, double dt, double gauss,
            double time_shift, const std::string &method = "time", 
            double water = 0.001,
            const std::string &rf_type = "P")
{
    int nlayer = thk.size();
    int int_rf_type = 0;
    dmat rcv_fun(nt);
    
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

    // compute receiver functions
    if(method == "time") {
        cal_rf_time_(
            thk.data(),vp.data(),vs.data(),
            rho.data(),qa.data(),qb.data(),
            nlayer,nt,dt,ray_p,gauss,time_shift,
            int_rf_type,rcv_fun.mutable_data()
        );
    }
    else {
        cal_rf_freq_(
            thk.data(),vp.data(),vs.data(),
            rho.data(),qa.data(),qb.data(),
            nlayer,nt,dt,ray_p,gauss,time_shift,
            water,int_rf_type,rcv_fun.mutable_data()
        );
    }

    return rcv_fun;
}

auto kernel(dmat &thk, dmat &rho, dmat &vp, 
                dmat &vs,dmat &qa,dmat &qb,
            double ray_p, int nt, double dt, double gauss,
            double time_shift,const std::string &method = "time",
            double water = 0.001,const std::string &rf_type = "P", 
            const std::string &par_type= "vs")
{
    int nlayer = thk.size();
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

    dmat rcv_fun(nt),rcv_fun_p({nlayer,nt});

    if(method == "time") {
        cal_rf_par_time_(
            thk.data(),vp.data(),vs.data(),
            rho.data(),qa.data(),qb.data(),
            nlayer,nt,dt,ray_p,gauss,time_shift,
            int_rf_type,int_par_type,rcv_fun.mutable_data(),
            rcv_fun_p.mutable_data()
        );
    }
    else {
        cal_rf_par_freq_(
            thk.data(),vp.data(),vs.data(),
            rho.data(),qa.data(),qb.data(),
            nlayer,nt,dt,ray_p,gauss,time_shift,water,
            int_rf_type,int_par_type,rcv_fun.mutable_data(),
            rcv_fun_p.mutable_data()
        );
    }

    return std::make_tuple(rcv_fun,rcv_fun_p);
       
}



auto kernel_all(dmat &thk, dmat &rho, dmat &vp, 
                dmat &vs,dmat &qa,dmat &qb,
            double ray_p, int nt, double dt, double gauss,
            double time_shift,const std::string &method = "time",
            double water = 0.001,const std::string &rf_type = "P"
)
{
    int nlayer = thk.size();
    int int_rf_type = 0;

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

    
    dmat rcv_fun(nt),rcv_fun_p({4,nlayer,nt});

    if(method == "time") {
        cal_rf_par_time_all_(
            thk.data(),vp.data(),vs.data(),
            rho.data(),qa.data(),qb.data(),
            nlayer,nt,dt,ray_p,gauss,time_shift,
            int_rf_type,rcv_fun.mutable_data(),
            rcv_fun_p.mutable_data()
        );
    }
    else {
        cal_rf_par_freq_all_(
            thk.data(),vp.data(),vs.data(),
            rho.data(),qa.data(),qb.data(),
            nlayer,nt,dt,ray_p,gauss,time_shift,water,
            int_rf_type,rcv_fun.mutable_data(),
            rcv_fun_p.mutable_data()
        );
    }

    return std::make_tuple(rcv_fun,rcv_fun_p);
       
}

PYBIND11_MODULE(librf,m){
    m.doc() = "Receiver function and partial derivative\n";
    m.def("forward",&forward,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"),arg("qa"),arg("qb"), arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"),
          arg("method") = "time", arg("water") = 0.001,
          arg("rf_type") = "P",
          "receiver function calculating c++ wrapper");

    m.def("kernel",&kernel,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"), arg("qa"),arg("qb"),arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"),
          arg("method") = "time", arg("water") = 0.001,
          arg("rf_type") = "P", arg("par_type") = "vs",
          "receiver function kernel calculating c++ wrapper");

    m.def("kernel_all",&kernel_all,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"), arg("qa"),arg("qb"),arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"),
          arg("method") = "time", arg("water") = 0.001,
          arg("rf_type") = "P"
          "receiver function kernel calculating c++ wrapper");
}
