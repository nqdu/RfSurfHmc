#include"rf_cal.hpp"
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include<iostream>
namespace py = pybind11;
using py::arg;

typedef py::array_t<double> dtensor;
typedef py::array_t<float> ftensor;
typedef std::tuple<dtensor,dtensor> tuple2;

dtensor forward(dtensor &thk, dtensor &rho, dtensor &vp, dtensor &vs,
            double ray_p, int nt, double dt, double gauss,
            double time_shift, double water_level, std::string rf_type)
{
    int nlayer = thk.size();
    int int_rf_type = 0;
    int int_par_type = 0;
    dtensor rcv_fun(nt);
    
    if( (rf_type == "P") || (rf_type == "p") )
    {
        int_rf_type = 1;
    }
    else if ( (rf_type == "S") || (rf_type == "s") )
    {
        int_rf_type = 2;
        //time_shift = - time_shift;
    }
    else
    {
        std::cout << "rf_type should be one of [P,p,S,s]" << std::endl;
        exit(-1);
    }
    
    receiverFunc_(thk.mutable_data(),vp.mutable_data(),vs.mutable_data(),
                 rho.mutable_data(),nlayer,nt,dt,ray_p,gauss,time_shift,
                 water_level,int_rf_type,rcv_fun.mutable_data());
    
    
    
    return rcv_fun;

}


tuple2 adjoint_kernel(dtensor &thk, dtensor &rho, dtensor &vp, dtensor &vs,
            double ray_p, int nt, double dt, double gauss,
            double time_shift, double water_level, std::string rf_type, 
            std::string par_type)
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
       // time_shift =  time_shift;
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

    dtensor rcv_fun_p(nt*nlayer), rcvfun(nt);

    receiverFunc_par_(thk.mutable_data(),vp.mutable_data(),vs.mutable_data(),
                    rho.mutable_data(),nlayer,nt,dt,ray_p,gauss,time_shift,
                    water_level,int_rf_type,int_par_type,rcvfun.mutable_data(),
                    rcv_fun_p.mutable_data());

    rcv_fun_p.resize({nlayer,nt});

    // maketuple
    auto t2 = std::make_tuple(rcvfun,rcv_fun_p);

    return t2;
   
}

PYBIND11_MODULE(librf,m){
    m.doc() = "Receiver function and partial derivative\n";
    m.def("forward",&forward,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"), arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"),
          arg("water_level"), arg("rf_type"),
          "receiver function calculating c++ wrapper");
    m.def("adjoint_kernel",&adjoint_kernel,arg("thk"), arg("rho"), arg("vp"),
          arg("vs"), arg("ray_p"), arg("nt"),
          arg("dt") ,arg("gauss"),arg("time_shift"),
          arg("water_level"), arg("rf_type"), arg("par_type"),
          "receiver function kernel calculating c++ wrapper");
}
