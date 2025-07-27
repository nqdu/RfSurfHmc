#include"surfdisp.hpp"
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
namespace py = pybind11;
using py::arg;

const auto FCST =  py::array::c_style | py::array::forcecast ;

typedef py::array_t<float,FCST> ftensor;
typedef py::array_t<double,FCST> dtensor;
typedef std::tuple<dtensor,dtensor,dtensor,dtensor,dtensor,bool> tupledt6;
typedef std::tuple<dtensor,bool> tuple2;

tuple2 forward(ftensor &thk,ftensor &vp,ftensor &vs,ftensor &rho,
                dtensor &t,std::string wavetype,
                int mode=0,bool sphere=false)
{
    // check input parameters
    bool flag = wavetype == "Rc" || wavetype == "Rg" || 
                wavetype == "Lc" || wavetype == "Lg";
    if(flag == false){
        std::cout << "cnew = _flat2sphere(double t,double c,std::string wavetp)\n";
        std::cout << "parameter wavetp should be in one of [Rc,Rg,Lc,Lg]\n ";
        exit(0);
    }

    // allocate space
    int nt = t.size(), n = thk.size();
    dtensor cg(nt);
    bool keep_flat = false;

    // forward computation
    int ierr;
    if(wavetype == "Rc"){
       ierr = _surfdisp(thk.mutable_data(),vp.mutable_data(),vs.mutable_data(),
                 rho.mutable_data(),n,t.mutable_data(),cg.mutable_data(),
                 nt,wavetype,mode,sphere,keep_flat);
    }
    else if (wavetype == "Rg"){
        ierr = _RayleighGroup(thk.mutable_data(),vp.mutable_data(),vs.mutable_data(),
                      rho.mutable_data(),n,t.mutable_data(),cg.mutable_data(),nt,
                       mode,sphere);
    }
    else if (wavetype=="Lc"){
        ierr =  _surfdisp(thk.mutable_data(),vp.mutable_data(),vs.mutable_data(),
                 rho.mutable_data(),n,t.mutable_data(),cg.mutable_data(),
                 nt,wavetype,mode,sphere,keep_flat);
    }
    else{
        ierr = _LoveGroup(thk.mutable_data(),vs.mutable_data(),rho.mutable_data(),
                    n,t.mutable_data(),cg.mutable_data(),nt,mode,sphere);
    }

    bool return_flag = true;
    if(ierr == 1) return_flag = false;
    auto tt = std::make_tuple(cg,return_flag);

    return tt;
}

tupledt6
adjoint_kernel(ftensor &thk,ftensor &vp,ftensor &vs,ftensor &rho,
                dtensor &t,std::string wavetype,int mode=0,
                bool sphere=false)
{
    // allocate space
    int nt = t.size(),n = thk.size();
    dtensor dcda(nt*n),dcdb(nt*n),dcdr(nt*n),dcdh(nt*n),c(nt); 
    dcda.resize({nt,n}), dcdb.resize({nt,n});
    dcdh.resize({nt,n}),dcdr.resize({nt,n});

   int ierr = _SurfKernel(thk.mutable_data(),vp.mutable_data(),vs.mutable_data(),
                rho.mutable_data(),n,t.mutable_data(),c.mutable_data(),
                nt,dcda.mutable_data(),dcdb.mutable_data(),dcdr.mutable_data(),
                dcdh.mutable_data(),wavetype,mode,sphere);
    
    bool return_flag = true;
    if(ierr == 1) return_flag = false;
    auto tt = std::make_tuple(c,dcda,dcdb,dcdr,dcdh,return_flag);

    return tt;
}

PYBIND11_MODULE(libsurf,m){
    m.doc() = "Surface wave dispersion and sensivity kernel\n";
    m.def("forward",&forward,arg("thk"),arg("vp"),arg("vs"),
          arg("rho"),arg("period"),arg("wavetype"),
          arg("mode")=0,arg("sphere")=false,
          "Surface wave dispersion c++ wrapper");
    m.def("adjoint_kernel",&adjoint_kernel,arg("thk"),arg("vp"),arg("vs"),
          arg("rho"),arg("period"),arg("wavetype"),
          arg("mode")=0,arg("sphere")=false,
          "Surface wave dispersion sensitivity kernel c++ wrapper");
}


