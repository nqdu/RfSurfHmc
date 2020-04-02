#include"celldisp.h"
#include"HMC.h"

int main(){
    LayerModel model;
    HamiltonMC<LayerModel> hmc;

    // read initial model
    std::string filename = "input/hmc.in",outdir="Models";
    hmc.read_input(model,filename);
    hmc.run(model,outdir);
    return 0;
}