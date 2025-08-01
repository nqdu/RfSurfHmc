import matplotlib.pyplot as plt 
import numpy as np 
import os
import yaml


from src.plot_results import Obs_data
from src.plot_results import plotFig

# load yaml files
with open("param.yaml","r") as f:
    param = yaml.safe_load(f)

# set para
plotFigure = plotFig(
    HMC_output_path=param['hmc']['OUTPUT_DIR'],
    chain_name=param['hmc']['name'],
    figure_output_path=param['hmc']['OUTPUT_DIR'])
plotFigure.obs_data.set_t(param['swd']['tRc'],"Rc")
plotFigure.obs_data.set_t(param['swd']['tRg'],"Rg")

t = np.arange(param['rf']['nt']) * param['rf']['dt'] - param['rf']['time_shift']
plotFigure.obs_data.set_t(t,"rf")
# set_t for setting the time series

real_data = np.load(f"{param['hmc']['OUTPUT_DIR']}/real_syn.npy")
plotFigure.update_real_data(real_data)


thk_real = np.array(param['true_model']['thk'])
vs_real = np.array(param['true_model']['vs'])
#thk_real = np.array([5.,10.,16.,0.0]) 
#vs_real = np.array([3.1,3.64,3.87,4.5])
plotFigure.update_refmodel(np.hstack((vs_real,thk_real)))




# set the search range of inversion
boundaries = np.ones((len(thk_real)*2,2))
for i in range(len(thk_real)):
    boundaries[i,0] = vs_real[i] - vs_real[i]*0.8
    boundaries[i,1] = vs_real[i] + vs_real[i]*0.8       
    if boundaries[i,0] < 1.5:
        boundaries[i,0] = 1.5
    if boundaries[i,1] > 5:
        boundaries[i,1] = 5
        
    boundaries[i + len(thk_real),0] = thk_real[i] - thk_real[i]*0.2
    boundaries[i + len(thk_real),1] = thk_real[i] + thk_real[i]*0.2     
boundaries[-1,:] = 0.0,2.0


#plotFigure.update_initmodel(np.hstack((vs_init,thk_init)))


figure = plotFigure.plot_bestmodels(np.array([0,80,200]))    
# np.array([0,65,200]) means the depth interpolation
# from 0 to 65km, 200 points

plotFigure.savefig(figure,"best_model.png")

figure = plotFigure.plot_best_fit()
plotFigure.savefig(figure,"best_syn_fit.png")

figure = plotFigure.plot_bestmodels_hit(np.array([2,5,100]),np.array([0,80,100]))
# np.array([1.8,4,100]) velocity 1.8km/s ~ 4km/s 100 points
plotFigure.savefig(figure,"best_model_hist.png")

figure = plotFigure.plot_chain_misfit()
plotFigure.savefig(figure,"misfit.png")



figure = plotFigure.plot_best_fit_hist(np.array([-0.2,0.45,100]),np.array([2.5,4,100]),np.array([2,4,100]))
# np.array([-0.1,0.2,100])  amplitude of rf 
# np.array([2,3,100])       velocity of Rc
# np.array([1.6,2.4,100])    rg 
plotFigure.savefig(figure,"plot_best_fit_hist.png")

figure = plotFigure.plot_var_statisc(boundaries, thk_real, vs_real)
plotFigure.savefig(figure,"static.png")

