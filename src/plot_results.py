import matplotlib.pyplot as plt 
import numpy as np 
import os
import os.path as op
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

EPS = 1e-10


class Obs_data(object):

    '''
    These class is used for storing the observation data, including the rf/swd, and
    also the time series.

    It is a little bit awkward, and it should be done by passing 
    the receiverFunc and SWD class
    '''
    def __init__(self):
        pass
        

    def set_obs(self, data, types):
        if types == "Rc":
            self.Rc = data
        elif types == "Rg":
            self.Rg= data
        elif types == "rf":
            self.rf = data
        else:
            print("wrong types")        

    def set_t(self, time, types):
        if types == "Rc":
            self.tRc = time
            self.ntRc = len(time)
        elif types == "Rg":
            self.tRg= time
            self.ntRg = len(time)
        elif types == "rf":
            self.trf = time
            self.ntrf = len(time)
        else:
            print("wrong types")
    



class plotFig(object):
    '''
    plot figure for hmc output
    '''

    def __init__(self, HMC_output_path = os.getcwd(), chain_num = 0, 
                figure_output_path = os.getcwd(), HMC_calculate_model = None):
        os.chdir(HMC_output_path)
        self.datapath = HMC_output_path
        self.figpath = figure_output_path
        self.chain_num = chain_num

        print('Current data path: %s' % self.datapath)

        if HMC_calculate_model is not None:
            #get parameter from hmc calculate model
            #TODO: to initialize the self.obs_data
            pass
        else:
            self.obs_data = Obs_data()
        
        self.refmodel = None
        self.initmodel = None
    
    def unpack_thk_vs(self, data):
        vs = data[:int(len(data)/2)]
        thk = data[int(len(data)/2):]
        thk[-1] = 0
        depth = np.cumsum(thk)
        assert len(depth)  == len(vs)
        return depth, vs

    def _unpack_thk_vs(self, data):
        vs = data[:int(len(data)/2)]
        thk = data[int(len(data)/2):]
        depth = np.cumsum(thk)
        assert len(depth)  == len(vs)
        return thk, vs

    def unpack_obs(self, data, ntrf, ntRc, ntRg):
        # TODO: it should be done for the input of different types of swd.

        rf = data[ : ntrf]
        Rc = data[ntrf : ntrf + ntRc ]
        Rg = data[ntrf + ntRc : ntrf + ntRc + ntRg] 
     

        return rf, Rc, Rg

    def savefig(self, fig, filename):
        if fig is not None:
            outfile = op.join(self.figpath, filename)
            fig.savefig(outfile, bbox_inches="tight")
            plt.close('all')
    
    def find_info_in_chain(self, chain_num, info_type, index):
        '''
        return the model(thk/vs) or synthetic data of specfic chain/index
        input:
            chain_num   (int):  the index number of chain
            info_type   (str):  model for thk/vs; synthetic for data(rf/swd)
            index       (int):  the index number of model

        return:
            self.depth / self.vs for model
            self.obs_data.rf / self.obs_data_Rc and Rg for syn data
        '''
        chain_data_folder = "chain_joint" + str(chain_num)
        if info_type == "model":
            dat_file_name = "model" + str(index) + ".dat"

            data = np.loadtxt(op.join(self.datapath,chain_data_folder,dat_file_name))
            self.depth , self.vs = self.unpack_thk_vs(data)

        if info_type == "model_init":
            dat_file_name = "initmodel"  + ".dat"

            data = np.loadtxt(op.join(self.datapath,chain_data_folder,dat_file_name))
            self.depth , self.vs = self.unpack_thk_vs(data)
            

        if info_type == "model_mean":
            dat_file_name = "model" + "mean" + ".dat"

            data = np.loadtxt(op.join(self.datapath,chain_data_folder,dat_file_name))
            self.depth , self.vs = self.unpack_thk_vs(data)
            
        if info_type == "synthetic":
            dat_file_name = "synthetic" + str(index) + ".dat"
            data = np.loadtxt(op.join(self.datapath,chain_data_folder,dat_file_name))
            #todo: why there is obs data in synthetic.dat

            self.obs_data.rf, self.obs_data.Rc, self.obs_data.Rg = self.unpack_obs(data[:,1],
            self.obs_data.ntrf, self.obs_data.ntRc, self.obs_data.ntRg)

        if info_type == "synthetic_mean":
            dat_file_name = "synthetic" + "mean" + ".dat"
            data = np.loadtxt(op.join(self.datapath,chain_data_folder,dat_file_name))
            #todo: why there is obs data in synthetic.dat

            self.obs_data.rf, self.obs_data.Rc, self.obs_data.Rg = self.unpack_obs(data[:,1],
            self.obs_data.ntrf, self.obs_data.ntRc, self.obs_data.ntRg)
   
        if info_type == "para":
            dat_file_name = "model" + str(index) + ".dat"

            data = np.loadtxt(op.join(self.datapath,chain_data_folder,dat_file_name))
            self.thk , self.vs = self._unpack_thk_vs(data)

    def update_refmodel(self,data):
        self.refmodel = data 
            
    def update_initmodel(self,data):
        self.initmodel = data 

    def update_real_data(self,data):
        self.real_data = data


    def plot_chain_misfit(self):
        '''
        retrun matplotlib fig

        Plot misfit of different chain.

        '''

        fig, ax = plt.subplots(figsize=(12,5))
        all_misfit = np.loadtxt(op.join(self.datapath,"misfit.dat"))
        if len(all_misfit.shape) == 1:
            all_misfit = all_misfit.reshape(1,len(all_misfit))
        # row for chains and col for model_index 
        max_chain = len(all_misfit[:,0])


        #plot best syn result in different chain
        for chain_iter in range(max_chain):
            misfit = np.min(all_misfit[chain_iter,:])
            
            print("misfit: {}".format(misfit))
            if misfit < EPS:
                print("skip chain {}".format(chain_iter))
                continue
            
            model_index = np.argmin(all_misfit[chain_iter,:])
            ax.plot(all_misfit[chain_iter,:] ,ls='-', lw=0.8, alpha=0.5, 
                    label = "chain %d" %(chain_iter) )

            
        ax.set_ylabel('misfit')
        ax.set_title('misfit curve from %d chains' %max_chain)
        ax.set_xlabel('iter')
        ax.legend()

        return fig

    def plot_para_stats(self):
        '''
        retrun matplotlib fig

        Plot statistic result of different chain.

        '''        



    def plot_refmodel(self, fig, mtype='model', **kwargs):
        if fig is not None and self.refmodel[mtype] is not None:
            if mtype == 'nlays':
                nlays = self.refmodel[mtype]
                fig.axes[0].axvline(nlays, color='red', lw=0.5, alpha=0.7)

            if mtype == 'model':
                dep, vs = self.unpack_thk_vs(self.refmodel['model'])
                assert len(dep) == len(vs)
                fig.axes[0].plot(vs, dep, **kwargs)
                if len(fig.axes) == 2:
                    deps = np.unique(dep)
                    for d in deps:
                        fig.axes[1].axhline(d, **kwargs)
        return fig

    def plot_best_fit(self):
        '''
        retrun matplotlib fig

        Plot the best fit of receiver function and surface wave dispersion.

        '''
        
        fig, ax = plt.subplots(3, 1, figsize=(7,10))


        all_misfit = np.loadtxt(op.join(self.datapath,"misfit.dat"))
        if len(all_misfit.shape) == 1:
            all_misfit = all_misfit.reshape(1,len(all_misfit))
        max_chain = len(all_misfit[:,0])

        thebestmisfit = 1e15
        thebestchain = 0
        thebestindex = 0

        #plot best syn result in different chain
        for chain_iter in range(max_chain):
            misfit = np.min(all_misfit[chain_iter,:])
            if misfit < EPS:
                continue
            model_index = np.argmin(all_misfit[chain_iter,:])
            
            if misfit < thebestmisfit:
                thebestmisfit = misfit
                thebestchain = chain_iter
                thebestindex = model_index
            
            self.find_info_in_chain(chain_iter,"synthetic",model_index)

            if chain_iter == 0:
                ax[0].plot(self.obs_data.trf, self.obs_data.rf , color='k', ls='-', 
                    lw=0.8, alpha=0.5)

                ax[2].plot(self.obs_data.tRg, self.obs_data.Rg , color='k', ls='-', 
                    lw=0.8, alpha=0.5)
                
                ax[1].plot(self.obs_data.tRc, self.obs_data.Rc , color='k', ls='-', 
                    lw=0.8, alpha=0.5)
            else:
                ax[0].plot(self.obs_data.trf, self.obs_data.rf , color='k', ls='-', 
                lw=0.8, alpha=0.5)
                ax[2].plot(self.obs_data.tRg, self.obs_data.Rg , color='k', ls='-', 
                lw=0.8, alpha=0.5)
                ax[1].plot(self.obs_data.tRc, self.obs_data.Rc , color='k', ls='-', 
                lw=0.8, alpha=0.5)

        #plot best syn result 
        self.find_info_in_chain(thebestchain,"synthetic",thebestindex)
        ax[0].plot(self.obs_data.trf, self.obs_data.rf , color='b', ls='-', 
        lw=1.5, alpha=0.5, label = 'best model')
        ax[2].plot(self.obs_data.tRg, self.obs_data.Rg , color='b', ls='-', 
        lw=1.5, alpha=0.5, label = 'best model')
        ax[1].plot(self.obs_data.tRc, self.obs_data.Rc , color='b', ls='-', 
        lw=1.5, alpha=0.5, label = 'best model')

        #plot real data
        if self.real_data is not None:
            self.obs_data.rf, self.obs_data.Rc,  self.obs_data.Rg = self.unpack_obs(
                self.real_data,self.obs_data.ntrf, self.obs_data.ntRc, self.obs_data.ntRg)
        ax[0].plot(self.obs_data.trf, self.obs_data.rf , color='r', ls='-', 
            lw=1.5, alpha=0.5, label = 'real model')  
        ax[2].plot(self.obs_data.tRg, self.obs_data.Rg , color='r', ls='-', 
            lw=1.5, alpha=0.5, label = 'real model')  
        ax[1].plot(self.obs_data.tRc, self.obs_data.Rc , color='r', ls='-', 
            lw=1.5, alpha=0.5, label = 'real model')  

        
        ax[0].set_ylabel('Amplitude')
        #ax[0].set_title('Best rf result from %d chains' %max_chain)
        ax[2].set_ylabel('Velocity')
        #ax[2].set_title('Best Rg result from %d chains' %max_chain)
        ax[1].set_ylabel('Velocity')
        #ax[1].set_title('Best Rc result from %d chains' %max_chain)

        for i in range(len(ax)):
            ax[i].set_xlabel('t(s)')
            ax[i].legend()
        return fig 


    def plot_var_statisc(self, boundaries, thk, vs):

        fig_num = int(np.ceil(np.sqrt(len(boundaries))))
        fig, ax = plt.subplots(fig_num, fig_num, figsize=(15,15))

        all_misfit = np.loadtxt(op.join(self.datapath,"misfit.dat"))
        if len(all_misfit.shape) == 1:
            all_misfit = all_misfit.reshape(1,len(all_misfit))
        thebestmisfit = 1e15
        thebestchain = 0
        thebestindex = 0


        max_chain = len(all_misfit[:,0])
        max_model = len(all_misfit[0,:])
        parameter_var = np.zeros((max_chain*max_model,len(boundaries)))
        for chain_iter in range(max_chain):
            for model_iter in range(max_model):        
                self.find_info_in_chain(chain_iter,"para",model_iter)
                
                if all_misfit[chain_iter,model_iter] < EPS:
                    continue
                
                if all_misfit[chain_iter,model_iter] < thebestmisfit:
                    thebestmisfit = all_misfit[chain_iter,model_iter] 
                    thebestchain = chain_iter
                    thebestindex = model_iter

                for i in range(len(self.thk)):
                    parameter_var[chain_iter*max_model+model_iter,i] = self.vs[i]
                    parameter_var[chain_iter*max_model+model_iter,i+len(self.thk)] = self.thk[i]

        self.find_info_in_chain(thebestchain,"para",thebestindex)

        
        for i in range(fig_num):
            for j in range(fig_num):
                iter_num_of_var = i*fig_num + j
                if (iter_num_of_var) < len(self.thk) * 2:
                
                    if (iter_num_of_var) == len(self.thk) * 2 - 1 :
                        continue                
                
                    ax[i][j].hist((parameter_var[:,iter_num_of_var])[parameter_var[:,iter_num_of_var] != 0])
                    


                    
                    if iter_num_of_var < len(self.vs):
                        ax[i][j].axvline(self.vs[iter_num_of_var],color = 'blue',label = "best")
                        ax[i][j].axvline(vs[iter_num_of_var],color = 'red',label = "real")
                    #    ax[i][j].axvline(boundaries[iter_num_of_var,0],color = 'black',label = "range")
                    #    ax[i][j].axvline(boundaries[iter_num_of_var,1],color = 'black')
                        ax[i][j].legend()
                        ax[i][j].set_xlabel("vs of layer {}".format(iter_num_of_var+1))
                    else:
                        ax[i][j].axvline(self.thk[iter_num_of_var - len(self.thk)],color = 'blue',label = "best")
                        ax[i][j].axvline(thk[iter_num_of_var - len(self.thk)],color = 'red',label = "real")
                    #    ax[i][j].axvline(boundaries[iter_num_of_var,0],color = 'black',label = "range")
                    #    ax[i][j].axvline(boundaries[iter_num_of_var,1],color = 'black')
                        ax[i][j].legend()
                        ax[i][j].set_xlabel("thk of layer {}".format(iter_num_of_var +1- len(self.thk)))
                else:
                    break

                

        return fig

    def plot_best_fit_hist(self, amplitude_range, velRc_range, velRg_range ):
        '''
        retrun matplotlib fig

        Plot the best fit of receiver function and surface wave dispersion.

        '''
        fig, ax = plt.subplots(3, 1, figsize=(7,10))

        velRc_vector = np.linspace(velRc_range[0],velRc_range[1],int(velRc_range[2]))
        velRg_vector = np.linspace(velRg_range[0],velRg_range[1],int(velRg_range[2]))
        amp_vector = np.linspace(amplitude_range[0],amplitude_range[1],int(amplitude_range[2]))
        

        

        hist_grid_swdRc = np.zeros((len(self.obs_data.tRc),len(velRc_vector)))
        hist_grid_swdRg = np.zeros((len(self.obs_data.tRg),len(velRg_vector)))
        hist_grid_rf = np.zeros((len(self.obs_data.trf),len(amp_vector)))

        thebestmisfit = 1e15
        thebestchain = 0
        thebestindex = 0

        all_misfit = np.loadtxt(op.join(self.datapath,"misfit.dat"))
        if len(all_misfit.shape) == 1:
            all_misfit = all_misfit.reshape(1,len(all_misfit))
        max_chain = len(all_misfit[:,0])


        for chain_iter in range(max_chain):
        
            temp = np.argwhere(all_misfit[chain_iter,:])
            mean_chain = np.mean(temp)
            cont = 0
            for model_iter in range(len(all_misfit[chain_iter,:])):

                if all_misfit[chain_iter,model_iter] > mean_chain:
                    cont += 1
                    continue    

                
                self.find_info_in_chain(chain_iter,"synthetic",model_iter)
                self.obs_data.rf, self.obs_data.Rc,  self.obs_data.Rg

                for i in range(len(self.obs_data.rf)):
                    amp_index = np.searchsorted(amp_vector,self.obs_data.rf[i])
                    if amp_index >= len(self.obs_data.rf) :
                        amp_index = len(self.obs_data.rf) - 1                    
                    hist_grid_rf[i,amp_index] += 1

                for i in range(len(self.obs_data.Rg)):
                    vel_index = np.searchsorted(velRg_vector,self.obs_data.Rg[i])
                    if vel_index >= len(velRg_vector) :
                        vel_index = len(velRg_vector) - 1                    
                    hist_grid_swdRg[i,vel_index] += 1
                
                for i in range(len(self.obs_data.Rc)):
                    vel_index = np.searchsorted(velRc_vector,self.obs_data.Rc[i])
                    if vel_index >= len(velRc_vector) :
                        vel_index = len(velRc_vector) - 1
                    hist_grid_swdRc[i,vel_index] += 1

            break


        for i in range(len(hist_grid_rf[:,0])):
            hist_grid_rf[i,:] = hist_grid_rf[i,:]  / np.max(hist_grid_rf[i,:])
        for i in range(len(hist_grid_swdRc[:,0])):
            hist_grid_swdRc[i,:] = hist_grid_swdRc[i,:]  / np.max(hist_grid_swdRc[i,:])
        for i in range(len(hist_grid_swdRg[:,0])):
            hist_grid_swdRg[i,:] = hist_grid_swdRg[i,:]  / np.max(hist_grid_swdRg[i,:])

        ax[0].pcolor(self.obs_data.trf, amp_vector, hist_grid_rf.T, cmap=plt.cm.jet, shading='auto'  )
        ax[1].pcolor(self.obs_data.tRc, velRc_vector, hist_grid_swdRc.T, cmap=plt.cm.jet, shading='auto'  )
        ax[2].pcolor(self.obs_data.tRg, velRg_vector, hist_grid_swdRg.T, cmap=plt.cm.jet, shading='auto'  )
        

        
        all_misfit = np.loadtxt(op.join(self.datapath,"misfit.dat"))
        if len(all_misfit.shape) == 1:
            all_misfit = all_misfit.reshape(1,len(all_misfit))
        max_chain = len(all_misfit[:,0])

        thebestmisfit = 1e15
        thebestchain = 0
        thebestindex = 0

        #plot best syn result in different chain
        for chain_iter in range(max_chain):
            misfit = np.min(all_misfit[chain_iter,:])
            model_index = np.argmin(all_misfit[chain_iter,:])

            if misfit < EPS:
                continue 
                
            if misfit < thebestmisfit :
                thebestmisfit = misfit
                thebestchain = chain_iter
                thebestindex = model_index
            
        #plot best syn result 
        self.find_info_in_chain(thebestchain,"synthetic",thebestindex)
        #ax[0].plot(self.obs_data.trf, self.obs_data.rf , color='b', ls='-', lw=1.5, alpha=0.5, label = "best model")
        #ax[2].plot(self.obs_data.tRg, self.obs_data.Rg , color='b', ls='-', lw=1.5, alpha=0.5, label = "best model")
        #ax[1].plot(self.obs_data.tRc, self.obs_data.Rc , color='b', ls='-', lw=1.5, alpha=0.5, label = "best model")

        self.find_info_in_chain(thebestchain,"synthetic_mean",thebestindex)
        ax[0].plot(self.obs_data.trf, self.obs_data.rf , color='black', ls='-', lw=2, alpha=1, label = "mean model")
        ax[2].plot(self.obs_data.tRg, self.obs_data.Rg , color='black', ls='-', lw=2, alpha=1, label = "mean model")
        ax[1].plot(self.obs_data.tRc, self.obs_data.Rc , color='black', ls='-', lw=2, alpha=1, label = "mean model")
        
        #plot real data
        if self.real_data is not None:
            self.obs_data.rf, self.obs_data.Rc,  self.obs_data.Rg = self.unpack_obs(self.real_data,
                    self.obs_data.ntrf, self.obs_data.ntRc, self.obs_data.ntRg)
        ax[0].plot(self.obs_data.trf, self.obs_data.rf , color='white', ls='-', lw=2, alpha=1, label = "real model")  
        ax[2].plot(self.obs_data.tRg, self.obs_data.Rg , color='white', ls='-', lw=2, alpha=1, label = "real model")  
        ax[1].plot(self.obs_data.tRc, self.obs_data.Rc , color='white', ls='-', lw=2, alpha=1, label = "real model")  

        ax[0].set_ylabel('Amplitude')
        #ax[0].set_title('Best rf result from %d chains' %max_chain)
        ax[2].set_ylabel('Velocity')
        #ax[2].set_title('Best Rg result from %d chains' %max_chain)
        ax[1].set_ylabel('Velocity')
        #ax[1].set_title('Best Rc result from %d chains' %max_chain)
        
        for i in range(len(ax)):
            ax[i].legend()

        
        return fig    
        ax[-1].set_xlabel('t(s)')    

    def plot_bestmodels(self, depth):
        """Return fig.

        Plot the best (fit) models ever discovered per each chain,
        ignoring outliers.
        """
        fig, ax = plt.subplots(figsize=(4, 6.5))


        depth_vector = np.linspace(depth[0],depth[1],int(depth[2]))

        thebestmisfit = 1e15
        thebestchain = 0
        thebestindex = 0

        all_misfit = np.loadtxt(op.join(self.datapath,"misfit.dat"))
        if len(all_misfit.shape) == 1:
            all_misfit = all_misfit.reshape(1,len(all_misfit))
        max_chain = len(all_misfit[:,0])



        #plot best model in different chain
        for chain_iter in range(max_chain):
            misfit = np.min(all_misfit[chain_iter,:])
            model_index = np.argmin(all_misfit[chain_iter,:])

            if misfit < EPS:
                continue 
                            
            if misfit < thebestmisfit:
                thebestmisfit = misfit
                thebestchain = chain_iter
                thebestindex = model_index
            
            self.find_info_in_chain(chain_iter,"model",model_index)
            self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)

            ax.step(self.vs, self.depth , color='k', ls='-', lw=0.8, alpha=0.5)

        # plot ref_model()
        if self.refmodel is not None:
            self.depth , self.vs = self.unpack_thk_vs(self.refmodel)
            self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)
            ref_legend = ax.step(self.vs, self.depth , color='r', ls='-', lw=1.5, alpha=0.5, label = 'real model')

        # plot init_model()
        '''
        if self.initmodel is not None:
            self.depth , self.vs = self.unpack_thk_vs(self.initmodel)
            self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)
            ax.step(self.vs, self.depth , color='y', ls='-', lw=1.5, alpha=0.5, label = 'init model')        
        '''
        #plot best model 
        self.find_info_in_chain(thebestchain,"model",thebestindex)
        self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)
        bst_legend = ax.step(self.vs, self.depth , color='b', ls='-', lw=1.5, alpha=0.5, label = 'best model')
        
        ax.legend()

        ax.invert_yaxis()
        ax.set_xlabel('$V_S$ in km/s')
        ax.set_ylabel('Depth in km')

        ax.set_title('Best fit models from %d chains' %max_chain)

        return fig

    def plot_bestmodels_hit(self, vel, depth):
        """Return fig.

        Plot the best (fit) models ever discovered per each chain,
        ignoring outliers.
        """
        fig, ax = plt.subplots(figsize=(8, 12.5))

        vel_vector = np.linspace(vel[0],vel[1],int(vel[2]))
        depth_vector = np.linspace(depth[0],depth[1],int(depth[2]))
        hist_grid = np.zeros((len(depth_vector),len(vel_vector)))


        thebestmisfit = 1e15
        thebestchain = 0
        thebestindex = 0

        all_misfit = np.loadtxt(op.join(self.datapath,"misfit.dat"))
        if len(all_misfit.shape) == 1:
            all_misfit = all_misfit.reshape(1,len(all_misfit))
        max_chain = len(all_misfit[:,0])


        for chain_iter in range(max_chain):
            mean_chain = np.mean(np.argwhere(all_misfit[chain_iter,:]))/0.1
            for model_iter in range(len(all_misfit[chain_iter,:])):
                self.find_info_in_chain(chain_iter,"model",model_iter)

                if all_misfit[chain_iter,model_iter] > mean_chain:
                    continue    

                new_depth = depth_vector
                self.depth, self.vs = self.interp(self.depth, self.vs, new_depth)
                
                # fill the gap of vs
              
                
               # print(self.depth, self.vs)

          
                
                for i in range(len(self.vs)):
                    vel_index = np.searchsorted(vel_vector,self.vs[i])
                    depth_index = np.searchsorted(depth_vector,self.depth[i])
                    
                    
                    if i >= 1:
                        if (self.vs[i] - self.vs[i-1] < 0.001 ):
                            # no gaps, discard
                            pass
                        else:
                            #print(np.searchsorted(vel_vector,self.vs[i]),np.searchsorted(vel_vector,self.vs[i-1]))
                            vs_index_min = min(np.searchsorted(vel_vector,self.vs[i]),
                            np.searchsorted(vel_vector,self.vs[i-1]))
                            vs_index_max = max(np.searchsorted(vel_vector,self.vs[i]),
                            np.searchsorted(vel_vector,self.vs[i-1]))
                            
                            for j in range(vs_index_min,vs_index_max):
                                #print("junliu debug : vs {}".format(j))
                                hist_grid[depth_index,j] += 1
                                #pass
                                
                                                        
                    
                    
                    if depth_index == int(vel[2]):
                        continue
                    if vel_index ==  int(depth[2]):
                        continue
                    
                    hist_grid[depth_index,vel_index] += 1

        #print("junliu debug: ",len(hist_grid[:0]))        
        for i in range(len(hist_grid[:,0])):
            hist_grid[i,:] = hist_grid[i,:] / np.max(hist_grid[i,:])
            #sprint(hist_grid[i,:])
                    
        clist = ['red','white','blue' ]
        newcmp = LinearSegmentedColormap.from_list('chao',clist)
        map = ax.pcolor(vel_vector, depth_vector, hist_grid, cmap=plt.cm.jet, shading='auto' )
        #colorbar_pos = fig.add_axes([0.15,0.25,0.35,0.03])
        colorbar_h = plt.colorbar(mappable=map,cax=None,ax=None,shrink=1,fraction=0.05)
        colorbar_h.ax.tick_params(labelsize = 16)
        colorbar_h.set_label("Probability",fontsize = 16)
        
        #plot best model in different chain
        for chain_iter in range(max_chain):
            misfit = np.min(all_misfit[chain_iter,:])
            model_index = np.argmin(all_misfit[chain_iter,:])

            
            if misfit < EPS:
                continue 
            
            if misfit < thebestmisfit and misfit > 1e-6: 
                thebestmisfit = misfit
                thebestchain = chain_iter
                thebestindex = model_index
            
            self.find_info_in_chain(chain_iter, "model", model_index)



        # plot ref_model()
        if self.refmodel is not None:
            self.depth , self.vs = self.unpack_thk_vs(self.refmodel)
            self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)
            ax.step(self.vs, self.depth , color='white', ls='-', lw=3, alpha=1, label = 'real model')


            
        # plot init_model()
        '''
        if self.initmodel is not None:
            self.depth , self.vs = self.unpack_thk_vs(self.initmodel)
            self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)
            ax.step(self.vs, self.depth , color='y', ls='-', lw=3, alpha=1, label = 'init model')        
        '''
        '''
        #plot best model 
        self.find_info_in_chain(thebestchain,"model",thebestindex)
        self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)
        ax.step(self.vs, self.depth , color='b', ls='-', lw=3, alpha=1, label = 'best model')
        '''
        
        print("best chains {}".format(thebestchain))
        self.find_info_in_chain(thebestchain,"model_mean",thebestindex)
        self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)
        ax.step(self.vs, self.depth , color='black', ls='-', lw=3, alpha=1, label = 'mean model')

        '''
        self.find_info_in_chain(thebestchain,"model_init",thebestindex)
        print("junliu debug: {}th chain".format(thebestchain))
        self.depth, self.vs = self.interp(self.depth, self.vs, depth_vector)
        ax.step(self.vs, self.depth , color='y', ls='-', lw=3, alpha=1, label = 'init model')
        '''
        
        #print(plt.rcParams['font.size'])
        plt.rcParams['font.size'] = 12.5
        
        ax.invert_yaxis()
        ax.set_xlabel('$V_S$ in km/s',fontdict = {'size':15})
        ax.set_ylabel('Depth in km',fontdict = {'size':15})
        ax.legend()
        #ax.set_title('Best fit models from %d chains' %max_chain)
        plt.xticks(size = 15)
        plt.yticks(size = 15)

        
        return fig

    def interp(self, x_before, y_before, x_after):

       # print(y_before)
        y_after = np.zeros((len(x_after)))

        x_before_index = 0

        for i in range(len(x_after)):
            x = x_before[x_before_index]
            if x_after[i] <= x:
                y_after[i] = y_before[x_before_index]
            else:
                if x_before_index >= len(x_before) - 1:
                    y_after[i] = y_before[x_before_index]
                else:
                    x_before_index += 1
                    y_after[i] = y_before[x_before_index]

        return x_after, y_after
