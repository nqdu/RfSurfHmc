import matplotlib.pyplot as plt 
import numpy as np 
import os
import os.path as op
from scipy.interpolate import interp1d

def next_power_of_2(x):
    return 1 if x == 0 else 2**(x-1).bit_length()

def get_rf_inv_para(real_rf_dobs, real_rf_time, t_start, t_end):
    if t_start < real_rf_time[0] or t_end > real_rf_time[-1]:
        raise Exception("wrong t_start or t_end")

    assert len(real_rf_dobs) == len(real_rf_time)
    obs_dt = real_rf_time[1] - real_rf_time[0]

    t_syn_len = t_end - t_start
    nt = int(t_syn_len / obs_dt)
    nt = next_power_of_2(nt)

    t = np.linspace(t_start, t_end, nt)
    
    dt = t[1] - t[0]
    time_shift = -t_start

    fun1 = interp1d(real_rf_time, real_rf_dobs, kind = 'cubic')
    interp_dobs = fun1(t)


    return interp_dobs, nt, dt, time_shift