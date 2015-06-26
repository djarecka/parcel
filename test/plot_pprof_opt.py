import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from libcloudphxx import common
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.io import netcdf
from parcel import parcel
import numpy as np

pprof_list = ["pprof_const_rhod", "pprof_const_th_rv", "pprof_piecewise_const_rhod"]

def test_pressure(dt=1):
    # running parcel model for different ways to solve for pressure  ...
    for pprof in pprof_list:
        parcel(dt=dt, outfreq = 10, pprof = pprof, outfile="test_" + pprof + ".nc")

    # ... plotting the results ...
    plt.figure(1, figsize=(18,10))
    plots    = []
    legend_l = []

    for i in range(6):
        plots.append(plt.subplot(2,3,i+1))

    plots[0].set_xlabel('p [hPa]')

    plots[1].ticklabel_format(useOffset=False) 
    plots[1].set_xlabel('th_d [K]')
    plots[2].set_xlabel('T [K]')
    # the different ways of solving for pressure come from different assumptions about the density profile
    # but those assumptions are just used when calculating the profile of pressure
    # later on the rho_d profile can be calculated (and is not the same as the one assumed)
    # so the kappa here is the actual profile of rho_d during the simulation (different than the one assumed)
    plots[3].set_xlabel('kappa(rho_d :)) [kg/m3]')  
    plots[4].set_xlabel('rv [g/kg]')
    plots[5].set_xlabel('RH')

    for ax in plots:
        ax.set_ylabel('z [m]')

    style = ["g.-", "b.-","r.-"]
    for i, pprof_val in enumerate(pprof_list):
        f = netcdf.netcdf_file("test_"+pprof_val+".nc", "r")
        z = f.variables["z"][:]
        plots[0].plot(f.variables["p"][:] / 100.   , z, style[i])
        plots[1].plot(f.variables["th_d"][:]       , z, style[i])
        plots[2].plot(f.variables["T"][:]          , z, style[i])
        plots[3].plot(f.variables["rhod"][:]       , z, style[i])
        plots[4].plot(f.variables["r_v"][:] * 1000 , z, style[i])
        plots[5].plot(
	  f.variables["RH"][:]                     , z, style[i], 
	  [f.variables["RH"][:].max()] * z.shape[0], z, style[i]
        )
        legend_l.append(pprof_val)
    plots[0].legend(legend_l, loc=1, prop = FontProperties(size=10))
    plt.savefig("plot_pressure.svg")

