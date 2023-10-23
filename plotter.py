import numpy as np 
import timeit
import matplotlib.pyplot as plt
from pathlib import Path
# import blobs_solver2 as pHyFlow
import pHyFlow
# import VorticityFoamPy as foam
import solvers.particle as particle

import os
import sys
import yaml
import re
import csv
import pandas

#---------------Current directory and paths-------------------
arg = sys.argv
if len(arg) > 2:
    raise Exception("More than two arguments inserted!")
if len(arg) <= 1:
    raise Exception("No config file specificed!")
configFile = arg[1]

#-----------------------Config the yaml file ------------------
loader = yaml.SafeLoader
loader.add_implicit_resolver(
    u'tag:yaml.org,2002:float',
    re.compile(u'''^(?:
     [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
    |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
    |\\.[0-9_]+(?:[eE][-+][0-9]+)?
    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
    |[-+]?\\.(?:inf|Inf|INF)
    |\\.(?:nan|NaN|NAN))$''', re.X),
    list(u'-+0123456789.'))
config = yaml.load(open(os.path.join(configFile)),Loader=loader)

start_index = configFile.index('_') + 1
end_index = configFile.index('.')

case = config['case']

case_dir = os.path.join(os.getcwd(), 'results', case)
data_dir = os.path.join(case_dir,config["data_folder"])
plots_dir = os.path.join(case_dir,config["plots_folder"])

Path(plots_dir).mkdir(parents=True, exist_ok=True)

nTimeSteps = config['nTimeSteps']
coreSize = config['coreSize']
gammaC = config['gammaC']
deltaTc = config['deltaTc']
writeInterval_plots = config['writeInterval_plots']
run_analytical_flag = config['run_analytical_flag']
N_initial = config['N_initial']


uxNorm = np.array([])
uyNorm = np.array([])
omegaNorm = np.array([])
t_norm = np.array([])
#Line plots
times_file = os.path.join(data_dir,"times_{}.csv".format(case))
times_data = pandas.read_csv(times_file)

time = times_data['Time']
noBlobs = times_data['NoBlobs']
evolution_time = times_data['Evolution_time']
circulation = times_data['Circulation']

fig, ax = plt.subplots(1,1,figsize=(6,6))
ax.plot(time,noBlobs, label='No of Particles')
plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.title('Total number of particles')
plt.ylabel('Particles')
plt.xlabel('time $(sec)$')
plt.legend()
plt.savefig("{}/number_of_particles_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

fig, ax = plt.subplots(1,1,figsize=(6,6))
ax.plot(time,circulation- gammaC, label='Circulation deficit')
plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.title('absolute error in circulation')
plt.ylabel('circulation')
plt.xlabel('time $(sec)$')
plt.legend()
plt.savefig("{}/circulation_error_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

fig = plt.subplots(figsize=(6,6))
index = np.arange(len(evolution_time))
width = deltaTc
lagrangian = plt.bar(index[1:]*deltaTc, evolution_time[1:], width)
plt.ylabel('Time (s)')
plt.xlabel('Simulation time (s)')
plt.title('Evolution time')
plt.savefig("{}/times_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

for timeStep in range(N_initial, nTimeSteps+1+N_initial):
    if timeStep+N_initial%writeInterval_plots == 0:
        ####Fields
        lagrangian_file = os.path.join(data_dir,'results_{}_{n:06d}.csv'.format(case,n=timeStep))
        lagrangian_data = np.genfromtxt(lagrangian_file)

        xplot = lagrangian_data[:,0]
        yplot = lagrangian_data[:,1]
        length = int(np.sqrt(len(xplot)))
        xPlotMesh = xplot.reshape(length,length)
        yPlotMesh = yplot.reshape(length,length)

        lagrangian_ux = lagrangian_data[:,2]
        lagrangian_uy = lagrangian_data[:,3]
        lagrangian_omega = lagrangian_data[:,4]

        analytical_file = os.path.join(data_dir,'results_analytical_{n:06d}.csv'.format(case,n=timeStep))
        analytical_data = np.genfromtxt(analytical_file)

        analytical_ux = analytical_data[:,2]
        analytical_uy = analytical_data[:,3]
        analytical_omega = analytical_data[:,4]

        xTicks = np.linspace(-2,2,5)
        yTicks = np.linspace(-2,2,5)

        fig, ax = plt.subplots(1,1,figsize=(6,6))
        ax.set_aspect("equal")
        ax.set_xticks(xTicks)
        ax.set_yticks(yTicks)
        plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
        plt.minorticks_on()
        plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        cax = ax.contourf(xPlotMesh,yPlotMesh,lagrangian_omega.reshape(length,length),levels=100,cmap='RdBu',extend="both")
        cbar = fig.colorbar(cax,format="%.4f")
        cbar.set_label("Vorticity (1/s)")
        plt.tight_layout()
        plt.savefig("{}/vorticity_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
        plt.close(fig)
        
        fig, ax = plt.subplots(1,1,figsize=(6,6))
        ax.set_aspect("equal")
        ax.set_xticks(xTicks)
        ax.set_yticks(yTicks)
        plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
        plt.minorticks_on()
        plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        cax = ax.contourf(xPlotMesh,yPlotMesh,lagrangian_ux.reshape(length,length),levels=100,cmap='RdBu',extend="both")
        cbar = fig.colorbar(cax,format="%.4f")
        cbar.set_label("Velocity (1/s)")
        plt.tight_layout()
        plt.savefig("{}/velocity_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
        plt.close(fig)

        if run_analytical_flag==True:
            fig, ax = plt.subplots(1,1,figsize=(6,6))
            ax.set_aspect("equal")
            ax.set_xticks(xTicks)
            ax.set_yticks(yTicks)
            plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
            plt.minorticks_on()
            plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            cax = ax.contourf(xPlotMesh,yPlotMesh,analytical_omega.reshape(length,length),levels=100,cmap='RdBu',extend="both")
            cbar = fig.colorbar(cax,format="%.4f")
            cbar.set_label("Vorticity (1/s)")
            plt.tight_layout()
            plt.savefig("{}/vorticity_analytical_{}.png".format(plots_dir,timeStep), dpi=300, bbox_inches="tight")
            plt.close(fig)

            fig, ax = plt.subplots(1,1,figsize=(6,6))
            ax.set_aspect("equal")
            ax.set_xticks(xTicks)
            ax.set_yticks(yTicks)
            plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
            plt.minorticks_on()
            plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            cax = ax.contourf(xPlotMesh,yPlotMesh,analytical_ux.reshape(length,length),levels=100,cmap='RdBu',extend="both")
            cbar = fig.colorbar(cax,format="%.4f")
            cbar.set_label("Velocity (1/s)")
            plt.tight_layout()
            plt.savefig("{}/velocity_analytical_{}.png".format(plots_dir,timeStep), dpi=300, bbox_inches="tight")
            plt.close(fig)

#### Errors
        omegaScale = np.max(np.abs(analytical_omega))

        omega_error = ((lagrangian_omega - analytical_omega)/(omegaScale))*100
        fig, ax = plt.subplots(1,1,figsize=(6,6))
        ax.set_aspect("equal")
        ax.set_xticks(xTicks)
        ax.set_yticks(yTicks)
        plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
        plt.minorticks_on()
        plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        cax = ax.contourf(xPlotMesh,yPlotMesh,omega_error.reshape(length,length),levels=100,cmap='jet',extend="both")
        cbar = fig.colorbar(cax,format="%.4f")
        cbar.set_label("Vorticity error (%) (1/s)")
        plt.tight_layout()
        plt.savefig("{}/vorticity_error_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
        plt.close(fig)


#### Blobs distribution

        blobs_file = os.path.join(data_dir,'blobs_{}_{n:06d}.csv'.format(case,n=timeStep))
        blobs_data = np.genfromtxt(blobs_file)

        if len(blobs_data.shape) == 1:
            blobs_x = blobs_data[0]
            blobs_y = blobs_data[1]
            blobs_g = blobs_data[2]
            if coreSize == 'variable':
                blobs_sigma = blobs_data[3]
        else:
            blobs_x = blobs_data[:,0]
            blobs_y = blobs_data[:,1]
            blobs_g = blobs_data[:,2]
            if coreSize == 'variable':
                blobs_sigma = blobs_data[:,3]
        if coreSize == 'variable':
            fig, ax = plt.subplots(1,1,figsize=(6,6))
            ax.scatter(blobs_x,blobs_y,c=blobs_g, s= blobs_sigma * 0.2 / np.min(blobs_sigma))
            plt.savefig("{}/blobs_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
            plt.close(fig)
        else:
            fig, ax = plt.subplots(1,1,figsize=(6,6))
            ax.scatter(blobs_x,blobs_y,c=blobs_g, s=0.2)
            plt.savefig("{}/blobs_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
            plt.close(fig)


#### L2 errors in vorticity and velocity


        uxNorm = np.append(uxNorm,np.linalg.norm(lagrangian_ux-analytical_ux)/np.linalg.norm(analytical_ux))
        uyNorm = np.append(uyNorm,np.linalg.norm(lagrangian_uy-analytical_uy)/np.linalg.norm(analytical_uy))
        omegaNorm = np.append(omegaNorm,np.linalg.norm(lagrangian_omega-analytical_omega)/np.linalg.norm(analytical_omega))
        t_norm = np.append(t_norm,deltaTc*timeStep)
fig, ax = plt.subplots(1,1,figsize=(6,6))
plt.plot(t_norm,uxNorm, label='ux-L2 error')
plt.plot(t_norm,uyNorm, label='uy-L2 error')
plt.plot(t_norm,omegaNorm, label='omega-L2 error')
plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.ylabel('L2-error')
plt.xlabel('time')
plt.legend()
plt.savefig("{}/L2_error_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")
plt.close(fig)