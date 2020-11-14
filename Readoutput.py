#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 13:56:17 2020

@author: ub
"""

import pandas as pd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import firn_animations as fa
import pathlib
import datetime


number_of_layers_int = 200
start_date = datetime.date(1989,1,1)

#Define dictionary with all the information needed to specify which animation 
#was analysed
simulation_inf = {'matlab_output_folder': 'KAN_U_1989-2019_200layers', 'initial_density_profile': 'DensityProfile_KAN-U_1989_metres', 'number_of_layers': '200layers', 'time_period': '1989-2019'}


outputfolder_matlabsim = simulation_inf['matlab_output_folder']
initial_density_profile=simulation_inf['initial_density_profile']
number_of_layers=simulation_inf['number_of_layers']
time_period=simulation_inf['time_period']


target_hr = 205631
target_hr_delta = datetime.timedelta(hours=target_hr)

#Define path to where the output will be placed
path_to_outputs='python_outputs/'+initial_density_profile+'/'+time_period+'/'+number_of_layers+'/'


#Create folder if they dont exstis for the plots and the animations
pathlib.Path(path_to_outputs+'plots/').mkdir(parents=True, exist_ok=True)
pathlib.Path(path_to_outputs+'Animations/').mkdir(parents=True, exist_ok=True)



#Save the simulation information (dictionary) as text, for easy change of analysis later on
f = open(path_to_outputs+"simulation_inf.txt","w")
f.write( str(simulation_inf) )
f.close()


#Adds this name to to output files (Animations and plots):
name='Kan_U_'+initial_density_profile+'_'+time_period+'_'+number_of_layers


#Load dataset rho_bin
fn_rho = 'Output/'+outputfolder_matlabsim+'/rho_bin_1.nc'
ds_rho = nc.Dataset(fn_rho)

#Put rho into a dataframe
rho = ds_rho['rho'][:]
rho_dataframe = pd.DataFrame(rho)




#Load dataset snowc
fn_snowc = 'Output/'+outputfolder_matlabsim+'/snowc_bin_1.nc'
ds_snowc = nc.Dataset(fn_snowc)


#Put snowc into a dataframe
snowc = ds_snowc['snowc'][:]
snowc_dataframe = pd.DataFrame(snowc)



#Load dataset snic
fn_snic = 'Output/'+outputfolder_matlabsim+'/snic_bin_1.nc'
ds_snic = nc.Dataset(fn_snic)


#Put snowc into a dataframe
snic = ds_snic['snic'][:]
snic_dataframe = pd.DataFrame(snic)

rho_all = (snowc_dataframe + snic_dataframe)/ (snowc_dataframe/rho_dataframe + snic_dataframe/917)


#%%
#Load dataset T_ice_bin_1
fn_T_ice = 'Output/'+outputfolder_matlabsim+'/T_ice_bin_1.nc'
ds_T_ice = nc.Dataset(fn_T_ice)

#Put T_ice into a dataframe
T_ice = ds_T_ice['T_ice'][:]
T_ice_dataframe = pd.DataFrame(T_ice)


#Put Depth into a dataframe
Depth = ds_rho['Depth'][:]
Depth_dataframe = pd.DataFrame(Depth)


#Load initial_density profile from data
init_density = pd.read_csv('Input/Initial state/density/'+initial_density_profile+'.csv', sep=',', header='infer')


#Load 2012_density profile from data
density_2012 = pd.read_csv('Input/Initial state/density/RetMIP_density_KAN-U_2012(deleted_NaN).csv', sep=',', header=None)


#%%



ani_rho = fa.density_animation(rho_all, Depth_dataframe)
ani_rho.save(path_to_outputs+'Animations/Densityprofile_animation'+name+'.mp4')

plt.clf() 


ani_T = fa.temperature_animation(T_ice_dataframe, Depth_dataframe)
ani_T.save(path_to_outputs+'Animations/Temperature_animation'+name+'.mp4')

plt.clf() 




#%%
# =============================================================================
# ani_T = fa.Hist_animation(Depth_dataframe)
# ani_T.save(path_to_outputs+'Animations/Depth_histogram_animation'+name+'.mp4')
# 
# plt.clf() 
# 
# =============================================================================
#%%


#Plot initial density_profile from Data
plt.plot(init_density.iloc[:,1], init_density.iloc[:,0])
plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
plt.ylabel(r'Depth ($m$)')
plt.gca().invert_yaxis()
plt.savefig(path_to_outputs+'plots/Densityprofile1989_data'+name+'.pdf')
plt.close()

#%%
#Plot initial density_profile from Simulation
plt.plot(rho_all.iloc[0], Depth_dataframe.iloc[0])
plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
plt.ylabel(r'Depth ($m$)')
plt.gca().invert_yaxis()
plt.savefig(path_to_outputs+'plots/Densityprofile1989_simulation'+name+'.pdf')
plt.close()

#%%
#Plot initial density_profile from Simulation and Data sumoultanously
plt.plot(init_density.iloc[:,1], init_density.iloc[:,0])
plt.plot(rho_all.iloc[0], Depth_dataframe.iloc[0])
plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
plt.ylabel(r'Depth ($m$)')
plt.gca().invert_yaxis()
plt.savefig(path_to_outputs+'plots/Densityprofile1989_simulation_and_data'+name+'.pdf')
plt.close()
#%%






if time_period=='1989-2019':
    #Plot 2012 density_profile from Data
    plt.plot(density_2012.iloc[:,1], density_2012.iloc[:,0])
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    plt.savefig(path_to_outputs+'plots/Densityprofile2012_data'+name+'.pdf')
    plt.close()


    #PLot density profile mid-May 2012
    plt.plot(rho_all.iloc[target_hr], Depth_dataframe.iloc[target_hr])
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    plt.savefig(path_to_outputs+'plots/density_profile_may_2012_simulation'+name+'.pdf')
    plt.close()


    #Plot initial density_profile from Simulation and Data sumoultanously
    plt.plot(density_2012[1], density_2012[0])
    plt.plot(rho_all.iloc[target_hr], Depth_dataframe.iloc[target_hr])
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    plt.text(400, 30, str(start_date+ target_hr_delta))
    plt.ax.legend([ ])
    plt.savefig(path_to_outputs+'plots/Densityprofile2012_simulation_and_data'+name+'.pdf')
    plt.close()

    
        
    plt.plot(density_2012[1], density_2012[0])
    plt.plot(rho_all.iloc[target_hr], Depth_dataframe.iloc[target_hr])
    plt.plot(diff_dens, density_2012[0])
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    plt.savefig(path_to_outputs+'plots/Difference_in_Densityprofile2012_simulation_and_data'+name+'.pdf')
    plt.close()
    
    
    
    
    
    plt.plot(np.gradient(density_2012[1]), density_2012[0])
    plt.plot(np.gradient(rho_all.iloc[target_hr]), Depth_dataframe.iloc[target_hr])
    #plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    plt.savefig(path_to_outputs+'plots/Gradients_of_Densityprofile2012_simulation_and_data'+name+'.pdf')
    plt.close()
    



    

