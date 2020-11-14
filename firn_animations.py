#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:10:49 2020

@author: ufbu
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig = plt.figure()


plt.close
def temperature_animation(T_ice_dataframe, Depth_dataframe, resize=100):
    ims_T=[]
    for i in range(int(len(T_ice_dataframe)/resize)):
        j=i*resize
        f_T= plt.plot(T_ice_dataframe.iloc[j], Depth_dataframe.iloc[j], color='blue')
        plt.xlabel(r'Temperature ($T$)')
        plt.ylabel(r'Depth ($m$)')
        plt.xlim(0,300)
        plt.ylim(-10, 80) 
        plt.gca().invert_yaxis()
        #plt.annotate('year='+ str(j/8765.813), (50,255))
        #plt.show(f)
        #im = plt.imshow(f, animated=True)
        ims_T.append(f_T)
    return animation.ArtistAnimation(fig, ims_T, interval=25, blit=True,
                                   repeat_delay=1000)
    
    
def density_animation(rho_dataframe, Depth_dataframe, resize=100):    
    ims_rho=[]
    resize=100
    for i in range(int(len(rho_dataframe)/resize)):
        j=i*resize
        f= plt.plot(rho_dataframe.iloc[j], Depth_dataframe.iloc[j], color='blue')
        plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
        plt.ylabel(r'Depth ($m$)')
        plt.xlim(0,1000)
        plt.ylim(-10, 80)        
        plt.gca().invert_yaxis()
        #plt.show(f)
        #im = plt.imshow(f, animated=True)
        ims_rho.append(f)
    return animation.ArtistAnimation(fig, ims_rho, interval=25, blit=True,
                                   repeat_delay=1000)
# =============================================================================
# 
# 
# def Hist_animation(Depth_dataframe, resize=100):    
#     ims_rho=[]
#     resize=100
#     for i in range(int(len(Depth_dataframe)/resize)):
#         j=i*resize
#         f= plt.hist(Depth_dataframe.iloc[j], bins='auto', alpha=0.7, rwidth=0.5, color='blue')
#         plt.xlabel(r'Depth ($m$)')
#         #plt.show(f)
#         #im = plt.imshow(f, animated=True)
#         ims_rho.append(f)
#     return animation.ArtistAnimation(fig, ims_rho, interval=25, blit=True,
#                                    repeat_delay=1000)
# =============================================================================
