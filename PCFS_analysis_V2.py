# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 10:47:29 2024

@author: colbu
"""

import numpy as np
import time as timing
import os, struct, scipy, re, glob
from Editied_Photons_for_any_TTBIN import Photons 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from Edited_PCFS import PCFS


#%%

pcfs2 = PCFS(r'C:\Data and Code\Data\24_10_10_SFRS_after_artifact\SFRS1', 2,simulation = False)

pcfs2.get_photons_all()

pcfs2.get_sum_signal_all()

#%%
pcfs2.get_intensity_correlations((1e3,1e13), 3)
#%%

pcfs2.get_blinking_corrected_PCFS()
#%%
white_fringe = 19.245790
pcfs2.plot_spectral_diffusion([1e6,1e9], 1)
print(19.245790-pcfs2.stage_positions)

    #%%

pcfs2.get_splev_mirror_spec_corr( 19.245790, 0,fit_interferogram = False)

aaaa = pcfs2.PCFS_interferogram

crossG2 = pcfs2.cross_corr_interferogram 
autoG2 = pcfs2.auto_corr_sum_interferogram

for i in range(50):
   if i%1 == 0:
    #plt.plot(pcfs2.tau,crossG2[:,i],label ='cross' )
    #plt.plot(pcfs2.tau,autoG2[:,i], label ='sum-auto' )
    plt.plot(pcfs2.tau,1-crossG2[:,i]/autoG2[:,i],label = f'stage pos {i+1}')
    plt.xscale('log')
    plt.ylim(-.15,.21)
plt.xlim(1e6,1e12)
#plt.legend()
plt.show()
    
print(pcfs2.mirror_stage_positions)
print(pcfs2.stage_positions)
   #%%
pcfs2.plot_splev_spec_corr([3e5,1e6,9e6,1e9,1e10], (-50,50))
    
from scipy.optimize import curve_fit

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))
def lorenzian(x, x0, gamma, A):
    return A * (gamma / ((x - x0) ** 2 + gamma ** 2) / np.pi)

fwhm_arr2 = []
x = pcfs2.splev_spec_corr['zeta']*8
plt.figure(figsize = (10,6),dpi = 200)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Assuming `x`, `pcfs2.splev_spec_corr['spectral_corr']`, `pcfs2.tau`, `pcfs2.test`, and `pcfs2.mirror_PCFS_interferogram` are already defined
print(pcfs2.tau[40]/1e9)
indices = [30,31,32,33,34,35,36,37,38,39,40]

# Normalize the iteration number for colormap
norm = Normalize(vmin=0, vmax=len(indices) - 1)
cmap = plt.get_cmap('Purples')

# Plot 1
plt.figure(figsize = (10,6), dpi = 200)
for iteration, i in enumerate(indices):
    color = cmap(norm(iteration))
    plt.plot(x, pcfs2.splev_spec_corr['spectral_corr'][i, :] / max(pcfs2.splev_spec_corr['spectral_corr'][i, :]), 
             label=f"{pcfs2.tau[i] / 1e6:.1f}", color=color)
plt.xlim(-40, 40)
plt.ylabel('Norm Spectral Corr')
plt.xlabel(r'$\zeta (cm^{-1})$')
plt.legend()
plt.show()

# Plot 2
plt.figure()
for iteration, i in enumerate(indices):
    color = cmap(norm(iteration))
    plt.plot(pcfs2.test, pcfs2.mirror_PCFS_interferogram[i, :], label=f"{pcfs2.tau[i] / 1e6:.1f}", color=color)
plt.ylabel('interferogram')
plt.xlabel(r'stage pos (mm')
plt.legend()
plt.show()


#%%
for i in range(len(pcfs2.tau)):
    
    
    y = pcfs2.splev_spec_corr['spectral_corr'][i,:]/max(pcfs2.splev_spec_corr['spectral_corr'][i,:])
    x = pcfs2.splev_spec_corr['zeta']
    
    try:
        params, covariance = curve_fit(lorenzian, x, y)
    except RuntimeError:
        fwhm_arr2.append(0)
    except ValueError:
        fwhm_arr2.append(0)
    else:
        #for lorenzian
        fwhm = 2 * params[1] 
        #for gaussian
        #fwhm = 2 * np.sqrt(2 * np.log(2)) * params[2]
        fwhm_arr2.append(fwhm)
        print(fwhm)
        #plt.plot(x, y, 'bo', label='Data')
        #plt.plot(x, lorenzian(x, *params), 'r-', label='Fit')
        #plt.legend()
       # plt.show()
#%%
fwhm_arr_cm2 = np.zeros(len(fwhm_arr2))
for i in range(len(fwhm_arr2)):
    fwhm_arr_cm2[i] = fwhm_arr2[i]*8
plt.figure(dpi = 200)
plt.scatter(pcfs2.tau,fwhm_arr_cm2, label='Scatter Plot on Log Scales')
plt.xscale('log')
plt.xlim(1e4,1e11)
plt.ylim(0,160)
plt.ylabel(f'fwhm $(cm^-1)$')
plt.xlabel('tau')
plt.show()
file_path = r'C:\Data and Code\Data\24_02_11_ensamble_PCFS_SIMS\FWHMensambleSize1.npy'

# Save the array to the specified file
#np.save(file_path, fwhm_arr)


#%%
 
    pcfs3 = PCFS(r"C:\Data and Code\Data\24_07_10_SFRS_SI\Run1", 2,simulation = False)
    
    pcfs3.get_photons_all()
    
    pcfs3.get_sum_signal_all()
    
    pcfs3.get_intensity_correlations((1e1,1e11), 3)
    
    pcfs3.get_blinking_corrected_PCFS()
    
    pcfs3.plot_spectral_diffusion([1e6,1e9], -4)
    
    pcfs3.get_splev_mirror_spec_corr( -38.802476, 0,fit_interferogram = True)
    
    pcfs3.plot_splev_spec_corr([3e5,1e6,1e7,1e9,1e10], (-2.5,2.50))
    
    from scipy.optimize import curve_fit
    
    fwhm_arr3 = []
    for i in range(len(pcfs3.tau)):
        
        
        y = pcfs3.splev_spec_corr['spectral_corr'][i,:]/max(pcfs3.splev_spec_corr['spectral_corr'][i,:])
        x = pcfs3.splev_spec_corr['zeta']
        
        try:
            params, covariance = curve_fit(lorenzian, x, y)
        except RuntimeError:
            fwhm_arr3.append(0)
        except ValueError:
            fwhm_arr3.append(0)
        else:
            #for lorenzian
            fwhm = 2 * params[1] 
            #for gaussian
            #fwhm = 2 * np.sqrt(2 * np.log(2)) * params[2]
            fwhm_arr3.append(fwhm)
            #plt.plot(x, y, 'bo', label='Data')
            #plt.plot(x, lorenzian(x, *params), 'r-', label='Fit')
            #plt.legend()
            #plt.show()
    
    fwhm_arr_cm3 = np.zeros(len(fwhm_arr3))
    for i in range(len(fwhm_arr3)):
        fwhm_arr_cm3[i] = fwhm_arr3[i]*8
    plt.scatter(pcfs3.tau,fwhm_arr_cm3, label='Scatter Plot on Log Scales')
    plt.xscale('log')
    plt.xlim(1e4,1e11)
    plt.ylim(0,8)
    plt.show()
    file_path = r'C:\Data and Code\Data\24_02_11_ensamble_PCFS_SIMS\FWHMensambleSize1.npy'
    
    # Save the array to the specified file
    #np.save(file_path, fwhm_arr)


#%%
if 1:  
    pcfs4 = PCFS(r"C:\Data and Code\Data\24_07_10_SFRS_SI\Run1", 2,simulation = False)
    
    pcfs4.get_photons_all()
    
    pcfs4.get_sum_signal_all()
    
    pcfs4.get_intensity_correlations((1e1,1e11), 3)
    
    pcfs4.get_blinking_corrected_PCFS()
    
    pcfs4.plot_spectral_diffusion([1e6,1e9], -4)
    
    pcfs4.get_splev_mirror_spec_corr( -38.802476, 0,fit_interferogram = True)
    
    pcfs4.plot_splev_spec_corr([3e5,1e6,1e7,1e9,1e10], (-2.5,2.50))
    
    from scipy.optimize import curve_fit
    
   
    fwhm_arr4 = []
    for i in range(len(pcfs4.tau)):
        
        
        y = pcfs4.splev_spec_corr['spectral_corr'][i,:]/max(pcfs4.splev_spec_corr['spectral_corr'][i,:])
        x = pcfs4.splev_spec_corr['zeta']
        
        try:
            params, covariance = curve_fit(lorenzian, x, y)
        except RuntimeError:
            fwhm_arr4.append(0)
        except ValueError:
            fwhm_arr4.append(0)
        else:
            #for lorenzian
            fwhm = 2 * params[1] 
            #for gaussian
            #fwhm = 2 * np.sqrt(2 * np.log(2)) * params[2]
            fwhm_arr4.append(fwhm)
            #plt.plot(x, y, 'bo', label='Data')
            #plt.plot(x, lorenzian(x, *params), 'r-', label='Fit')
            #plt.legend()
            #plt.show()
    
    fwhm_arr_cm4 = np.zeros(len(fwhm_arr4))
    for i in range(len(fwhm_arr4)):
        fwhm_arr_cm4[i] = fwhm_arr4[i]*8
    plt.scatter(pcfs4.tau,fwhm_arr_cm4, label='Scatter Plot on Log Scales')
    plt.xscale('log')
    plt.xlim(1e4,1e11)
    plt.ylim(0,8)
    plt.show()
    file_path = r'C:\Data and Code\Data\24_02_11_ensamble_PCFS_SIMS\FWHMensambleSize1.npy'
    
    # Save the array to the specified file
    #np.save(file_path, fwhm_arr)



#%%


# Create a scatter plot with logarithmic scales and open circles
plt.figure(figsize=(10, 6), dpi=600)  # Adjust the figure size and DPI as needed

scatter = plt.scatter(pcfs1.tau, fwhm_arr, c='blue', marker='o', edgecolors='none', label='FWHM vs. Tau')
plt.xscale('log')
plt.xlim(1e4, 1e11)
plt.ylim(0.5, 3)

# Labeling and titles with larger font sizes
plt.xlabel('Tau', fontsize=16,  color='black')
plt.ylabel(rf'FWHM of $P(\zeta)$', fontsize=16, color='black')
plt.title('Spectral Diffusion', fontsize=18, fontweight='bold',color='black')

# Remove gridlines
plt.grid(False)

# Customize ticks and tick labels
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=12)

# Legend
plt.legend()

# Increase the marker size
scatter.set_sizes([50])
custom_dash = [10, 7]
# Connect the dots with very thin lines
plt.plot(pcfs1.tau, fwhm_arr, linestyle=(0, (custom_dash[0], custom_dash[1])), color='black', linewidth=0.5)
plt.savefig(r"C:\Data and Code\Data\23_10_06_4ATPsimulations\Plots\my_plot.png", dpi=300)
plt.show()


#%%





# Create instances of PCFS and store them in a list
pcfs_list = [
    PCFS(r"C:\Data and Code\Data\24_07_10_SFRS_SI\Run1", 2, simulation=False),
    PCFS(r"C:\Data and Code\Data\24_07_10_SFRS_SI\Run2", 2, simulation=False),
    PCFS(r"C:\Data and Code\Data\24_07_08_Si_raman_PCFS\run2", 2, simulation=False),
    PCFS(r"C:\Data and Code\Data\24_07_08_Si_raman_PCFS\Run1", 2, simulation=False)
]
splev_params = [
    -38.802476,
    -38.902476,
    -38.808600,
    -38.808600
]

# Loop over each instance and perform the analysis
for i,pcfs in enumerate(pcfs_list):
    pcfs.get_photons_all()

    pcfs.get_sum_signal_all()

    pcfs.get_intensity_correlations((1e1, 1e11), 3)

    pcfs.get_blinking_corrected_PCFS()

    pcfs.plot_spectral_diffusion([1e6, 1e9], -4)

    pcfs.get_splev_mirror_spec_corr(splev_params[i], 0, fit_interferogram=True)


    pcfs.plot_splev_spec_corr([3e5, 1e6, 1e7, 1e9, 1e10], (-2.5, 2.50))
