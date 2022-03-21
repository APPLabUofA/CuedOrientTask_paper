# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 13:38:10 2021

@author: ssshe
"""


# General imports
import numpy as np
# from scipy.io import loadmat, savemat
import scipy.io as sio
# from os.path import dirname, join as pjoin
import os
import matplotlib.pyplot as plt

# Import the FOOOF object, and custom helper & utility functions
from fooof import FOOOF, FOOOFGroup
from fooof.analysis import get_band_peak_fm, get_band_peak_fg
from fooof.utils import trim_spectrum
from fooof.data import FOOOFSettings
# Import the Bands object, which is used to define frequency bands
from fooof.bands import Bands
# Import plotting functions
from fooof.plts.spectra import plot_spectrum, plot_spectra
from fooof.plts.spectra import plot_spectrum_shading, plot_spectra_shading

# Import spectral power functions
from neurodsp.spectral import compute_spectrum
from neurodsp.plts.spectral import plot_power_spectra



# Change directory
os.chdir("C:\\Users\\ssshe\\Documents\\MathLab\\Analysis\\OrientWheel_Exo")

# Name of processing settings
set_name = 'byTargets_v3'

# Location of saved data files
# data_dir = os.path.join(os.getcwd(), 'fooof_trl','fix', set_name) #single trial
# data_dir = os.path.join(os.getcwd(), 'fooof_trl','mid', set_name) #single trial
# data_dir = os.path.join(os.getcwd(), 'fooof_trl','late', set_name) #single trial
data_dir = os.path.join(os.getcwd(), 'fooof_trl','post', set_name) #single trial


# Location to save fooof results
# save_dir = os.path.join(os.getcwd(), 'fooof_trl','fix', set_name, 'fooof_results') #single trial
# save_dir = os.path.join(os.getcwd(), 'fooof_trl','mid', set_name, 'fooof_results') #single trial
# save_dir = os.path.join(os.getcwd(), 'fooof_trl','late', set_name, 'fooof_results') #single trial
save_dir = os.path.join(os.getcwd(), 'fooof_trl','post', set_name, 'fooof_results') #single trial


# Load the participant names
partname = sio.loadmat(os.path.join(data_dir, 'participants.mat'))
partlist = partname['nametmp']
# partlist = partname['nametmp'][2:3] #selecting a couple participants

# Define frequency bands of interest
bands = Bands({'delta' : [1, 3.5],
               'theta' : [3.5, 7.5],
               'alpha' : [7.5, 14.15],
               'beta1' : [14.5, 22.5]})

# List electrode labels
elect = ['Oz','Pz','Cz','FCz','Fz','O1','O2','PO3','PO4','P7','P8','P5','P6',
         'P3','P4','CP5','CP6','CP1','CP2','C3','C4','FC5','FC6','FC1','FC2',
         'F7','F8','F3','F4','Fp1','Fp2']

# List of condition labels
conds = ['L_iv','L_v','R_iv','R_v']

fs = 1000 #sampling rate


for i, part in enumerate(partlist[:]):
    
    for ee, ichan in enumerate(elect):
        
        for cc, icond in enumerate(conds):
            # Get list of data files in folder
            mat_fname = os.listdir(os.path.join(data_dir, part, ichan, icond))
            
            # Make new folder is one does not exist
            if not(os.path.isdir(os.path.join(save_dir, part))):
                os.mkdir(os.path.join(save_dir, part))
            # Make new folder is one does not exist    
            if not(os.path.isdir(os.path.join(save_dir, part, ichan))):
                os.mkdir(os.path.join(save_dir, part, ichan))  
            # Make new folder is one does not exist    
            if not(os.path.isdir(os.path.join(save_dir, part, ichan, icond))):
                os.mkdir(os.path.join(save_dir, part, ichan, icond))    
   
            
            for j, ifile in enumerate(mat_fname):    
                # Load the mat file 
                data = sio.loadmat(os.path.join(data_dir, part, ichan, icond, ifile))
                
                sig = np.squeeze(data['sig']).astype('float')
                sig = sig.T # Transpose signal
        
                # Median of spectrogram ("median Welch")
                freqs, psds = compute_spectrum(sig, fs, method='welch', avg_type='mean', nperseg=sig.shape[0])
        
                # Unpack data from dictionary, and squeeze numpy arrays
                # freqs = np.squeeze(data['freqs']).astype('float')
                # psds = np.squeeze(data['psd']).astype('float')
                # # ^Note: this also explicitly enforces type as float (type casts to float64, instead of float32)
                # #  This is not strictly necessary for fitting, but is for saving out as json from FOOOF, if you want to do that
        
                # # Transpose power spectra, to have the expected orientation for FOOOF
                # psds = psds.T
                
                # Initialize FOOOFGroup object
                fm = FOOOF(max_n_peaks=8,peak_width_limits=[1,5])
                # fg = FOOOFGroup(max_n_peaks=8,peak_threshold=0.2,min_peak_height=0.22,peak_width_limits=[1,5])
                
                # Run FOOOF across all power spectra
                fm.fit(freqs, psds, [1, 40])
                
                # Fit the FOOOF model on all PSDs, and report
                # fm.report(freqs, psds, [1, 40])
                
                # Get all band peaks from a group of power spectrum models
                beta1 = get_band_peak_fm(fm, bands.beta1)
                alphas = get_band_peak_fm(fm, bands.alpha)
                thetas = get_band_peak_fm(fm, bands.theta)
                delta = get_band_peak_fm(fm, bands.delta)
                
                
                # Save out fooof results to json file
                #  There is a utility file to load this json file directly into Matlab
                savename = 'f_results_' + ifile[6:-4]
                fm.save(os.path.join(save_dir, part, ichan, icond, savename), 
                        save_results=True)
                
                # Save out a specific FOOOF measure of interest
                savenameb = 'bands_' + ifile[6:-4] + '.mat'
                sio.savemat(os.path.join(save_dir, part, ichan, icond, savenameb), 
                            {'beta1' : beta1, 'alphas' : alphas, 'thetas' : thetas,
                             'delta' : delta})
                
                # Save out a specific FOOOF measure of interest
                # spect_flat = fm._spectrum_flat
                # spect_peakrm = fm._spectrum_peak_rm
                # spect_freqs = fm.freqs
                # savenamec = 'spect_new_' + ifile[6:-4] + '.mat'
                # sio.savemat(os.path.join(save_dir, part, ichan, icond, savenamec), 
                #             {'spect_flat' : spect_flat, 'spect_peakrm' : spect_peakrm,
                #              'spect_freqs' : spect_freqs})
                
                # Clear variables
                del data, freqs, psds, savename, savenameb, ifile, fm, alphas, thetas, beta1, delta,
                sig
            del mat_fname, j 
        del cc, icond
    del ichan, ee
del i       
        



# Plot the power spectra
plot_power_spectra([freqs[:200]],
                   [psds[:200,0]],
                   ['Median Welch'])



# Find the index of the worst model fit from the group
worst_fit_ind = np.argmax(fg.get_params('error'))
# Extract this model fit from the group
fm = fg.get_fooof(worst_fit_ind, regenerate=True)     

# Check results and visualize the extracted model
fm.plot() 
fm.print_results()
