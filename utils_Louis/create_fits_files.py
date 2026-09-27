import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, text
from pathlib import Path
from shutil import copy
import subprocess
import os
from astropy.io import fits
import h5py
import scipy.interpolate

def read_fits(filename):
    with fits.open(filename, memmap=True) as hdul:
        T = hdul[1].data
        P = hdul[2].data
        wl = hdul[3].data * 1e4

        # Shape: (20 g-points, 5119 wavelengths)
        k = hdul[0].data[:, :, :, :].copy()
        print('test', k.shape)
        #let's make k bigger along the pressure/0th dimesion(for now with 0's)
        k = np.concatenate((k, np.zeros((4, k.shape[1], k.shape[2], k.shape[3]))), axis=0)
        print('test', k.shape)
    
    return T, P, wl, k

def read_hdf5(filename):
    f = h5py.File(filename,'r+') 
    P = 10**f['P'][:] * 1e-5
    T = f['T'][:]
    wl = f['wave'][:] * 1e6
    xsec = 10**f['cross_sec'][:] * 1e4
    
    return T, P, wl, xsec

def change_wl(wl_fits, wl_hdf5, xsec_hdf5, P_hdf5, P_target):

    Pmax = max(P_hdf5)
    if Pmax < 1000:
        print('Warning: We do not have data until 1000 bar (happens at T <= 200K)')
        return None
    else:
        idx = np.argmin(np.abs(P_hdf5 - P_target))
        print('P_target:', P_target, 'P_hdf5[idx]:', P_hdf5[idx])

    interp = scipy.interpolate.interp1d(wl_hdf5, xsec_hdf5[:, idx, 0], kind='linear', bounds_error=False, fill_value=1e-240) 
    xsec = interp(wl_fits)
    return xsec


def write_fits(filename, T, P, wl, k):
    nlam=len(wl)
    nP=len(P)
    nT=len(T)
    ng=20

    lmin = wl[0]*wl[0]/np.sqrt(wl[0]*wl[1])
    lmax =wl[nlam-1]*wl[nlam-1]/np.sqrt(wl[nlam-2]*wl[nlam-1])

    hdr = fits.Header()
    hdr['TMIN']=T[0]
    hdr['TMAX']=T[nT-1]
    hdr['PMIN']=P[0]
    hdr['PMAX']=P[nP-1]
    hdr['L_MIN']=lmin
    hdr['L_MAX']=lmax
    hdr['NT']=nT
    hdr['NP']=nP
    hdr['NLAM']=nlam
    hdr['NG']=ng

    

    primary_hdu = fits.PrimaryHDU(k,header=hdr)
    image_hdu = fits.ImageHDU(T)
    image_hdu2 = fits.ImageHDU(P)
    image_hdu3 = fits.ImageHDU(wl)

    hdul = fits.HDUList([primary_hdu, image_hdu, image_hdu2,image_hdu3])

    hdul.writeto(filename)


if __name__ == "__main__":
    file_name = "/Users/louissiebenaler/ARCiS/Data/Opacities/opacity_K.fits"
    T_fits, P_fits, wl_fits, k_fits = read_fits(file_name)

    P_values = [300, 500, 700, 1000]

    for T in T_fits:
        T = int(T)
        file_name = f'/Volumes/L7aler_HD/PhD/Snellius/cross_sections/K_Allard_full/combined/cross_T{T}.0.hdf5'
        T_hdf5, P_hdf5, wl_hdf5, xsec_hdf5 = read_hdf5(file_name)

        idx = np.argmin(np.abs(T_fits - T))

        for idxP, P in enumerate(P_values):

            xsec = change_wl(wl_fits, wl_hdf5, xsec_hdf5, P_hdf5, P)
            if xsec is None:
                xsec = k_fits[-5, idx, :, :] #just set it to the cross-section of 100 bar that already exists in Arcis
                k_fits[idxP - 4, idx, :, :] = xsec
            else:
                #every g-point should have the same cross-section, so we just take the first g-point (0) and set it to all g-points
                k_fits[idxP - 4, idx, :, :] = xsec
            
    
    
    P_fits = np.concatenate((P_fits, np.array(P_values)), axis=0)
    write_fits("K_test.fits", T_fits, P_fits, wl_fits * 1e-4, k_fits)
    pass