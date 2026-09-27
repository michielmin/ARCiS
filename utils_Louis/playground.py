import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, text
from pathlib import Path
from shutil import copy
import subprocess
import os
from arcis_wrapper import ArcisWrapper
import pyARCiS




if __name__ == "__main__":

    


    # 'Al2O3[s]', 'TiO2[s]', 'Ti4O7[s]', 'MgTi2O5[s]', 'CaTiO3[s]', 'MgO[s]', 'MgAl2O4[s]', 'FeO[s]', 'Fe2O3[s]', 'Fe3O4[s]', 'Fe[s]', 'Zn[s]', 'Cr[s]', 'Ni[s]', 'W[s]', 'SiO[s]', 'SiO2[s]', 'MgSiO3[s]', 'ENSTATITE[s]', 'Mg2SiO4[s]', 'FORSTERITE[s]', 'FeSiO3[s]', 'FERROSILITE[s]', 'Fe2SiO4[s]', 'FAYALITE[s]', 'NaAlSi3O8[s]', 'CaSiO3[s]', 'ZnS[s]', 'FeS[s]', 'Na2S[s]', 'MnS[s]', 'NaCl[s]', 'KCl[s]', 'H2O[s]', 'NH3[s]', 'NH4SH[s]', 
    '''
    case = ArcisWrapper()
    path = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/WJupiter_syoung_clouds/'
    c = 'black'
    case.plot_tp_profile(path, show = False, color = c, ls = '--')
    print('Teq for clouds:', case.get_Teq(path))

    path = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/WJupiter_syoung_clouds_K10000000000.0/'
    c = 'red'
    case.plot_tp_profile(path, show = True, color = c, ls = '--')
    print('Teq for gas:', case.get_Teq(path))
    '''


    
    case = ArcisWrapper()

    '''
    print(case.get_Teq(output_path = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Dplanet_{1.5}_clouds_correct2/', bond_albedo = 0.022))
    print(case.get_Teq(output_path = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Dplanet_{0.8}_clouds_correct2/', bond_albedo = 0.341))

    #Tint = [25, 100, 200, 400, 600, 1000, 1500]
    D = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 2]
    D = [0.2, 0.5, 1, 2]
    #logg = [1.5, 3, 4, 5]

    colors = ['black', 'red', 'green', 'orange', 'blue', 'purple', 'brown']
    for i, d in enumerate(D):
        path_gas = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Dplanet_{d}_gas/'
        path_clouds = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Dplanet_{d}_clouds/'
        c = colors[i]
        if i < len(D)-1:
            case.plot_tp_profile(path_gas, show = False, color = c, ls = '--', label = f'Gas-only, D={d} au')
            case.plot_tp_profile(path_clouds, show = False, color = c, ls = '-', label = f'Clouds, D={d} au')
        else:
            plt.title(r'Jupiter-mass/radius-composition:$K_{\rm zz}=10^{12}$ cm$^2$/s, $\dot{\Sigma}_{\rm n} = 10^{-11}$, $T_{\rm int} = 200K$', fontsize = 10)
            case.plot_tp_profile(path_gas, show = False, color = c, ls = '--', label = f'Gas-only, D={d} au')
            case.plot_tp_profile(path_clouds, show = True, color = c, ls = '-', label = f'Clouds, D={d} au')

    d = 0.8
    path_clouds = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Dplanet_{d}_clouds_correct2/'
    case.plot_tp_profile(path_clouds, show = False, color = colors[0], ls = '-', label = f'Clouds, D={d} au')
    d = 1
    path_clouds = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Dplanet_{d}_clouds_correct2/'
    case.plot_tp_profile(path_clouds, show = False, color = colors[1], ls = '-', label = f'Clouds, D={d} au')
    d = 1.5
    path_clouds = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Dplanet_{d}_clouds_correct2/'
    case.plot_tp_profile(path_clouds, show = True, color = colors[2], ls = '-', label = f'Clouds, D={d} au')
    '''
    


    

    
    
   


   

    
    '''
    #compare the gas mixing rations between wrapper_hotJupiter_clouds_old_noZn and wrapper_hotJupiter_noclouds_old
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Tint_300_clouds_specresLR100/'
    #path = '/Users/louissiebenaler/ARCiS/Examples_Lumen/test_specresLR_200/'
    case.plot_tp_profile(path, show = False, color = 'orange', ls = '-', label = f'R = 100')
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Tint_300_clouds_specresLR60/'
    #path = '/Users/louissiebenaler/ARCiS/Examples_Lumen/test_specresLR_80/'
    #case.plot_tp_profile(path, show = False, color = 'magenta', ls = '-', label = f'R = 60')
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Tint_300_clouds_specresLR40/'
    #path = '/Users/louissiebenaler/ARCiS/Examples_Lumen/test_specresLR_40/'
    case.plot_tp_profile(path, show = False, color = 'black', ls = '-', label = f'R = 40')
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Tint_300_clouds_specresLR20/'
    #path = '/Users/louissiebenaler/ARCiS/Examples_Lumen/test_specresLR_20/'
    case.plot_tp_profile(path, show = False, color = 'red', ls = '-', label = f'R = 20')
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Tint_300_clouds_specresLR10/'
    #path = '/Users/louissiebenaler/ARCiS/Examples_Lumen/test_specresLR_10/'
    case.plot_tp_profile(path, show = True, color = 'green', ls = '-', label = f'R = 10')

    
   

  
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/Sonara_comparison/test_clouds_Sonora_lowerKzz_sigmadot-15/'
    case.plot_tp_profile(path, show = False, color = 'black', ls = '-', label = r'ARCiS: $\dot{\Sigma}_{\rm n} = 10^{-15}$, $K_{\rm zz} = 10^8$ cm$^2$/s')
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/test_clouds_Sonora_og/'
    case.plot_tp_profile(path, show = False, color = 'grey', ls = '--', label = r'ARCiS: $\dot{\Sigma}_{\rm n} = 10^{-15}$, $K_{\rm zz} = 10^8$ cm$^2$/s')
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/test_clouds_Sonora_AM_fsed1_DIFFUSE6/'
    case.plot_tp_profile(path, show = False, color = 'green', ls = '--', label = r'ARCiS: $\dot{\Sigma}_{\rm n} = 10^{-15}$, $K_{\rm zz} = 10^8$ cm$^2$/s')
                
    file = '/Users/louissiebenaler/Downloads/pressure-temperature_profiles 2/t900g31f1_m0.0_co1.0.pt'
    tmp = np.loadtxt(file, skiprows = 2)
    plt.plot(tmp[:, 2], tmp[:, 1], label = r'Sonora: clouds ($f_{\rm sed} = 1$)', color = 'blue', ls = '-')
    file = '/Users/louissiebenaler/Downloads/pressure-temperature_profiles 2/t900g31nc_m0.0_co1.0.pt'
    tmp = np.loadtxt(file, skiprows = 2)
    plt.plot(tmp[:, 2], tmp[:, 1], label = 'Sonora: gas', color = 'blue', ls = '--')
    plt.title('Files: t1800g31f1_m0.0_co1.0.pt, t1800g31nc_m0.0_co1.0.pt', fontsize = 10)
    plt.yscale('log')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (bar)')
    plt.legend(frameon = False)
    plt.gca().invert_yaxis()
    plt.show()


    case.get_Teq(output_path = path, bond_albedo = 0.1)


    
    #Al2O3[s], TiO2[s], CaTiO3[s], MgO[s], MgAl2O4[s], Fe[s], Cr[s], Ni[s], S[s], SiO[s], MgSiO3[s], Mg2SiO4[s], FeSiO3[s], Fe2SiO4[s]
    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/test_clouds_Sonora_AM_fsed1_compact_spheres_Teff1800K/'
    case.plot_cloud_mixing_ratios(output_path = path, specie = 'Fe[s]', show = False, label = '', color = 'black', ls = '-')
    case.plot_cloud_mixing_ratios(output_path = path, specie = 'Al2O3[s]', show = False, label = '', color = 'red', ls = '-')
    case.plot_cloud_mixing_ratios(output_path = path, specie = 'MgSiO3[s]', show = False, label = '', color = 'blue', ls = '-')

    path = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/test_clouds_Sonora_AM_fsed2_compact_spheres_Teff1800K/'
    case.plot_cloud_mixing_ratios(output_path = path, specie = 'Fe[s]', show = False, label = '', color = 'black', ls = '--')
    case.plot_cloud_mixing_ratios(output_path = path, specie = 'Al2O3[s]', show = False, label = '', color = 'red', ls = '--')
    case.plot_cloud_mixing_ratios(output_path = path, specie = 'MgSiO3[s]', show = True, label = '', color = 'blue', ls = '--')
    '''


    case.tp_sanity(output_path = '/Users/louissiebenaler/ARCiS/output_gas_EOSZ/case_88_Dplanet_0.01_Tint_1500_logg_1.5_met_-0.5/', plot = True)
    