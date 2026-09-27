from __future__ import annotations
import numpy as np 
import pandas as pd
import glob
import h5py
from scipy.constants import N_A
from scipy.interpolate import RegularGridInterpolator
import scipy.interpolate 
from scipy.constants import Boltzmann, erg, N_A
from matplotlib import pyplot as plt, text
from pathlib import Path
from typing import Union
k_B_cgs = Boltzmann / erg  # Convert from J/K to erg/K 


def _fortran_scientific(value: float, width: int = 29, precision: int = 15) -> str:
    """Format a number with the three-digit exponent used by the CIA files."""
    mantissa, exponent = f"{value:.{precision}E}".split("E")
    result = f"{mantissa}E{int(exponent):+04d}"
    if len(result) > width:
        raise ValueError(f"Value {value!r} does not fit in a {width}-character field")
    return result.rjust(width)


def write_cia_file(
    data: np.ndarray,
    temperatures: np.ndarray,
    output_file: Union[str, Path],
    collision_pair: str = "H2-H2",
) -> None:
    """Write an ARCiS-compatible CIA file.

    Parameters
    ----------
    data
        Array with shape ``(n_temperature, 2, n_wavenumber)``.
    temperatures
        One temperature in kelvin for each element of the first data axis.
    output_file
        Destination CIA filename.
    collision_pair
        Species label placed in the header, for example ``"H2-H2"``.
    """
    data = np.asarray(data, dtype=np.float64)
    temperatures = np.asarray(temperatures, dtype=np.float64).reshape(-1)
    output_file = Path(output_file)

    if data.ndim != 3 or data.shape[1] != 2:
        raise ValueError(
            "data must have shape (n_temperature, 2, n_wavenumber); "
            f"received {data.shape}"
        )
    if data.shape[0] != temperatures.size:
        raise ValueError(
            f"data contains {data.shape[0]} temperature blocks, but "
            f"{temperatures.size} temperatures were supplied"
        )
    if data.shape[2] < 2:
        raise ValueError("Each temperature block must contain at least two points")
    if not collision_pair or len(collision_pair) > 20:
        raise ValueError("collision_pair must contain between 1 and 20 characters")
    if not np.all(np.isfinite(data)):
        raise ValueError("data contains NaN or infinite values")
    if not np.all(np.isfinite(temperatures)):
        raise ValueError("temperatures contains NaN or infinite values")
    if np.any(temperatures <= 0.0):
        raise ValueError("All temperatures must be positive")
    if np.any(data[:, 1, :] < 0.0):
        raise ValueError("CIA cross sections must be non-negative")

    # ARCiS stops if temperature blocks are not strictly increasing.
    temperature_order = np.argsort(temperatures, kind="stable")
    temperatures = temperatures[temperature_order]
    data = data[temperature_order]
    if np.any(np.diff(temperatures) <= 0.0):
        raise ValueError("Temperatures must be unique")

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="ascii", newline="\n") as handle:
        for block, temperature in zip(data, temperatures):
            wavenumber = block[0]
            cross_section = block[1]

            point_order = np.argsort(wavenumber, kind="stable")
            wavenumber = wavenumber[point_order]
            cross_section = cross_section[point_order]

            if np.any(wavenumber <= 0.0):
                raise ValueError("All wavenumbers must be positive")
            if np.any(np.diff(wavenumber) <= 0.0):
                raise ValueError(
                    f"Wavenumbers must be unique at T={temperature:g} K"
                )

            n_points = wavenumber.size
            minimum = f"{wavenumber[0]:.3f}"
            maximum = f"{wavenumber[-1]:.3f}"
            count = str(n_points)
            temp = f"{temperature:.1f}"

            if len(minimum) > 10 or len(maximum) > 10:
                raise ValueError("A wavenumber does not fit in the CIA header")
            if len(count) > 7:
                raise ValueError("The number of points does not fit in the CIA header")
            if len(temp) > 7:
                raise ValueError("A temperature does not fit in the CIA header")

            # Matches the Fortran format (a20, a10, a10, a7, a7) in CIA.f.
            handle.write(
                f"{collision_pair:>20}"
                f"{minimum:>10}"
                f"{maximum:>10}"
                f"{count:>7}"
                f"{temp:>7}\n"
            )

            handle.writelines(
                f"{nu:19.13f}{_fortran_scientific(sigma)}\n"
                for nu, sigma in zip(wavenumber, cross_section)
            )






if __name__ == "__main__":
    folder = '/Volumes/L7aler_HD/PhD/Snellius/cross_sections/CIA_H2H2/combined/'
    files = glob.glob(folder + '/*.hdf5')
    files.sort(key=lambda x: float(x.split('_T')[-1].split('.hdf5')[0]))  # Sort by temperature

    data = np.zeros((len(files), 2, 19981))
    T = []
    for i, file in enumerate(files):
        f = h5py.File(file, 'r+')
        t = f['T'][:]
        print('Temperature:', t)
        p = 10**f['P'][0] * 1e1 #pressure in cgs
        k = 1 / (f['wave'][:] * 1e2)
        xsec = 10**f['cross_sec'][:, 0, 0] * 1e4 #units are in cm^2 now
        tmp = xsec * k_B_cgs * t / p #units that the HITRAN data is in 

        #interpolate the data to a smaller wavenumber grid (to save some memory)
        interp_func = scipy.interpolate.interp1d(k, tmp, kind='linear', bounds_error=False, fill_value=0)
        k_new = np.arange(min(k), 20001, 1)
        xsec_new = interp_func(k_new)

        data[i, 0, :] = k_new
        data[i, 1, :] = xsec_new
        T.append(t)
        




        '''
        plt.plot(1 / k_new * 1e4, xsec_new)
        file = '/Users/louissiebenaler/ARCiS/Data/CIA/tmp.dat'
        d = np.loadtxt(file, skiprows=1)
        plt.plot(1 / d[:, 0] * 1e4, d[:, 1], ls = '--', color = 'black')
        plt.xscale('log')
        plt.yscale('log')
        plt.show()
        '''

        break



    #write_cia_file(data, np.array(T), '/Users/louissiebenaler/ARCiS/Data/CIA/H2-H2_Louis.cia', collision_pair='H2-H2')

