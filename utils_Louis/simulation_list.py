import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, text
from pathlib import Path
from shutil import copy
import subprocess
import os
from itertools import product


if __name__ == "__main__":
    Dplanet = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2]
    Z = [-0.5, 0, 0.3, 0.5, 0.7, 1, 1.5, 1.7, 2]
    Tint = [25,  50, 100, 150, 200, 300, 400, 500, 750, 1000, 1250, 1500]
    #g = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    g = [2.4, 2.8, 3.0, 3.2, 3.4, 3.6, 4, 4.4]
    Kzz = [1e6, 1e8, 1e10]
    Sigma_n = [1e-19, 1e-15, 1e-11]



    atmospheres = list(product(Dplanet, Z, Tint, g))
    cloud_cases = list(product(Kzz, Sigma_n))


    #produce the grid of gas atmosphere cases first
    rows = []
    for case_id, atmosphere in enumerate(atmospheres):
        dplanet, metallicity, tint, gravity = atmosphere
        rows.append({
                        "case_id": case_id,
                        "Dplanet": dplanet,
                        "[M/H]": metallicity,
                        "Tint": tint,
                        "logg": gravity,
                        "attempt": 0,
                        "converged": 0
                    })

    grid = pd.DataFrame(rows)
    grid.to_csv("arcis_noclouds_grid_new.dat", index=False, sep = '\t')



    #produce the grid of cloud atmosphere cases now
    rows = []

    for base_id, atmosphere in enumerate(atmospheres):
        dplanet, metallicity, tint, gravity = atmosphere

        for cloud_id, cloud in enumerate(cloud_cases):
            kzz, sigmadot = cloud
            case_id = base_id * len(cloud_cases) + cloud_id

            rows.append({
                "case_id": case_id,
                "Dplanet": dplanet,
                "[M/H]": metallicity,
                "Tint": tint,
                "logg": gravity,
                "Kzz": kzz,
                "Sigmadot": sigmadot,
                "attempt": 0,
                "converged": 0
            })

    grid = pd.DataFrame(rows)

    grid.to_csv("arcis_clouds_grid_new.dat", index=False, sep = '\t')


    #still want to add the following columns to the grid:
    #attempt column: 0 mean not attempted, 1 means attempted
    #convergence column: 0 means failed, 1 means successful

  