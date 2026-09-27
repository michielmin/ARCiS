import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, text
from pathlib import Path
import shutil
import subprocess
from multiprocessing import Pool
import os
from arcis_wrapper import ArcisWrapper
import pyARCiS


num_processes = 1


def run_simulation(process_id):
    file = f"arcis_noclouds_grid_new.dat"  
    grid = pd.read_csv(file, sep='\t')
    grid = grid.iloc[process_id::num_processes]  # Split the grid into chunks for each process
    input_file = '/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/arcis_gas.dat'
    output_dir = f'/Users/louissiebenaler/ARCiS/ARCIS_2026/src/Example/pyArcis_{process_id}/'
    wrapper = ArcisWrapper(input_file = input_file, output_path = output_dir)

    # Initialize once: opacity FITS files are read here.
    pyARCiS.pyinit(input_file, output_dir)
    pyARCiS.pyverbose(True)

    
    for sim in grid.to_dict(orient='records'):
        case_id, Dplanet, metallicity, Tint, logg, attempt, converged = sim['case_id'], sim['Dplanet'], sim['[M/H]'], sim['Tint'], sim['logg'], sim['attempt'], sim['converged']
        if case_id > 1528:
            if converged == 0:
                print(f"Running simulation for case_id: {case_id}, Dplanet: {Dplanet}, [M/H]: {metallicity}, Tint: {Tint}, logg: {logg}")

                #setting the parameters for the simulation
                pyARCiS.pysetvalue("Dplanet", Dplanet)
                pyARCiS.pysetvalue("TeffP", Tint)
                pyARCiS.pysetvalue("metallicity", metallicity)
                R = wrapper.define_radius(1, logg)
                pyARCiS.pysetvalue("Rp", R)
                #todo: could add something on how to initialize the TP profile based on Tint, to help find a solution faster and convergence
        
                # Recalculate the model using the in-memory opacity tables.
                #pyARCiS.pycomputemodel()
                #pyARCiS.pywritefiles()
                pyARCiS.pyrunarcis() #check if this also works
        
                final_dir = Path('./output_gas_failed_cases') / f"case_{case_id}_Dplanet_{Dplanet}_Tint_{Tint}_logg_{logg}_met_{metallicity}"
        
                # Avoid accidentally mixing old and new results.
                if final_dir.exists():
                    raise FileExistsError(
                        f"{final_dir} already exists; rename or remove it first."
                    )
                # Save all files generated for this run.
                shutil.copytree(output_dir, final_dir)

                #check if the simulation converged by looking for the presence of the output files
                convergence_file = final_dir / 'temperature_convergence.dat'
                with open(convergence_file, 'r') as f:
                    for line in f:
                        key, separator, value = line.strip().partition("=")

                        if key.strip().lower() == "converged":
                            converged =  value.strip().lower() == "true" #a boolean value indicating whether the simulation converged or not
                            if converged:
                                print(f"Simulation for case_id: {case_id} converged successfully.")
                                converged = 1
                            if not converged:
                                print(f"Simulation for case_id: {case_id} did not converge.")
                                converged = 0


                # Update the grid to indicate that this simulation has been attempted and converged.
                file = f"arcis_noclouds_grid_new.dat"  
                tmp = pd.read_csv(file, sep='\t')
                tmp.loc[tmp['case_id'] == case_id, 'attempt'] = 1
                tmp.loc[tmp['case_id'] == case_id, 'converged'] = converged
                tmp.to_csv(file, index=False, sep='\t')

            else:
                print(f"Skipping simulation for case_id: {case_id}, Dplanet: {Dplanet}, [M/H]: {metallicity}, Tint: {Tint}, logg: {logg} (already converged)")


            #break
        pass



if __name__ == "__main__":
    #num_processes = 2#16  # Number of parallel processes to run
    with Pool(processes=num_processes) as pool:
        # Run the simulations in parallel
        pool.map(run_simulation, range(num_processes))  # Adjust the range based on the number of simulations you want to run