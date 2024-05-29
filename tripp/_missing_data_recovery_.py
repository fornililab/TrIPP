import MDAnalysis as mda 
from tripp._edit_pdb_ import edit_pdb 
from tripp._edit_pdb_ import mutate 
from propka import run 
from tripp._extract_pka_file_data_ import extract_pka_file_data 
import os 
import pandas as pd 

def missing_data_recovery(pka_file, surf_file, universe, mutation, chain): 

    if pka_file != None: 
        df_pka = pd.read_csv(pka_file) 
        timesteps = df_pka['Time [ps]'].to_list() 
        for index, ts in enumerate(universe.trajectory):
            time = ts.time 
            if time not in timesteps: 
                with mda.Writer(f'temp.pdb') as w:
                    w.write(universe)
                edit_pdb(f'temp.pdb')
                if mutation != None: 
                    mutate('temp', mutation) 
                run.single('temp.pdb')
                #Writing pKa csv
                data = ','.join(extract_pka_file_data('temp.pka', chain=chain, time=time)[0][1]).replace('*', '')
                f = open(pka_file, "a")
                f.write(data+'\n')
                f.close() 
                os.remove('temp.pdb') 
                os.remove('temp.pka')
    
    if surf_file != None: 
        df_surf = pd.read_csv(surf_file) 
        timesteps = df_surf['Time [ps]'].to_list() 
        for index, ts in enumerate(universe.trajectory):
            time = ts.time 
            if time not in timesteps: 
                with mda.Writer(f'temp.pdb') as w:
                    w.write(universe)
                edit_pdb(f'temp.pdb')
                if mutation != None: 
                    mutate('temp', mutation) 
                run.single('temp.pdb') 
                data = ','.join(extract_pka_file_data('temp.pka', chain=chain, time=time)[1][1]).replace('*', '') 
                f = open(surf_file, "a")
                f.write(data+'\n')
                f.close() 
                os.remove('temp.pdb') 
                os.remove('temp.pka') 