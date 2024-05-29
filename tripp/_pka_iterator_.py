import MDAnalysis as mda 
from tripp._edit_pdb_ import edit_pdb 
from tripp._edit_pdb_ import mutate 
from propka import run 
from tripp._extract_pka_file_data_ import extract_pka_file_data 
import os 

def pka_iterator(trajectory_slices, core, universe, temp_name, mutation, chain, out, extract_surface_data): 
            
    start = trajectory_slices[core][0] 
    end = trajectory_slices[core][1]

    for index, ts in enumerate(universe.trajectory[start:end]):
        with mda.Writer(f'{temp_name}.pdb') as w:
            w.write(universe)
        edit_pdb(f'{temp_name}.pdb')
        if mutation != None: 
            mutate(temp_name, mutation) 
        run.single(f'{temp_name}.pdb')
        time = ts.time
        #Writing pKa csv
        header = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[0][0]) 
        data = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[0][1]).replace('*', '')
        if index == 0 and core == 0: 
            f = open(f'{out}_pka.csv', "w")
            f.write(header+'\n')
            f.write(data+'\n')
            f.close() 
        else: 
            f = open(f'{out}_pka.csv', "a")
            f.write(data+'\n')
            f.close()
        #Writing buridness csv if extract_surface_data set to true.
        if extract_surface_data == True: 
            header = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[1][0]) 
            data = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[1][1]).replace('*', '') 
            if index == 0 and core == 0: 
                f = open(f'{out}_surf.csv', "w")
                f.write(header+'\n')
                f.write(data+'\n')
                f.close() 
            else: 
                f = open(f'{out}_surf.csv', "a")
                f.write(data+'\n')
                f.close() 
    
    os.remove(f'{temp_name}.pdb') 
    os.remove(f'{temp_name}.pka')