import MDAnalysis as mda 
from tripp._edit_pdb_ import edit_pdb 
from tripp._edit_pdb_ import mutate 
from propka import run 
from tripp._extract_pka_file_data_ import extract_pka_file_data 
import os 

def pka_iterator(trajectory_slices, core, universe, temp_name, mutation, chain, extract_surface_data, optargs): 
            
    start = trajectory_slices[core][0] 
    end = trajectory_slices[core][1]
    if os.path.isfile(f'temp_pka_worker{core}.csv'):
        os.remove(f'temp_pka_worker{core}.csv')
        
    for index, ts in enumerate(universe.trajectory[start:end]):
        with mda.Writer(f'{temp_name}.pdb') as w:
            w.write(universe)
        if mutation != None: 
            mutate(temp_name, mutation) 
        run.single(f'{temp_name}.pdb',optargs=optargs)
        time = ts.time
        #Writing pKa csv
        header = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[0][0])
        data = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[0][1]).replace('*', '')
        with open(f'temp_pka_worker{core}.csv','a') as output:
            if index == 0 and core == 0:
                output.write(header+'\n')
            output.write(data+'\n')
            
        #Writing buridness csv if extract_surface_data set to true.
        if extract_surface_data == True: 
            header = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[1][0]) 
            data = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[1][1]).replace('*', '') 
            with open(f'temp_surf_worker{core}.csv','a') as output:
                if index == 0 and core == 0:
                    output.write(header+'\n')
                output.write(data+'\n')
    
    os.remove(f'{temp_name}.pdb') 
    os.remove(f'{temp_name}.pka')