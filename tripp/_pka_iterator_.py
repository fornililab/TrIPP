import MDAnalysis as mda 
from tripp._edit_pdb_ import mutate 
from propka import run 
from tripp._extract_pka_file_data_ import extract_pka_file_data 
import os


def pka_iterator(trajectory_slices, core, universe, output_directory, mutation, chain, extract_buriedness_data, optargs=[]):

    temp_name = f'{output_directory}/.temp_{core}' 
        
    start = trajectory_slices[core][0] 
    end = trajectory_slices[core][1]
    
    if os.path.isfile(f'{output_directory}/.temp_pka_worker{core}.csv'):
        os.remove(f'{output_directory}/.temp_pka_worker{core}.csv')
    if extract_buriedness_data == True:
        if os.path.isfile(f'{output_directory}/.temp_buriedness_worker{core}.csv'):
            os.remove(f'{output_directory}/.temp_buriedness_worker{core}.csv')
            
    cwd = os.getcwd()
    for index, ts in enumerate(universe.trajectory[start:end]):
        with mda.Writer(f'{temp_name}.pdb') as w:
            w.write(universe)
        if mutation != None: 
            mutate(temp_name, mutation)
            
        os.chdir(output_directory)
        run.single(f'.temp_{core}.pdb',optargs=optargs)
        os.chdir(cwd)
        
        time = ts.time
        #Writing pKa csv
        header = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[0][0])
        data = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[0][1]).replace('*', '')
        with open(f'{output_directory}/.temp_pka_worker{core}.csv','a') as output:
            if index == 0 and core == 0:
                output.write(header+'\n')
            output.write(data+'\n')

        #Writing buriedness csv if extract_burideness_data set to true.
        if extract_buriedness_data == True:
            header = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[1][0]) 
            data = ','.join(extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)[1][1]).replace('*', '') 
            with open(f'{output_directory}/.temp_buriedness_worker{core}.csv','a') as output:
                if index == 0 and core == 0:
                    output.write(header+'\n')
                output.write(data+'\n')

    os.remove(f'{temp_name}.pdb') 
    os.remove(f'{temp_name}.pka')