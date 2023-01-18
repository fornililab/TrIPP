import MDAnalysis as mda 
import propka 
from propka import run
import numpy as np 
import os 
import multiprocessing as mp 

class Trajectory: 

    def __init__(self, trajectory_file, topology_file, cpu_number): 

        self.trajectory_file = trajectory_file 
        self.topology_file = topology_file 
        self.cpu_number = cpu_number 
        self.universe = mda.Universe(self.topology_file, self.trajectory_file) 
        frames_nr = len(self.universe.trajectory) 
        slices_nr = self.cpu_number 
        slice_length = round(frames_nr/slices_nr) 
        slices = [] 
        start_frame = 0 
        end_frame = slice_length 
        for i in range(slices_nr): 
            if end_frame > frames_nr: 
                slices.append([start_frame, frames_nr]) 
            else: 
                slices.append([start_frame, end_frame]) 
            start_frame+=slice_length 
            end_frame+=slice_length 
        self.trajectory_slices = slices 
        

    def calculate_pka(self, output_file, extract_surface_data=False, chain='A', mutation=None, thread=None): 

        def edit_pdb(pdbfile): 
            with open(pdbfile, 'r') as file: 
                data = file.read() 
                corrected_amino_acids = { 
                    'HSD' : 'HIS', 
                    'HSP' : 'HIS', 
                    'GLUP' : 'GLU ', 
                    'ASPP' : 'ASP ' 
                    }
                for correction in corrected_amino_acids: 
                    data = data.replace(correction, corrected_amino_acids[correction])
            with open(pdbfile, 'w') as file: 
                file.write(data) 
        
        def mutate(temp_name): 
            lines = [] 
            changed_lines = [] 
            deleted_lines = [] 
            with open(f'{temp_name}.pdb', 'r') as f: 
                lines = f.readlines() 
            for index, item in enumerate(lines): 
                line = item.split()
                if line[0] == 'ATOM' and line[5] == str(mutation): 
                    if line[2] == 'N' or line[2] == 'HN' or line[2] == 'CA' or line[2] == 'HA' or line[2] == 'CB' or line[2] == 'O' or line[2]  == 'C': 
                        changed_lines.append(index) 
                    else: 
                        deleted_lines.append(index) 
                    residue_type = line[3] 
            
            with open(f'{temp_name}.pdb', 'w') as fp:
                for number, line in enumerate(lines): 
                    if number not in deleted_lines:
                        if number in changed_lines: 
                            line = line.replace(residue_type, 'ALA') 
                        fp.write(line)
            return residue_type 
        
        def extract_data(file, chain, time): 

            pkafile = open(file, 'r') 
            data_pka_list = [] 
            data_surf_list = [] 
            for line in pkafile: 
                line_processed = line.rstrip() 
                line_list = line_processed.strip().split() 
                if len(line_list) > 15 and line_list[2] == chain: 
                    data_pka_list.append([line_list[0]+line_list[1], line_list[3]]) 
                    data_surf_list.append([line_list[0]+line_list[1], line_list[4]]) 
            
            time_array = np.array([['Time [ps]'], [time]]) 
            data_pka_array = np.array(data_pka_list).T 
            data_pka = np.concatenate((time_array, data_pka_array), axis=1).tolist() 
            
            data_surf_array = np.array(data_surf_list).T 
            data_surf = np.concatenate((time_array, data_surf_array), axis=1).tolist() 
            
            return data_pka, data_surf

        def pka_iterator(thread): 
            
            start = self.trajectory_slices[thread][0] 
            end = self.trajectory_slices[thread][1]
            
            if mutation != None: 
                out = f'{output_file}_{mutation}' 
                temp_name = f'temp_{mutation}_{thread}' 
            else: 
                out = output_file 
                temp_name = f'temp_{thread}' 

            for index, ts in enumerate(self.universe.trajectory[start:end]):
                with mda.Writer(f'{temp_name}.pdb') as w: 
                    w.write(self.universe) 
                edit_pdb(f'{temp_name}.pdb')
                if mutation != None: 
                    residue_type = mutate(temp_name) 
                run.single(f'{temp_name}.pdb')
                time = ts.time
                header = ','.join(extract_data(f'{temp_name}.pka', chain=chain, time=time)[0][0]) 
                data = ','.join(extract_data(f'{temp_name}.pka', chain=chain, time=time)[0][1]).replace('*', '')
                if index == 0 and thread == 0: 
                    f = open(f'{out}_pka.csv', "w")
                    f.write(header+'\n')
                    f.write(data+'\n')
                    f.close() 
                else: 
                    f = open(f'{out}_pka.csv', "a")
                    f.write(data+'\n')
                    f.close() 
                if extract_surface_data == True: 
                    header = ','.join(extract_data(f'{temp_name}.pka', chain=chain, time=time)[1][0]) 
                    data = ','.join(extract_data(f'{temp_name}.pka', chain=chain, time=time)[1][1]).replace('*', '') 
                    if index == 0 and thread == 0: 
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
        
        pka_iterator(thread) 

    

directory = '/Volumes/chris_drive/Simulations/GPR68/MD/Inactive/M0498/Mutations/I12A/MD1/'
t = Trajectory(f'{directory}GPR68I_I12A_M0498_POPC_eq_md_450ns_p_fit_skip100.xtc', f'{directory}GPR68I_I12A_M0498_POPC_min_sd.pdb', 8) 
def loop_function(index): 
    t.calculate_pka('./temp', extract_surface_data=True, chain='A', thread=index, mutation=64) 

if __name__ == '__main__':
    pool = mp.Pool(t.cpu_number)
    # Create jobs
    jobs = []
    for index, item in enumerate(t.trajectory_slices):
        # Create asynchronous jobs that will be submitted once a processor is ready
        job = pool.apply_async(loop_function, args=(index, ))
        jobs.append(job)
    # Submit jobs
    results = [job.get() for job in jobs]
    pool.close() 



