import numpy as np 

def extract_pka_file_data(file, chain, time): 

    pkafile = open(file, 'r') 
    data_pka_list = [] 
    data_buriedness_list = [] 
    for line in pkafile: 
        line_processed = line.rstrip() 
        line_list = line_processed.strip().split() 
        if len(line_list) > 15 and line_list[2] == chain: 
            data_pka_list.append([line_list[0]+line_list[1], line_list[3]]) 
            data_buriedness_list.append([line_list[0]+line_list[1], line_list[4]]) 
    
    time_array = np.array([['Time [ps]'], [time]]) 
    data_pka_array = np.array(data_pka_list).T 
    data_pka = np.concatenate((time_array, data_pka_array), axis=1).tolist() 
    
    data_buriedness_array = np.array(data_buriedness_list).T 
    data_buriedness = np.concatenate((time_array, data_buriedness_array), axis=1).tolist() 
    
    return data_pka, data_buriedness