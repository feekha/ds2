
import paramiko
import pysftp
from credentials import HOST, USERNAME, PASSWORD
import pandas as pd
import numpy as np
from typing import List
import os
from tqdm import tqdm
import time
import pickle

def get_filenames_dict(sftp, type_list: List[str]=["sta","csv","dat"]):
    # wdir = os.getcwd()
    print('get filenames...')
    wdir = sftp.pwd
    print(wdir)
    filenames_dict = {}
    print(sftp.listdir(wdir))

    i = 1
    for loc in sftp.listdir(wdir):
        # mast data
        filenames_mast = sftp.listdir(wdir + "/" + loc + "/Mast")
        print(i , 'loc' , loc ,  'filenames_mast' , filenames_mast)
        extensions = [x.rsplit(".")[-1] for x in filenames_mast]
        print(i , 'loc' , loc , 'extensions' , extensions)
        for k in range(0,len(type_list)):
            if type_list[k] in extensions:
                files=[x for x in filenames_mast if x.endswith("."+type_list[k])]
                filenames_dict[wdir + "/" + loc + "/Mast"] = files
                break

        print(i , filenames_dict)

        # rsd data
        ip_folders = sftp.listdir(wdir + "/" + loc + "/RSD")

        print(i , 'ip_folder' , ip_folders)
        for folder in ip_folders:
            filenames = sftp.listdir(wdir + "/" + loc + "/RSD/" + folder)

            print(i , 'folder' , folder , 'filenames' , filenames)

            folder_dir = wdir + "/" + loc + "/RSD/" + folder

            print(i , 'folder' , folder , 'folder_dir' , folder_dir)

            if "STA" in filenames:
                filenames = sftp.listdir(wdir + "/" + loc + "/RSD/" + folder + "/STA")

                print(i , 'folder' , folder , 'filename' , filenames)

                folder_dir = wdir + "/" + loc + "/RSD/" + folder + "/STA"

                print(i , 'folder' , folder , 'folder_dir' , folder_dir)

            extensions = [x.rsplit(".")[-1] for x in filenames] 

            print(i , 'folder' , folder , 'extension' , extensions)

            for k in range(0,len(type_list)):

                print('k' , k , 'type_list' , type_list)

                if type_list[k] in extensions:

                    print('yes')

                    for x in filenames:
                        print('x' , x)
                        if x.endswith('.' + type_list[k]):
                            print(1 , 'yes')
                            if sftp.isfile(folder_dir + '/' + x):
                                print(2 , 'yes')
                                print(x)
                    
                    files=[x for x in filenames if x.endswith("."+type_list[k]) & sftp.isfile(folder_dir + '/' + x)] # this for loop will need some time, because sftp.isfile needs 0.03s
                    
                    print(i , 'folder' , folder , 'files' , files)

                    filenames_dict[folder_dir] = files

                    print(i , 'folder' , folder , 'filenames_dict' , filenames_dict)

                    break

        i = i + 1

    print(filenames_dict)
    return filenames_dict


def load_files(filenames_dict, sftp):
    inputdata_raw_dict = {}
    print('load files...')
    for folder in tqdm(filenames_dict):
        rdat = list()
        params = {}
        for filename in filenames_dict[folder]:
            dir_str = folder + "/" + filename
            # print(dir_str)
            extension = filename.rsplit(".")[-1]           
            if 'HAW/Mast' in folder:
                # put all possible mastdata formats here and later put it into funtions
                if extension == "csv":
                    with sftp.open(dir_str,'r',bufsize=-1) as file:
                        dat = pd.read_csv(file).loc[0:143]                  
                    rdat.append(dat)
                    device_id = "Mast HAW"
            elif 'Janneby/Mast' in folder:
                if extension == "dat":
                    with sftp.open(dir_str,'r',bufsize=-1) as file:
                        dat = pd.read_csv(file, skiprows=[0,2,3], header=0, dtype=str)
                    dat.iloc[:,1:] = np.float64(dat.iloc[:,1:])
                    rdat.append(dat)
                    device_id = "Mast Janneby"
            else:
                # put all possible RSD data formats here and later into funtions                
                if extension == "csv": # ZX 
                    continue
                    # df = pd.read_csv(dir_str)
                if extension == "sta":
                    with sftp.open(dir_str,'r',bufsize=-1) as file:
                        file_bytes = file.read()
                    file_str = file_bytes.decode(errors='replace')
                    file_str = file_str.replace('\ufffd', '°')
                    file_lines = file_str.splitlines()
                    
                    if len(file_lines) <= 42: # filters empty files
                        print('datafile ' + dir_str + ' is empty')
                        continue
                    
                    if 'Header' not in file_lines[0]:
                        print('datafile ' + dir_str + ' has no header')
                        continue
                    
                    file_cells = [row.split('\t') for row in file_lines[41:]]

                    for i in range(len(file_cells)-1,-1,-1): # filters rows with invalid number of columns --> backwards because of del 
                        if len(file_cells[0]) != len(file_cells[i]):
                            print('datafile ' + dir_str + ' has invalid number of columns')
                            del file_cells[i]
                            
                    # with open('file_cells.pickle', 'wb') as handle:
                    #     pickle.dump(file_cells, handle, protocol=pickle.HIGHEST_PROTOCOL)  
                    file_df = pd.DataFrame(file_cells[1:],columns = file_cells[0])
                    device_id = file_lines[2][10:]
                    coordinates = file_lines[5]
                    
                    if 'Lat' not in coordinates: # in case of hidden coordinates
                        coordinates = {'Lat': 54.0, 'Long': 9.0} # maybe build something that gives correct coordinates 
                    else:
                        coordinates = {'Lat': float(coordinates.split('Lat:')[1][:5]),
                                       'Long': float(coordinates.split('Long:')[1][:5])}  
                    hei = file_lines[39].split('=')[1].split('\t')
                    hei = [float(h) for h in hei if h.isdigit()]
                   
                    if 'hei' in params:
                        if hei != params['hei']:
                            print('detected change in heights')
                            rdat = list() # only picks last heights, filters change in heights, faults when height changes in the end                  
                    
                    wd_offset = float(file_lines[29].split('=')[1])
                    params = {'hei': hei, 'wd_offset': float(wd_offset), 'coordinates': coordinates}
                    
                    rdat.append(file_df)   
        
        if len(rdat) == 0: # in case of empty directories or only invalid data
            continue       
        # put some functions here to take care of invalid data
        # with open('rdat.pickle', 'wb') as handle:
        #     pickle.dump(rdat, handle, protocol=pickle.HIGHEST_PROTOCOL)     
        df = pd.concat(rdat)
        df.reset_index(inplace=True, drop=True)
        inputdata_raw_dict[device_id]= {'dat': df, 'param': params}
    return inputdata_raw_dict


def get_sftp():
    paramiko.util.log_to_file("paramiko.log") # exclude this after debugging    
    cnopts = pysftp.CnOpts()
    cnopts.hostkeys = None           
    with pysftp.Connection(host=HOST, username=USERNAME, password=PASSWORD, cnopts=cnopts) as sftp:
        print('successfully connected to sftp')
        filenames_dict = get_filenames_dict(sftp)
        # test_files = {}
        # test_files['/lidardashboard/Janneby/RSD/Local_IP_16/STA'] = filenames_dict['/lidardashboard/Janneby/RSD/Local_IP_16/STA']
        inputdata_raw_dict = load_files(filenames_dict, sftp)
        print('data loaded successfully, connection will be closed now')       
    return inputdata_raw_dict
        
    
if __name__ == "__main__":
    
    paramiko.util.log_to_file("paramiko.log") # exclude this after debugging
    
    cnopts = pysftp.CnOpts()
    cnopts.hostkeys = None    
        
    with pysftp.Connection(host=HOST, username=USERNAME, password=PASSWORD, cnopts=cnopts) as sftp:
        print('successfully connected to sftp')
        
        filenames_dict = get_filenames_dict()
        
        # test_files = {}
        # test_files['/lidardashboard/Janneby/RSD/Local_IP_16/STA'] = filenames_dict['/lidardashboard/Janneby/RSD/Local_IP_16/STA']
        
        inputdata_raw_dict = load_files(filenames_dict)
        
        print('data loaded successfully, connection will be closed now')


    
    
    
    