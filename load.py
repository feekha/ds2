# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 14:54:17 2022

@author: HENROG
"""
from typing import List
import os
import pandas as pd
import numpy as np
import pickle

def get_filenames_dict(type_list: List[str]=["sta","csv", "dat"]):
    wdir = os.getcwd()
    filenames_dict = {}
    for i in os.listdir(wdir + "\\data"):     
        filenames = os.listdir(wdir + "\\data" + "\\" + i)
        extensions = [x.rsplit(".")[-1] for x in filenames]
        for k in range(0,len(type_list)):
            if type_list[k] in extensions:
                files=[x for x in filenames if x.endswith("."+type_list[k])]
                filenames_dict[i] = files
                break
    return filenames_dict

def load_files(filenames_dict):
    inputdata_raw_dict = {}
    for folder in filenames_dict:
        rdat = list()
        for filename in filenames_dict[folder]:
            dir_str = "data\\" + folder + "\\" + filename
            extension = filename.rsplit(".")[-1]
            if folder == "mastdata_HAW":
                # put all possible mastdata formats here and later put it into funtions
                if extension == "csv":
                    dat = pd.read_csv(dir_str).loc[0:143]                  
                    rdat.append(dat)
                    device_id = "Mast HAW"
                    params = {}
            elif folder == "mastdata_Janneby":
                if extension == "dat":
                    dat = pd.read_csv(dir_str, skiprows=[0,2,3], header=0, dtype=str)
                    dat.iloc[:,1:] = np.float64(dat.iloc[:,1:])
                    rdat.append(dat)
                    device_id = "Mast Janneby"
                    params = {}
            else:
                # put all possible RSD data formats here and later into funtions                
                if extension == "csv":
                    df = pd.read_csv(dir_str)
                if extension == "sta":
                    array_sta = np.loadtxt(dir_str,
                                           dtype=str,
                                           skiprows=42,
                                           delimiter="\t"
                                           )                    
                    header_sta = np.loadtxt(dir_str,
                                            dtype=str,
                                            skiprows=41,
                                            max_rows=1,
                                            delimiter="\t"
                                            )
                    device_id = np.loadtxt(dir_str,
                                           dtype=str,
                                           skiprows=2,
                                           max_rows=1
                                           )[1][7:]
                    hei = np.loadtxt(dir_str,
                                           dtype=str,
                                           skiprows=39,
                                           max_rows=1
                                           )
                    wd_offset = np.loadtxt(dir_str,
                                           dtype=str,
                                           skiprows=29,
                                           max_rows=1,
                                           delimiter='='
                                           )[1]
                    coordinates = np.loadtxt(dir_str,
                                           dtype=str,
                                           skiprows=5,
                                           max_rows=1,
                                           delimiter='='
                                           )[1]
                    
                    coordinates = {'Lat': float(coordinates.split('Lat:')[1][:5]),
                                'Long': float(coordinates.split('Long:')[1][:5])}                                
                    hei = [float(h) for h in hei if h.replace('.','',1).isdigit()] # WC 2.0 and 2.1 different syntax --> this line of code will filter
                    params = {'hei': hei, 'wd_offset': float(wd_offset), 'coordinates': coordinates}
                    # put function here to check if wd_offset or hei settings change

                    if len(array_sta) == 0:
                        continue # in case there is no data in file
                    elif not isinstance(array_sta[0], np.ndarray):
                        array_sta = [array_sta] # in case ther is only one line of data in file
                    
                    array_sta_lst = []
                    for array in array_sta:
                         array_sta_lst.append(array.tolist()) # fixing errors with np.str_ data type
                    rdat.append(pd.DataFrame(array_sta_lst,columns=header_sta))   
        
        # put some functions here to take care of invalid data: different column names due to hei change
        # with open('rdat.pickle', 'wb') as handle:
        #     pickle.dump(rdat, handle, protocol=pickle.HIGHEST_PROTOCOL)     
        df = pd.concat(rdat)
        df.reset_index(inplace=True, drop=True)
        inputdata_raw_dict[device_id]= {'dat': df, 'param': params}
    return inputdata_raw_dict

                 
if __name__ == "__main__":
    # inputs
    type_list = ["sta", "csv", "dat"]  
    
    filenames_dict = get_filenames_dict(type_list)    
    inputdata_raw_dict = load_files(filenames_dict)
    
                    
                
            
    
