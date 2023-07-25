# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 11:27:03 2022

@author: HENROG
"""
import load
import pandas as pd
import numpy as np
import collections

# later make device and entity as base class
class Entity:
    
    def remove_duplicate_timestamps(self, inputdata):
        inputdata.drop_duplicates(inputdata.columns[0], inplace=True, keep=False, ignore_index=True)
        return inputdata
    
    def get_timestamps(self, inputdata):
        jdi = inputdata.iloc[:,0]
        return jdi
    
    def adjust_period_general(self, start_idx, end_idx):
        # self.inputdata = self.inputdata.loc[start_idx[0]:end_idx[0]]
        self.jdR = self.jdR.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.ws = self.ws.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.ws_std = self.ws_std.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.ti = self.ti.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.wd = self.wd.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
                
    
    def calculate_TI(self, ws_std, ws):
        ti = ws_std/ws.to_numpy()
        for i in ti.columns:
            ti.rename(columns={str(i): str(i.rsplit(" ")[0]) + " TI"}, inplace=True)
        return ti 

class RSD(Entity):
    
    def apply_time_shift(self, **kwargs):
        '''
        Adds or subtracts a time delta to/from windcube time series.

        Parameters
        ----------
        **kwargs : int, float
            Available kwargs: 
                * ‘W’, ‘D’, ‘T’, ‘S’, ‘L’, ‘U’, or ‘N’
                * ‘days’ or ‘day’
                * ‘hours’, ‘hour’, ‘hr’, or ‘h’
                * ‘minutes’, ‘minute’, ‘min’, or ‘m’
                * ‘seconds’, ‘second’, or ‘sec’
                * ‘milliseconds’, ‘millisecond’, ‘millis’, or ‘milli’
                * ‘microseconds’, ‘microsecond’, ‘micros’, or ‘micro’
                * ‘nanoseconds’, ‘nanosecond’, ‘nanos’, ‘nano’, or ‘ns’
        '''          
        self.jdR = self.jdR.apply(lambda x: x+pd.Timedelta(**kwargs))
            
    def apply_wd_shift(self, wd_offset):
        self.wd = self.wd.applymap(lambda x: (x-wd_offset)%360, na_action='ignore')


class Mast(Entity):
    
    def format_data(self, inputdata_raw):
        inputdata_raw = inputdata_raw.replace(["",999.00,9990,9991,9992,9993,9994,9995,9996,9997,9998,9999,"#N/A","NA","#NA"],float("NaN"))
        inputdata_raw.iloc[:,0] = pd.to_datetime(inputdata_raw.iloc[:,0])
        return inputdata_raw
    
    def dataset_assignment(self, inputdata, idata):
        data_dict = {}
        for item in idata:    
            data = []
            for i in idata[item]:
                data.append(inputdata.iloc[:,i])
            data_dict[item] = pd.DataFrame(data).transpose()
        return data_dict
    
    def opposite_cup_avg(self, wsA, wsB, wd, wdLowA, wdUpA, wdLowB, wdUpB):
        j_av = np.logical_and(np.logical_or(wd.to_numpy() < wdLowA, wd.to_numpy() > wdUpA),np.logical_or(wd.to_numpy() < wdLowB, wd.to_numpy() > wdUpB))
        j_A = np.logical_and(wd.to_numpy() >= wdLowA, wd.to_numpy() <= wdUpA)
        j_B = np.logical_and(wd.to_numpy() >= wdLowB, wd.to_numpy() <= wdUpB)
        
        wsRes = pd.DataFrame(np.full((len(wsA),2), np.nan))
        wsDiff = pd.DataFrame(np.full((len(wsA),1), 0))
        
        wsRes.iloc[j_av,0] = (wsA.iloc[j_av,0]+wsB.iloc[j_av,0])/2
        wsRes.iloc[j_A,0] = wsA.iloc[j_A,0]
        wsRes.iloc[j_B,0] = wsB.iloc[j_B,0]
        wsRes.iloc[j_av,1] = (wsA.iloc[j_av,1]+wsB.iloc[j_av,1])/2
        wsRes.iloc[j_A,1] = wsA.iloc[j_A,1]
        wsRes.iloc[j_B,1] = wsB.iloc[j_B,1]
        wsDiff[j_av] = wsA.iloc[j_av,0]-wsB.iloc[j_av,0]
        return wsRes, wsDiff
    
    def combine_oca(self, oca):
        WS_Mast = pd.concat([oca[j][0].iloc[:,0] for j in range(0,len(oca))],axis=1)
        WS_std_Mast = pd.concat([oca[j][0].iloc[:,1] for j in range(0,len(oca))],axis=1)
        WSD_Mast = pd.concat([oca[j][1] for j in range(0,len(oca))],axis=1)
        return WS_Mast, WS_std_Mast, WSD_Mast
    
    def adjust_period(self, start_idx, end_idx):
        super().adjust_period_general(start_idx, end_idx)
        self.wsDiff = self.wsDiff.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.wd_std = self.wd_std.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.humi = self.humi.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.press = self.press.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.temp = self.temp.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
    
    def initial_preprocessing(self, idata):
        self.inputdata = super().remove_duplicate_timestamps(self.inputdata)      
        self.data = self.dataset_assignment(self.inputdata,idata)               
        self.jdR = super().get_timestamps(self.inputdata)
               
    def create_variables(self):
        # create rest of df
        self.ti = super().calculate_TI(self.ws_std, self.ws)
        self.wd = self.data["wd"][self.data["wd"].columns[::2]]
        self.wd.columns = [str(h)+"m WD" for h in self.vane_hei]
        self.wd_std = self.data["wd"][self.data["wd"].columns[1::2]]
        self.wd_std.columns = [str(h)+"m WD_std" for h in self.vane_hei]
        self.humi = self.data["humi"]
        self.humi.columns = [str(h)+"m Humidity" for h in self.humi_hei]
        self.press = self.data["press"]
        self.press.columns = [str(h)+"m Pressure" for h in self.press_hei]
        self.temp = self.data["temp"]
        self.temp.columns = [str(h)+"m Temperature" for h in self.temp_hei]


class Mast_HAW(Mast):
    
    def __init__(self, inputdata_raw):
        self.sn = "Mast HAW"
        self.coordinates = {'Lat': 53.47098, 'Long': 10.20683}
        
        self.inputdata = super().format_data(inputdata_raw['dat'])
        # inputs: ----------
        self.hei = [120.5, 80, 65, 45]
        self.vane_hei = [112, 62.5]
        self.temp_hei = [113.1, 10]
        self.humi_hei = [113.1, 10]
        self.press_hei = [113.1, 10]
        self.cup_orient = [313, 133, 313, 133, 313, 133, 313, 133]
        vane_idx = [0,0,2,2]
        idata = {}
        iws = np.array([1,6,21,26,36,41,46,51])
        idata["ws"] = iws
        idata["wsMsd"] = iws + 3
        idata["wd"] = np.array([59,60,68,69])
        idata["temp"] = np.array([86,96])
        idata["humi"] = np.array([71,76])
        idata["press"] = np.array([101,106])
        # ------------------
        
        super().initial_preprocessing(idata)
        
        # opposite cup average
        oca = []
        for i in range(0,len(self.cup_orient),2):
            k = int(i/2)
            oca.append(self.opposite_cup_avg(
                pd.concat([self.data["ws"].iloc[:,i],self.data["wsMsd"].iloc[:,i]],axis=1),
                pd.concat([self.data["ws"].iloc[:,i+1],self.data["wsMsd"].iloc[:,i+1]],axis=1),
                self.data["wd"].iloc[:,vane_idx[k]],
                self.cup_orient[i]-30,
                self.cup_orient[i]+30,
                self.cup_orient[i+1]-30,
                self.cup_orient[i+1]+30))
            oca[k][0].columns = [str(self.hei[k]) + "m WS", str(self.hei[k]) + "m WS_std"]
            oca[k][1].columns = [str(self.hei[k]) + "m WSD"]
        self.ws, self.ws_std, self.wsDiff = self.combine_oca(oca)
        
        super().create_variables()              
            
    def adjust_period(self, start_idx, end_idx):
        super().adjust_period(start_idx, end_idx)


class Mast_Janneby(Mast):
    
    def __init__(self, inputdata_raw):
        self.sn = "Mast Janneby"
        self.coordinates = {'Lat': 54.637616, 'Long': 9.316778}
        
        self.inputdata = super().format_data(inputdata_raw['dat'])
        # inputs: ----------
        self.hei = [103, 75, 57, 29] # 29m excluded ??
        self.vane_hei = [96, 54]
        self.temp_hei = [95, 10]
        self.humi_hei = [95, 10]
        self.press_hei = [95, 10]
        # self.cup_orient = [330, 150, 330, 150, 330, 150, 330] # 2x 103m, 2x 75m, 2x 57m, 1x 29m --- before refit
        self.cup_orient = [150, 330, 150, 330, 150, 330] # 1x 103m, 2x 75m, 2x 57m, 1x 29m --- after refit
        # cup_ignore_range_oca = [(30,30),(20,25),(15,25),(20,20),(15,25),(15,20)] # --- before refit
        cup_ignore_range_oca = [(15,25),(20,20),(15,25),(15,20)] # --- after refit
        vane_idx = [0,0,2,2]
        idata = {}
        # idata["ws"] = np.array([16,17,18,19,20,21,22]) # --- before refit
        idata["ws"] = np.array([17,18,19,20,21,22]) # --- after refit
        # idata["wsMsd"] = np.array([34,35,36,37,38,39,40]) # --- before refit
        idata["wsMsd"] = np.array([35,36,37,38,39,40]) # --- after refit
        idata["wd"] = np.array([8,9,10,11])
        idata["temp"] = np.array([29,30])
        idata["humi"] = np.array([31,32])
        idata["press"] = np.array([27,28])
        
        # before refit
        # slopes = [0.98835, 0.98846, 0.98941, 0.98854, 0.98847, 0.98950, 0.98814]
        # offsets = [0.0625, 0.0622, 0.0615, 0.0604, 0.0585, 0.0487, 0.067]
        
        # after refit
        slopes = [0.98783, 0.99002, 0.98742, 0.98938, 0.98876, 0.98852]
        offsets = [0.0525, 0.0270, 0.0645, 0.0189, 0.0463, 0.0374]
        # ------------------
        
        super().initial_preprocessing(idata)
        self.ws_correction(slopes, offsets)
        self.wd_correction()

        # opposite cup average
        oca = []
        # for i in range(0,len(cup_ignore_range_oca),2): # --- before refit
        for i in range(1,len(cup_ignore_range_oca),2): # --- after refit
            k = int(i/2+1)
            oca.append(self.opposite_cup_avg(
                pd.concat([self.data["ws"].iloc[:,i],self.data["wsMsd"].iloc[:,i]],axis=1),
                pd.concat([self.data["ws"].iloc[:,i+1],self.data["wsMsd"].iloc[:,i+1]],axis=1),
                self.data["wd"].iloc[:,vane_idx[k]],
                self.cup_orient[i]-cup_ignore_range_oca[i-1][0],
                self.cup_orient[i]+cup_ignore_range_oca[i-1][1],
                self.cup_orient[i+1]-cup_ignore_range_oca[i][0],
                self.cup_orient[i+1]+cup_ignore_range_oca[i][1]))
            # oca[k][0].columns = [str(self.hei[k]) + "m WS", str(self.hei[k]) + "m WS_std"] # --- before refit
            # oca[k][1].columns = [str(self.hei[k]) + "m WSD"] # --- before refit
            oca[k-1][0].columns = [str(self.hei[k]) + "m WS", str(self.hei[k]) + "m WS_std"] # --- after refit
            oca[k-1][1].columns = [str(self.hei[k]) + "m WSD"] # --- after refit
        self.ws, self.ws_std, self.wsDiff = self.combine_oca(oca)
        
        # --- before refit
        # self.ws['29m WS'] = self.data['ws'].iloc[:,6]
        # self.ws_std['29m WS_std'] = self.data['wsMsd'].iloc[:,6]
        # self.wsDiff['29m WSD'] = [0 for x in range(len(self.wsDiff))]
        
        # --- after refit
        self.ws['29m WS'] = self.data['ws'].iloc[:,5]
        self.ws_std['29m WS_std'] = self.data['wsMsd'].iloc[:,5]
        self.wsDiff['29m WSD'] = [0 for x in range(len(self.wsDiff))]
        ws_cols = self.ws.columns
        ws_std_cols = self.ws_std.columns
        wsDiff_cols = self.wsDiff.columns
        self.ws['103m WS'] = self.data['ws'].iloc[:,0]
        self.ws_std['103m WS_std'] = self.data['wsMsd'].iloc[:,0]
        self.wsDiff['103m WSD'] = [0 for x in range(len(self.wsDiff))]
        self.ws = self.ws[['103m WS']+list(ws_cols)]
        self.ws_std = self.ws_std[['103m WS_std']+list(ws_std_cols)]
        self.wsDiff = self.wsDiff[['103m WSD']+list(wsDiff_cols)]
        
        # wire disturbances        
        self.disturbance_correction(ws='57m WS',wsd='57m WSD',wsB='WS057SE_Avg',wd=self.data["wd"].iloc[:,2],wdLow=200,wdUp=240)
        self.disturbance_correction(ws='57m WS',wsd='57m WSD',wsB='WS057NW_Avg',wd=self.data["wd"].iloc[:,2],wdLow=260,wdUp=300)
        self.disturbance_correction(ws='75m WS',wsd='75m WSD',wsB='WS075NW_Avg',wd=self.data["wd"].iloc[:,0],wdLow=180,wdUp=220)
        self.disturbance_correction(ws='75m WS',wsd='75m WSD',wsB='WS075SE_Avg',wd=self.data["wd"].iloc[:,0],wdLow=260,wdUp=305)
        
        super().create_variables()    
        
    def ws_correction(self, slopes, offsets):
        for i in range(self.data["ws"].shape[1]):
            self.data["ws"].iloc[:,i] = self.data["ws"].iloc[:,i].mul(slopes[i]/100)
            self.data["ws"].iloc[:,i] = self.data["ws"].iloc[:,i].add(offsets[i])
        for i in range(self.data["wsMsd"].shape[1]):
            self.data["wsMsd"].iloc[:,i] = self.data["wsMsd"].iloc[:,i].div(100)
            
    def wd_correction(self):
        self.data["wd"].iloc[:,0] = self.data["wd"].iloc[:,0].add(-15.5)%360
        self.data["wd"].iloc[:,2] = self.data["wd"].iloc[:,2].add(4)%360
        
    def disturbance_correction(self, ws, wsd, wsB, wd, wdLow, wdUp):
        j = np.logical_and(wd.to_numpy() >= wdLow, wd.to_numpy() <= wdUp)
        self.ws[ws][j] = self.data["ws"][wsB][j]
        self.wsDiff[wsd][j] = self.wsDiff[wsd][j].mul(0)

    
    def adjust_period(self, start_idx, end_idx):
        super().adjust_period(start_idx, end_idx)


class Windcube(RSD):
    
    def __init__(self, inputdata_raw, sn):     
        self.sn = sn
        self.type_rsd = self.type_device(sn)
        self.coordinates = inputdata_raw['param']['coordinates']
        
        self.hei = inputdata_raw['param']['hei']
        self.wd_offset = inputdata_raw['param']['wd_offset'] # just saved here, not needed, because corrected in device itself, but to check if applied correctly
        self.inputdata = self.format_data(inputdata_raw['dat'])
        if self.type_rsd == "V2.0":
            hei_offset = 1
            self.hei = [h+hei_offset for h in self.hei]
            self.inputdata = self.apply_height_offset(self.inputdata,hei_offset)
        self.inputdata = super().remove_duplicate_timestamps(self.inputdata)
        self.jdR = super().get_timestamps(self.inputdata)
        # self.hei = self.get_heights(self.inputdata)
        self.ws, self.ws_std, self.ti, self.wd, self.cnr, self.ava = self.dataset_assignment(self.inputdata)
        self.ws, self.ws_std, self.ti, self.wd, self.cnr, self.ava = self.internal_QC(self.ws, self.ws_std, self.ti, self.wd, self.cnr, self.ava)
                
    def format_data(self, inputdata_raw):
        inputdata_raw = inputdata_raw.replace(["","999.00","9990","9991","9992","9993","9994","9995","9996","9997","9998","9999","#N/A","NA","#NA"],"NaN")
        inputdata_raw.iloc[:,1:] = inputdata_raw.iloc[:,1:].astype(np.float64, errors="raise")
        inputdata_raw.iloc[:,0] = pd.to_datetime(inputdata_raw.iloc[:,0])
        return inputdata_raw
    
    def type_device(self, sn):
        if "WLS866-" in sn:
            type_rsd = "Offshore"
        elif len(sn.rsplit("-")[-1]) == 3:
            type_rsd = "V2.0"
        else:
            type_rsd = "V2.1"
        return type_rsd
    
    # will just affect column names
    def apply_height_offset(self, inputdata, hei_offset):
        for i in inputdata.columns:
            y = i.split(" ",1)
            if len(y[0])<3 or len(y[0])>4:
                continue
            if y[0][-1] != "m":
                continue
            new_hei = str(int(y[0][:-1])+hei_offset) + "m " + y[1]
            inputdata.rename(columns={str(i): new_hei}, inplace=True)
        return inputdata
    
    # possibility to get height from dataframe columns 
    def get_heights(self, inputdata):
        hei = []
        for i in inputdata.columns:
            y = i.split(" ",1)[0]
            if len(y)<3 or len(y)>4:
                continue
            if (y[-1] != "m") or y in hei:
                continue
            hei.append(y)
        return hei
    
    def dataset_assignment(self, inputdata):
        ws = inputdata.iloc[:,7::12]
        ws_std = inputdata.iloc[:,8::12]
        ti = super().calculate_TI(ws_std, ws)
        wd = inputdata.iloc[:,11::12]
        cnr = inputdata.iloc[:,14::12]
        ava = inputdata.iloc[:,17::12]
        return ws, ws_std, ti, wd, cnr, ava
    
    def internal_QC(self, ws, ws_std, ti, wd, cnr, ava):
        iikeep = (ava.to_numpy() >= 80.0) & (cnr.to_numpy() >= -23.0) & (cnr.to_numpy() <= 18.0) 
        ws = ws[pd.DataFrame(iikeep,columns=ws.columns)]
        ws_std = ws_std[pd.DataFrame(iikeep,columns=ws_std.columns)]
        ti = ti[pd.DataFrame(iikeep,columns=ti.columns)]
        wd = wd[pd.DataFrame(iikeep,columns=wd.columns)]
        cnr = cnr[pd.DataFrame(iikeep,columns=cnr.columns)]
        ava = ava[pd.DataFrame(iikeep,columns=ava.columns)]
        return ws, ws_std, ti, wd, cnr, ava        
    
    def adjust_period(self, start_idx, end_idx):
        super().adjust_period_general(start_idx, end_idx)
        self.cnr = self.cnr.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
        self.ava = self.ava.loc[start_idx[0]:end_idx[0]].reset_index(drop=True)
    

class Zephir(RSD):
    pass

            
    
if __name__ == "__main__":
    # inputs
    type_list = ["sta", "csv", "dat"]  
    
    filenames_dict = load.get_filenames_dict(type_list)    
    inputdata_raw_dict = load.load_files(filenames_dict)
    
    windcubes = {}
    for i in inputdata_raw_dict:
        if "WLS" in i: 
            windcubes[i] = Windcube(inputdata_raw_dict[i],i)
                
    masts = {}
    for i in inputdata_raw_dict:
        if i == "Mast HAW":
            masts[i] = Mast_HAW(inputdata_raw_dict[i])
        if i == "Mast Janneby":
            masts[i] = Mast_Janneby(inputdata_raw_dict[i])