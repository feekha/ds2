# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 15:42:59 2023

@author: HENROG
"""
import load
from entities import Windcube, Mast_HAW, Mast_Janneby
import numpy as np
import pandas as pd
from copy import deepcopy
from math import floor
from scipy.signal import correlate
from scipy.stats import linregress
import time
import plotly.express as px
import plotly.graph_objects as go
from dash import html, dcc
import pickle

class Verification:
    
    def __init__(self, RSD, Ref):
        # objects of RSD and Mast
        self.RSD = deepcopy(RSD)
        self.Ref = deepcopy(Ref)
        
        self.do_initial_verification()
        
        print('Reference SN: ' + self.Ref.sn)
        print('RSD SN: ' + self.RSD.sn)

    def get_comparison_levels(self):
        ws_hei = []
        self.Ref.hei = [round(h) if h%1 == 0 else h for h in self.Ref.hei]
        # exclude 29 m Janneby
        if self.Ref.sn == 'Mast Janneby':
            del self.Ref.hei[3] 
        for hei in self.Ref.hei:
            ws_hei.append((hei, min(self.RSD.hei, key = lambda x: abs(x-hei))))
            if ws_hei[-1][1]%1 == 0:
                ws_hei[-1] = (hei, round(ws_hei[-1][1]))
        self.ws_hei = ws_hei
        
        wd_hei = []
        if hasattr(self.Ref,'vane_hei'):
            for hei in self.Ref.vane_hei:
                wd_hei.append((hei,min(self.RSD.hei, key = lambda x: abs(x-hei))))
        else:
            wd_hei = ws_hei
        self.wd_hei = wd_hei

    def auto_sync(self, x_dat, jdX, ref_dat, jdRef):
        '''
        Calculates the time shift of a data series X (e.g. measured RSD wind speed) compared to a Reference dataset.
        
    
        Parameters
        ----------
        x_dat : array-like or pandas.Series
            Data from X to be compared.
        jdX : array-like or pandas.Series
            Time series from X.
        ref_dat : array-like or pandas.Series
            Data from Reference to be compared.
        jdRef : array-like or pandas.Series
            Time series from Reference.
    
        Returns
        -------
        minutes_shift : float
            Resulting time shift in minutes. If negativ: X is ahead the Reference. If positive: X is behind reference.
            --> Time shift that needs to be applied to the X time series.
        '''
        dfX = pd.DataFrame({'timestamp':jdX,'data':x_dat})
        dfRef = pd.DataFrame({'timestamp':jdRef,'data':ref_dat})
        merged = pd.merge(dfRef, dfX, on='timestamp', how='outer').sort_values('timestamp')
        merged = merged.fillna(0) 
    
        # calculate the cross-correlation of datasets, set to full, to extract idx of max correlation
        corr = correlate(merged['data_x'],merged['data_y'],mode='full',method='auto')
        self.corr = corr
        self.merged = merged
        # calculate timeshift from difference to middle (at perfect correlation max value would be in the middle of corr array)
        minutes_shift = ((len(corr)-1)/2 - np.argmax(corr))*10 # maybe exchange the 10 with a step size for timestamps, but usually its 10 minutes
        return -minutes_shift
    
    def check_date_idx(self,date,df):
        if date in df.values:
            idx = df.loc[df==date].index
        else:         
            aux = [abs((x-date).total_seconds()) for x in df]
            for x in df:
                if abs((x-date).total_seconds()) == min(aux):
                    a = x
            idx = df[df==a].index
            # idx = df.iloc[(df-date).abs().argsort()[:1]].index 
        return idx
    
    def filter_time_series(self, **buffer):
        # buffer is to not delete valid data that before time shift
        self.startdate = max(self.RSD.jdR.iloc[0], self.Ref.jdR.iloc[0]) - pd.Timedelta(**buffer)
        self.enddate = min(self.RSD.jdR.iloc[-1], self.Ref.jdR.iloc[-1]) + pd.Timedelta(**buffer)
        # print(self.startdate)
        # print(self.enddate)
        self.RSD.adjust_period(self.check_date_idx(self.startdate, self.RSD.jdR),self.check_date_idx(self.enddate, self.RSD.jdR))
        self.Ref.adjust_period(self.check_date_idx(self.startdate, self.Ref.jdR),self.check_date_idx(self.enddate, self.Ref.jdR))
    
    def get_verification_period(self):
        # pre filter datasets, because crosscorelation will falsly correlate with higher WS later or earlier
        # maybe put some linregression correlation instead: is it efficient? 
        # self.filter_time_series(days=1)
        self.filter_time_series(hours=2)
        # calculate and apply time shift
        shift = self.auto_sync(
            x_dat = self.RSD.ws[[col for col in self.RSD.ws.columns if str(self.ws_hei[0][1]) == col.split('m')[0]][0]],
            jdX = self.RSD.jdR,
            ref_dat = self.Ref.ws[[col for col in self.Ref.ws.columns if str(self.ws_hei[0][0]) == col.split('m')[0]][0]],
            jdRef = self.Ref.jdR,
            )
        self.RSD.apply_time_shift(minutes=shift)
        # filter data sets with exact start and enddate
        self.filter_time_series(minutes=0)
        print('timeshift: '+str(shift))
    
        
    def hei_from_col(self, name):
        num = ''
        for i in name:
            if i.isdigit():
                num += i
            else: 
                break
        return int(num)
    
    def get_comparison_data(self):
        h0_ws = [floor(h[0]) for h in self.ws_hei]
        h1_ws = [floor(h[1]) for h in self.ws_hei]
        self.RSD.ws = self.RSD.ws.drop(columns=[col for col in self.RSD.ws.columns if self.hei_from_col(col) not in h1_ws])
        self.Ref.ws = self.Ref.ws.drop(columns=[col for col in self.Ref.ws.columns if self.hei_from_col(col) not in h0_ws])
        self.RSD.ws_std = self.RSD.ws_std.drop(columns=[col for col in self.RSD.ws_std.columns if self.hei_from_col(col) not in h1_ws])
        self.Ref.ws_std = self.Ref.ws_std.drop(columns=[col for col in self.Ref.ws_std.columns if self.hei_from_col(col) not in h0_ws])
        self.RSD.ti = self.RSD.ti.drop(columns=[col for col in self.RSD.ti.columns if self.hei_from_col(col) not in h1_ws])
        self.Ref.ti = self.Ref.ti.drop(columns=[col for col in self.Ref.ti.columns if self.hei_from_col(col) not in h0_ws])
        
        h0_wd = [floor(h[0]) for h in self.wd_hei]
        h1_wd = [floor(h[1]) for h in self.wd_hei]
        self.RSD.wd = self.RSD.wd.drop(columns=[col for col in self.RSD.wd.columns if self.hei_from_col(col) not in h1_wd])
        self.Ref.wd = self.Ref.wd.drop(columns=[col for col in self.Ref.wd.columns if self.hei_from_col(col) not in h0_wd])
        
    def timeline_raw(self):
        max_hei_ws_rsd = self.RSD.ws[[col for col in self.RSD.ws.columns if str(self.ws_hei[0][1]) == col.split('m')[0]][0]]
        max_hei_ws_ref = self.Ref.ws[[col for col in self.Ref.ws.columns if str(self.ws_hei[0][0]) == col.split('m')[0]][0]]
        
        dfX = pd.DataFrame({'timestamp':self.RSD.jdR,'data':max_hei_ws_rsd})
        dfRef = pd.DataFrame({'timestamp':self.Ref.jdR,'data':max_hei_ws_ref})
        merged = pd.merge(dfRef, dfX, on='timestamp', how='outer').sort_values('timestamp')
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=merged['timestamp'],
                                 y=merged['data_x'],
                                 name=self.Ref.sn,
                                 line_color = '#00FF00',))
        fig.add_trace(go.Scatter(x=merged['timestamp'],
                                 y=merged['data_y'],
                                 name= self.RSD.sn,
                                 line_color = '#FF0000',))
        fig.update_layout(plot_bgcolor = '#FFFFFF',
                          font_family = 'Arial',
                          font_color = 'black',
                          title = self.Ref.sn + ' @ ' + str(self.ws_hei[0][0]) + ' m and ' + self.RSD.sn + ' @ ' + str(self.ws_hei[0][1]) + ' m',
                          title_x = 0.5,
                          legend_yanchor="top", 
                          legend_y=1.0,
                          legend_xanchor="left", 
                          legend_x=0.4,
                          legend_orientation='h',
                          margin_t = 25) 
        fig.update_xaxes(
            showgrid = True,
            gridwidth = 1,
            gridcolor = '#DDDDDD',
            showline = True,
            linecolor='#000000',
            mirror = True,       
            )
        fig.update_yaxes(
            range = [0,20],
            showgrid = True,
            gridwidth = 1,
            gridcolor = '#DDDDDD',
            showline = True,
            linecolor='#000000',
            mirror = True,
            title_text = 'WS in m/s',
            )
        
        self.timeline_plot = dcc.Graph(figure=fig)

        
        
    def alpha_coef(self,):
        pass
    
    def shear_correction(self,):
        pass
    
    def QC_filtering(self):
        if self.Ref.sn == 'Mast HAW':
            print('filter for HAW')
            self.QC_filtering_HAW()
        elif self.Ref.sn == 'Mast Janneby':
            print('filter for Janneby')
            self.QC_filtering_Janneby()
        else:
            a_ref = self.Ref.jdR.dt.strftime('%y-%m-%d-%H-%M')
            a_rsd = self.RSD.jdR.dt.strftime('%y-%m-%d-%H-%M')
            
            self.qc_idx_ref = {}
            self.qc_idx_rsd = {}
            for hei in self.ws_hei:
                idx_list = []
                ws_id = [col for col in self.Ref.ws.columns if str(hei[0]) == col.split('m')[0]][0]
                # nan in ws
                idx_list.append(self.Ref.ws[ws_id].loc[self.Ref.ws[ws_id].notnull()].index.to_list())

                ws_id_rsd = [col for col in self.RSD.ws.columns if str(hei[1]) == col.split('m')[0]][0]
                # convert rsd idx to ref idx
                self.idx_ws_nan_rsd = self.RSD.ws[ws_id_rsd].loc[self.RSD.ws[ws_id_rsd].notnull()].index.to_list()       
                b = set(self.RSD.jdR[self.idx_ws_nan_rsd].dt.strftime('%y-%m-%d-%H-%M'))
                idx_list.append([i for i, x in enumerate(a_ref) if x in b])    #!!!        
                
                valid_idx_ref =  list(set(idx_list[0]).intersection(*idx_list[1:]))
                self.qc_idx_ref[hei[0]] = valid_idx_ref
                b = set(self.Ref.jdR[valid_idx_ref].dt.strftime('%y-%m-%d-%H-%M'))
                self.qc_idx_rsd[hei[1]] = [i for i, x in enumerate(a_rsd) if x in b]
        # put some qc idx for lidar vs lidar --> filter nan and stuff
    
    def QC_filtering_HAW(self,):
        ws_thres = 3.0
        plausibility = 0.3
        rh_thres = 80.0
        temp_thres = 2.0
        
        wd_lb1=20
        wd_ub1=120
        
        # !! introduce a and b as helper var for loops !!
        # check later if multiple conditions can be summarized in one step to safe computational power
        a_ref = self.Ref.jdR.dt.strftime('%y-%m-%d-%H-%M')
        a_rsd = self.RSD.jdR.dt.strftime('%y-%m-%d-%H-%M')
        
        self.qc_idx_ref = {}
        self.qc_idx_rsd = {}
        for hei in self.ws_hei:
            
            # nearest wind vane --> should be different?
            nearest_wd = min(self.wd_hei, key = lambda x: abs(x[0]-hei[0]))
            # temp hei only the highest
            nearest_temp = self.Ref.temp_hei[0]
            nearest_humi = self.Ref.humi_hei[0]

            idx_list = []
            # wd filtering with nearest wd vane
            wd_id = [col for col in self.Ref.wd.columns if str(nearest_wd[0]) == col.split('m')[0]][0]
            idx_list.append(self.Ref.wd[wd_id].loc[(self.Ref.wd[wd_id]<=wd_lb1) | (self.Ref.wd[wd_id]>=wd_ub1)].index.to_list())

            # ws filtering
            ws_id = [col for col in self.Ref.ws.columns if str(hei[0]) == col.split('m')[0]][0]
            idx_list.append(self.Ref.ws[ws_id].loc[self.Ref.ws[ws_id] >= ws_thres].index.to_list())
            # add WS filtering for Lidar ???? or wrong in Matlab??
            # nan in ws
            idx_list.append(self.Ref.ws[ws_id].loc[self.Ref.ws[ws_id].notnull()].index.to_list())

            ws_id_rsd = [col for col in self.RSD.ws.columns if str(hei[1]) == col.split('m')[0]][0]
            # convert rsd idx to ref idx
            self.idx_ws_nan_rsd = self.RSD.ws[ws_id_rsd].loc[self.RSD.ws[ws_id_rsd].notnull()].index.to_list()       
            b = set(self.RSD.jdR[self.idx_ws_nan_rsd].dt.strftime('%y-%m-%d-%H-%M'))
            idx_list.append([i for i, x in enumerate(a_ref) if x in b])    #!!!        
            # print(time.time()-starttime)
            
            # ws difference filtering
            wsDiff_id = [col for col in self.Ref.wsDiff.columns if str(hei[0]) == col.split('m')[0]][0]
            idx_list.append(self.Ref.wsDiff[wsDiff_id].loc[abs(self.Ref.wsDiff[wsDiff_id]) <= plausibility].index.to_list())

            # icing filter
            temp_id = [col for col in self.Ref.temp.columns if str(nearest_temp) == col.split('m')[0]][0]
            humi_id = [col for col in self.Ref.humi.columns if str(nearest_humi) == col.split('m')[0]][0]
            idx_ice = []           
            for i, (temp, humi) in enumerate(zip(self.Ref.temp[temp_id],self.Ref.humi[humi_id])):
                if (temp >= temp_thres) | (humi <= rh_thres):
                    idx_ice.append(i)
            idx_list.append(idx_ice)

            valid_idx_ref =  list(set(idx_list[0]).intersection(*idx_list[1:]))
            self.qc_idx_ref[hei[0]] = valid_idx_ref
            b = set(self.Ref.jdR[valid_idx_ref].dt.strftime('%y-%m-%d-%H-%M'))
            self.qc_idx_rsd[hei[1]] = [i for i, x in enumerate(a_rsd) if x in b]

            # Leider super wenige Daten im Besipiel Datensatz, weil Icing und WS alles filtert.
        
        
    
    def QC_filtering_Janneby(self,): # combine HAW and Janneby into functions later
        ws_thres = 3.0
        plausibility = 0.3
        rh_thres = 80.0
        temp_thres = 2.0
        
        wd_lb1=50
        wd_ub1=110
        wd_lb2=130
        wd_ub2=170
        
        # !! introduce a and b as helper var for loops !!
        # check later if multiple conditions can be summarized in one step to safe computational power
        a_ref = self.Ref.jdR.dt.strftime('%y-%m-%d-%H-%M')
        a_rsd = self.RSD.jdR.dt.strftime('%y-%m-%d-%H-%M')
        
        self.idx_list = {}
        
        self.qc_idx_ref = {}
        self.qc_idx_rsd = {}
        for hei in self.ws_hei:
            # nearest wind vane --> should be different?
            nearest_wd = min(self.wd_hei, key = lambda x: abs(x[0]-hei[0]))
            print('WD bei ' + str(hei) + 'm: ' + str(nearest_wd))
            # temp hei only the 10, humi the heighest
            nearest_temp = self.Ref.temp_hei[1]
            nearest_humi = self.Ref.humi_hei[0]
            
            idx_list = []
            # wd filtering with nearest wd vane
            wd_id = [col for col in self.Ref.wd.columns if str(nearest_wd[0]) == col.split('m')[0]][0]
            idx_list.append(self.Ref.wd[wd_id].loc[(self.Ref.wd[wd_id]<=wd_lb1) | (self.Ref.wd[wd_id]>=wd_ub1)].index.to_list())          
            idx_list.append(self.Ref.wd[wd_id].loc[(self.Ref.wd[wd_id]<=wd_lb2) | (self.Ref.wd[wd_id]>=wd_ub2)].index.to_list())

            # ws filtering
            ws_id = [col for col in self.Ref.ws.columns if str(hei[0]) == col.split('m')[0]][0]
            idx_list.append(self.Ref.ws[ws_id].loc[self.Ref.ws[ws_id] >= ws_thres].index.to_list())
            # nan in ws
            idx_list.append(self.Ref.ws[ws_id].loc[self.Ref.ws[ws_id].notnull()].index.to_list())

            ws_id_rsd = [col for col in self.RSD.ws.columns if str(hei[1]) == col.split('m')[0]][0]
            # convert rsd idx to ref idx
            self.idx_ws_nan_rsd = self.RSD.ws[ws_id_rsd].loc[self.RSD.ws[ws_id_rsd].notnull()].index.to_list()       
            b = set(self.RSD.jdR[self.idx_ws_nan_rsd].dt.strftime('%y-%m-%d-%H-%M'))
            idx_list.append([i for i, x in enumerate(a_ref) if x in b])    #!!!        
            
            # ws difference filtering
            wsDiff_id = [col for col in self.Ref.wsDiff.columns if str(hei[0]) == col.split('m')[0]][0]
            idx_list.append(self.Ref.wsDiff[wsDiff_id].loc[abs(self.Ref.wsDiff[wsDiff_id]) <= plausibility].index.to_list())

            # icing filter
            temp_id = [col for col in self.Ref.temp.columns if str(nearest_temp) == col.split('m')[0]][0]
            humi_id = [col for col in self.Ref.humi.columns if str(nearest_humi) == col.split('m')[0]][0]
            idx_ice = []           
            for i, (temp, humi) in enumerate(zip(self.Ref.temp[temp_id],self.Ref.humi[humi_id])):
                if (temp >= temp_thres) | (humi <= rh_thres):
                    idx_ice.append(i)
            idx_list.append(idx_ice)

            valid_idx_ref =  list(set(idx_list[0]).intersection(*idx_list[1:]))
            self.qc_idx_ref[hei[0]] = valid_idx_ref
            b = set(self.Ref.jdR[valid_idx_ref].dt.strftime('%y-%m-%d-%H-%M'))
            self.qc_idx_rsd[hei[1]] = [i for i, x in enumerate(a_rsd) if x in b]
            
            self.idx_list[hei[0]] = idx_list
        
    
    def build_corr_fig(self,df,slope,const,R2,rsd_mean,ref_mean,diff_mean,diff_std,num,hei):
        fig = px.scatter(df,
                        x=df.columns[0],
                        y=df.columns[1],
                        title= self.Ref.sn + ' @ ' + str(hei[0]) + ' m vs. ' + self.RSD.sn + ' @ ' + str(hei[1]) + ' m'  ,
                        )
        fig.update_layout(
            plot_bgcolor='#FFFFFF',
            font_family='Arial',
            font_color='#000000',
            title_xanchor = 'center',
            title_yanchor = 'top',
            title_y = 0.9,
            title_x = 0.5,
            )
        fig.update_xaxes(
            range=[0,25],
            gridcolor='#DDDDDD',
            showline=True,
            linecolor='#000000',
            mirror=True)
        fig.update_yaxes(
            range=[0,25],
            gridcolor='#DDDDDD',
            showline=True,
            linecolor='#000000',
            mirror=True)
        fig.update_traces(marker_color='#0000FF')
        # add incline line 
        fig.add_shape(
            type='line',
            x0=0, y0=0,
            x1=25, y1=25,
            line_color='#000000',
            xref='x', yref='y'
            )
        
        # build regress line
        tl_x = [x for x in range(0,26)]
        tl_y = [x*slope+const for x in tl_x]
        fig.add_trace(go.Scatter(
            x=tl_x,
            y=tl_y,
            mode='lines',
            hovertemplate=None,
            hoverinfo='skip',
            marker_color='#00EEFF',
            showlegend=False,
            ))
      
        # # get OLS results
        # fit_results = px.get_trendline_results(fig).px_fit_results.iloc[0]
        # R2 = fit_results.rsquared
        # const = fit_results.params[0]
        # slope = fit_results.params[1]
        
        # # build new ols line
        # tl_x = [x for x in range(0,26)]
        # tl_y = [x*slope+const for x in tl_x]
        
        # # find ols line in traces of figure
        # for  k, trace  in enumerate(fig.data):
        #         if trace.mode is not None and trace.mode == 'lines':
        #             tr_line = k
        # # redifine ols line
        # # print(fig.data[tr_line])
        # fig.data[tr_line].update(
        #     x=tl_x,
        #     y=tl_y,
        #     hovertemplate=None,
        #     hoverinfo='skip',
        #     marker_color='#00EEFF')
        
        # add textbox
        fig.add_trace(go.Scatter(
            x=[1,1,1,1,1,1,1,17],
            y=[24,22.75,21.5,20.25,19,17.75,16.5,1],
            text=['Mean Reference WS: %s m/s'% round(ref_mean,2), # text als input bauen
                  'Mean RSD WS: %s m/s'% round(rsd_mean,2),
                  'Mean WS-Diff.: %s m/s (%s %%)'% (round(diff_mean,2),round(diff_mean/ref_mean*100,2)),
                  'Std. WS-Diff.: %s m/s (%s %%)'% (round(diff_std,2),round(diff_std/ref_mean*100,2)),
                  'Slope m = %s'% round(slope,4),
                  'Offset b = %s'% round(const,4),
                  'R\U000000B2 = %s'% round(R2,4),
                  '10-min-values # %s'% num],
            mode='text',
            textposition='middle right',
            showlegend=False,
            ))
        return fig
    
    def generate_table(self, header, data, classname='table'):
        return html.Table([
            html.Thead([
                html.Tr([html.Th(header[col][i]) for col in range(0,len(header))
                         ]) for i in range(0,len(header[0]))
                
            ]),
            html.Tbody([
                html.Tr([
                    html.Td(data[col][i][0], className=data[col][i][1]) for col in range(0,len(header))
                ]) for i in range(0,len(data[0]))
            ])
        ],
        className=classname)

    def WS_correlation(self,):
        
        # prepare for table as well later
        
        self.ws_corr = {}
        for hei in self.ws_hei:
            # get ws at given height with qc filtering
            rsd_ws = self.RSD.ws[[col for col in self.RSD.ws.columns if str(hei[1]) == col.split('m')[0]][0]][self.qc_idx_rsd[hei[1]]]
            ref_ws = self.Ref.ws[[col for col in self.Ref.ws.columns if str(hei[0]) == col.split('m')[0]][0]][self.qc_idx_ref[hei[0]]]
            # set date as idx to match ws
            rsd_ws.index = self.RSD.jdR[self.qc_idx_rsd[hei[1]]].to_list()
            ref_ws.index = self.Ref.jdR[self.qc_idx_ref[hei[0]]].to_list()
            
            df = pd.concat([ref_ws, rsd_ws], axis=1)
            df.columns = ['Reference Wind Speed (m/s)', 'RSD Wind Speed (m/s)']
            
            # see here for table: https://dash.plotly.com/layout
            # function above copied from there
            # later maybe put all dash graph/table generation functions into one script (out of class)
            
            header = [
                ['%s m level'% hei[0],' '],
                ['# values','-'],
                ['slope','-'],
                ['R\U000000B2','-'],
                ['WS-avg Ref','[m/s]'],
                ['WS-avg RSD','[m/s]'],
                ['mean diff.','[m/s]'],
                ['rel. mean diff.','%']
                ]
            data = [[] for i in range(8)]
            data[0] = [('All >= 3 m/s',None),('4 - 16 m/s',None)]
            for b in [False,True]:
                if b == True:
                    ref_ws = ref_ws[(ref_ws >= 3.75) & (ref_ws <= 16.25)]
                    rsd_ws = rsd_ws.loc[ref_ws.index]
                    
                    df = pd.concat([ref_ws, rsd_ws], axis=1)
                    df.columns = ['Reference Wind Speed (m/s)', 'RSD Wind Speed (m/s)']
                
                
                # RSD WS average
                rsd_mean = rsd_ws.mean()
                data[5].append((round(rsd_mean,2),None))
                
                # Reference (cup) average
                ref_mean = ref_ws.mean()
                data[4].append((round(ref_mean,2),None))
                
                # mean ws difference
                diff_mean = rsd_ws.sub(ref_ws).mean()
                data[6].append((round(diff_mean,2),None))
                
                # relative mean ws diff
                rel_diff_mean = diff_mean/ref_mean*100
                if (rel_diff_mean>=-1) & (rel_diff_mean<=1):
                    rel_diff_mean_color = 'green-cell'
                elif (rel_diff_mean>=-1.5) & (rel_diff_mean<=1.5):
                    rel_diff_mean_color = 'yellow-cell'
                else:
                    rel_diff_mean_color = 'red-cell'
                data[7].append((round(rel_diff_mean,2),rel_diff_mean_color))
                
                # std of mean ws diff
                diff_std = rsd_ws.sub(ref_ws).std()
                
                # number of values
                num = rsd_ws.sub(ref_ws).count()
                if num >= 600:
                    num_color = 'green-cell'
                else:
                    num_color = 'red-cell'
                data[1].append((num,num_color))
                
                # regression stats
                # put something here for the case that everything is nan
                # with open('rdat.pickle', 'wb') as handle:
                #     pickle.dump(df, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    
                slope, const, R2, p, std_err = linregress(df.dropna().astype(float))
                
                # slope
                if (slope>=0.98) & (slope<=1.02):
                    slope_color = 'green-cell'
                elif (slope>=0.97) & (slope<=1.03):
                    slope_color = 'yellow-cell'
                else:
                    slope_color = 'red-cell'
                data[2].append((round(slope,4),slope_color))
                    
                # R2
                if R2>=0.98:
                    R2_color = 'green-cell'
                elif R2>=0.97:
                    R2_color = 'yellow-cell'
                else:
                    R2_color = 'red-cell'
                data[3].append((round(R2,4),R2_color))
                
                if b == False:
                    fig = self.build_corr_fig(df,slope,const,R2,rsd_mean,ref_mean,diff_mean,diff_std,num,hei)
                    
            tab = self.generate_table(header, data)
            # also save R2 slope and const into dict here
            self.ws_corr[hei[0]] = {}
            self.ws_corr[hei[0]]['table'] = tab
            self.ws_corr[hei[0]]['figure'] = fig
        
    
    def TI_correlation(self,):
        pass
    
    def WD_correlation(self,):
        pass
    
    
    def bin_df(self, df, bins):
        
        if type(df) == type(pd.Series(0)):
            df = pd.DataFrame(df,columns=[df.name])
        
        binned_dfs = []
        for col in df.columns:
            df[str(col)+'_binned'] = pd.cut(df[col], bins=bins) # Use cut function to create the bins
            df_binned = df[[col,str(col)+'_binned']].groupby(str(col)+'_binned').count()
            binned_dfs.append(df_binned)

        df_binned = pd.concat(binned_dfs, axis = 1)
        return df_binned
    
    def WS_bins(self):
        # Create bins
        bins = pd.interval_range(start=3.75, end=16.25, freq=0.5, closed='right') 
        
        header = [[' ']]+[[str(b.mid)] for b in bins]
        data = [[] for i in range(len(bins)+1)]
        
        for hei in self.ws_hei:
            # get ws at given height with qc filtering
            # rsd_ws = self.RSD.ws[[col for col in self.RSD.ws.columns if str(hei[1]) == col.split('m')[0]][0]][self.qc_idx_rsd[hei[1]]]
            ref_ws = self.Ref.ws[[col for col in self.Ref.ws.columns if str(hei[0]) == col.split('m')[0]][0]][self.qc_idx_ref[hei[0]]]
            ws_ref_binned = self.bin_df(ref_ws,bins)
            
            data[0].append(('# values at '+str(hei[0])+' m',None))
            i = 1
            for val in ws_ref_binned.iloc[:,0]:
                if val >= 3:
                    color=None
                else:
                    color='red-cell'
                data[i].append((val,color))
                i += 1
        
        self.ws_ref_binned_table = self.generate_table(header, data, classname='table iec-bins')
        
    
    def absolute_error(self,):
        pass
    
    def Availability(self):
        # max number of 10-Min steps
        max_num = ((self.enddate-self.startdate).total_seconds()+600)/600
        # self.max_num_corr = max_num
        self.sys_ava_num = len(self.RSD.jdR.loc[~self.RSD.jdR.isin(list(self.outages_unique))])
        self.sys_ava = self.sys_ava_num/self.max_num_corr
        rsd_ws_valid = self.RSD.ws.loc[~self.RSD.jdR.isin(list(self.outages_unique))].count()
        
        header_ava = [['Height']]
        header_num = [['WS-Range']]
        
        data_ava = [[]]
        data_ava[0] = [('Max # 10-min points',None),
                       ('After accounting outages',None),
                   ('Data present',None),
                   ('System availability',None),
                   ('Total # 10-min valid data',None),
                   ('Data availability',None),
                   ('# after external filtering',None),
                   ('Data availability for comparison',None)
                   ]
        
        data_num = [[]]
        data_num[0] = [('All >= 3 m/s',None),
                       ('4-8 m/s',None),
                       ('8-12 m/s',None),
                       ('4-16 m/s',None)
                       ]
        
        for hei in self.ws_hei:
            header_ava.append(['%s m'% hei[1]])
            header_num.append(['%s m'% hei[1]])
            rsd_ws_valid_num = rsd_ws_valid[[col for col in self.RSD.ws.columns if str(hei[1]) == col.split('m')[0]][0]]
            valid_sys_ava = rsd_ws_valid_num/self.max_num_corr
            valid_comp_num = len(self.qc_idx_rsd[hei[1]])
            valid_comp_ava = valid_comp_num/self.max_num_corr
                       
            data_ava.append([(max_num,None),
                             (self.max_num_corr,None),
                       (self.sys_ava_num,None),
                       ('%s %%'% (round(self.sys_ava*100,1)),'red-cell' if self.sys_ava < 0.95 else 'green-cell'),
                       (rsd_ws_valid_num,None),
                       ('%s %%'% (round(valid_sys_ava*100,1)),'red-cell' if valid_sys_ava < 0.9 else 'green-cell'),
                       (valid_comp_num,None),
                       ('%s %%'% (round(valid_comp_ava*100,1)),None)
                       ])
        
            rsd_ws = self.Ref.ws[[col for col in self.Ref.ws.columns if str(hei[0]) == col.split('m')[0]][0]][self.qc_idx_ref[hei[0]]]
            rsd_ws_4_8 = rsd_ws[(rsd_ws >= 4) & (rsd_ws < 8)]
            rsd_ws_8_12 = rsd_ws[(rsd_ws >= 8) & (rsd_ws < 12)]
            rsd_ws_4_16 = rsd_ws[(rsd_ws >= 3.75) & (rsd_ws < 16.25)]
            
            data_num.append([(valid_comp_num,'red-cell' if valid_comp_num < 600 else 'green-cell'),
                             (len(rsd_ws_4_8),'red-cell' if len(rsd_ws_4_8) < 200 else 'green-cell'),
                             (len(rsd_ws_8_12),'red-cell' if len(rsd_ws_8_12) < 200 else 'green-cell'),
                             (len(rsd_ws_4_16),'red-cell' if len(rsd_ws_4_16) < 600 else 'green-cell')
                             ])
        
        self.table_availabilities = self.generate_table(header_ava, data_ava, classname='table ava') # create new css for availability tables
        self.table_num_datapoints = self.generate_table(header_num, data_num, classname='table ava')
    
    
    
    def find_outages(self,jdR):
        outages = []
        for t in range(1,len(jdR)):                
            if jdR[t]-jdR[t-1] != pd.Timedelta(10,'minutes'):
                num_10_min = int((jdR[t]-jdR[t-1])/pd.Timedelta(10,'minutes'))-1
                outages.append('From ' + str(jdR[t-1]+pd.Timedelta(10,'minutes')) + ' to ' + str(jdR[t]-pd.Timedelta(10,'minutes')) + ' ......' + str(num_10_min))
        return outages
    
    def timestamp_outages(self, outages_valid):
        all_timestamps = []
        if outages_valid == None:
            return all_timestamps
        else:
            for outage in outages_valid:
                outage = str(outage)
                starttime = ' '.join(outage.split()[1:3])
                pd.to_datetime(starttime)
                endtime =  ' '.join(outage.split()[4:6])
                pd.to_datetime(endtime)
                all_timestamps.append((starttime, endtime))
            return all_timestamps
    
    def apply_outages(self, rsd_outages_str, ref_outages_str):
        rsd_outages_valid = self.timestamp_outages(rsd_outages_str)
        ref_outages_valid = self.timestamp_outages(ref_outages_str)
        all_outage_timestamps = []
        for outage in rsd_outages_valid + ref_outages_valid:
            all_outage_timestamps = all_outage_timestamps + (
                                                            pd.date_range(
                                                                start=outage[0],
                                                                end=outage[1],
                                                                freq='10T').to_list())
        outages_unique = np.unique(all_outage_timestamps)
        max_num = ((self.enddate-self.startdate).total_seconds()+600)/600
        self.max_num_corr = max_num - len(outages_unique)
        self.outages_unique = outages_unique
        self.Availability()

    def do_initial_verification(self):
        # starttime = time.time()
        self.get_comparison_levels()
        # print('get com: '+str(time.time()-starttime))
        self.get_verification_period()
        # print('get period: '+str(time.time()-starttime))
        self.get_comparison_data()
        # print('get data: '+str(time.time()-starttime))
        self.timeline_raw()
        
        self.QC_filtering()
        # print('QC filter: '+str(time.time()-starttime))
        self.WS_correlation()
        # print('WS corr: '+str(time.time()-starttime))
        self.WS_bins()
        
        self.apply_outages([], [])
        self.rsd_outages = self.find_outages(self.RSD.jdR)
        self.ref_outages = self.find_outages(self.Ref.jdR)
        self.Availability()



if __name__ == "__main__":
    # inputs
    type_list = ["sta", "csv", "dat"]  
    starttime = time.time()
    filenames_dict = load.get_filenames_dict(type_list)    
    inputdata_raw_dict = load.load_files(filenames_dict)
    
    RSDs = {}
    for i in inputdata_raw_dict:
        if "WLS" in i: 
            RSDs[i] = Windcube(inputdata_raw_dict[i],i)
                
    masts = {}
    for i in inputdata_raw_dict:
        if i == "Mast HAW":
            masts[i] = Mast_HAW(inputdata_raw_dict[i])
        if i == "Mast Janneby":
            masts[i] = Mast_Janneby(inputdata_raw_dict[i])
    print(time.time()-starttime)
    
    starttime1 = time.time()
    cV = Verification(RSDs['WLS71160'],masts['Mast HAW'])
    print(time.time()-starttime1)
            
    # minutes_shift = auto_sync(windcubes['WLS7-346'].ws.iloc[:,7],windcubes['WLS7-346'].jdR,
    #                              masts['Mast HAW'].WS_Mast.iloc[:,0],masts['Mast HAW'].jdM)
    
    # windcubes['WLS7-346'].apply_time_shift(minutes=minutes_shift)