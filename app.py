# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 13:05:03 2022

@author: HENROG
"""

import load
import sftp_connect
from verification import Verification
from entities import Windcube, Mast_HAW, Mast_Janneby
from weather import get_weatherdata, daily_mean # import weather data funtions
import numpy as np
import pandas as pd
from copy import deepcopy
from dash import html, Dash, Output, Input, dcc, State, ctx
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pickle
from numbers import Number
import time

def dump_data(df, filename):
    '''
    Dumps excel-df as pickle file --> cache
    '''
    with open('cache/%s.pickle'% filename, 'wb') as handle:
        pickle.dump(df, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_data(filename):
    '''
    Loads excel-df from cache
    '''
    with open('cache/%s.pickle'% filename, 'rb') as handle:
        df = pickle.load(handle)
    return df
        
def find_nearest_mast(sn):
    RSDs = load_data('RSDs')
    masts = load_data('masts')
    rsd_coord = RSDs[sn].coordinates
    mast_coords = [(m,masts[m].coordinates) for m in masts]
    mast_coords.sort(key=lambda tup: abs(tup[1]['Lat']-rsd_coord['Lat'])+abs(tup[1]['Long']-rsd_coord['Long']), reverse=False)
    mast_sn = mast_coords[0][0]
    return mast_sn    

def build_WS_tabs(cV):
    return dcc.Tabs([
                dcc.Tab(
                    label = str(hei[0])+' m',
                    children = [
                        html.Br(),
                        html.P('Wind Speed Regression: '),
                        dcc.Graph(figure=cV.ws_corr[hei[0]]['figure'],className='dcc-graph'),
                        cV.ws_corr[hei[0]]['table']
                        ]
                    ) for hei in cV.ws_hei
                ],
            className='Tabs')

def get_wd_sign(wd_list):
    '''
    Bins wind directions
    '''
    wds = []
    dirs = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'] # provides order for windrose
    for wd in wd_list:
        if wd > 337:
            wds.append('\U00002B61') # N
        elif wd > 292:
            wds.append('\U00002B66') # NW
        elif wd > 247:
            wds.append('\U00002B60') # W
        elif wd > 202:
            wds.append('\U00002B69') # SW
        elif wd > 157:
            wds.append('\U00002B63') # S
        elif wd > 112:
            wds.append('\U00002B68') # SE
        elif wd > 67:
            wds.append('\U00002B62') # E
        elif wd > 22:
            wds.append('\U00002B67') # NE
        elif wd <= 22:
            wds.append('\U00002B61') # N
        elif wd in dirs:
            wds.append(wd)
    return wds, dirs

def build_weather_figure(location,weather_df):
    '''
    Builds weather forecast figure 
    '''
    condition_icons = {
        'dry-cloudy': '\U00002601',
        'dry-partly-cloudy-day': '\U0001F324',
        'dry-partly-cloudy-night': '\U0000263E'+'\U00002601',
        'rain-rain': '\U0001F327',
        'rain-cloudy': '\U0001F327',
        'rain-partly-cloudy-day': '\U0001F326',
        'rain-partly-cloudy-night': '\U0000263E'+'\U00002601',
        'fog-fog': '\U0001F32B',
        'sleet-sleet': '\U0001F327'+'Hagel',
        'snow-snow': '\U0001F327',
        'rain-wind': '\U0001F32C'+'\U000026C6',
        'dry-wind': '\U0001F32C',
        'dry-clear-day': '\U0000263C',
        'dry-clear-night': '\U0000263E',
        'hail-hail': 'Hagel',
        'thunderstorm-thunderstorm': '\U000026C8',
        }
    

    daily_weather_df = daily_mean(weather_df, 'temperature', 'wind_direction') # calcuate daily mean for given paramemeters
    
    days = []
    for t in weather_df['timestamp']:
        if str(t).split(' ')[0] not in days:
            days.append(str(t).split(' ')[0])
    day_12h = [day+' 12:00:00+00:00' for day in days]

    wd_12h = [weather_df[weather_df['timestamp'] == t]['wind_direction'].iloc[0] for t in day_12h]

    wd_signs = get_wd_sign(wd_12h)[0]
    
    
    cond_text = []
    temp_text = []
    # create text to display and append to list
    for t in day_12h:
        cond = weather_df[weather_df['timestamp'] == t]['condition'].iloc[0]
        icon = weather_df[weather_df['timestamp'] == t]['icon'].iloc[0]
        # put something here for the case of None in cond or icon
        cond_icon = cond + '-' + icon
        if cond_icon in condition_icons:
            # t_md = dcc.Markdown(f'''#### {condition_icons[cond_icon]} ''') # Markdown gives different size
            cond_text.append(f"{condition_icons[cond_icon]}") # Use css size here          
        else:
            cond_text.append(f"{cond} {icon}") # when no condition is given
        temp_text.append(f"{round(weather_df[weather_df['timestamp'] == t]['temperature'].iloc[0])} °C")
    
    weather_figure = make_subplots(specs=[[{'secondary_y': True}]]) # to get axis on left and right side
    weather_figure.add_trace(go.Scatter(
        x=day_12h,
        y=[weather_df['wind_speed'].max()+1 for x in wd_signs],
        text=wd_signs,
        mode='text',
        textposition='middle center',
        name='WD',
        showlegend=False,
        hoverinfo='skip',
        textfont_size=30,
        ))
    weather_figure.add_trace(go.Scatter(
        x=day_12h,
        y=[weather_df['wind_speed'].max()+2 for x in wd_signs],
        text=cond_text,
        mode='text',
        textposition='middle center',
        name='Condition',
        showlegend=False,
        hoverinfo='skip',
        textfont_size=30,
        ))
    weather_figure.add_trace(go.Scatter(
        x=day_12h,
        y=[0 for x in wd_signs],
        text=temp_text,
        mode='text',
        textposition='middle center',
        name='Temperature',
        showlegend=False,
        hoverinfo='skip',
        textfont_size=14,
        textfont_family = 'Arial Black',
        ))
    weather_figure.add_trace(
        go.Scatter(
            x = weather_df['timestamp'],
            y = np.round(weather_df['wind_speed'],2),
            name = 'WS',
            mode = 'lines',
            line_color = 'rgb(153,217,240)',
            ),
        secondary_y = False, # y axis on left side
        )
    # weather_figure.add_trace(
    #     go.Scatter(
    #         x = weather_df['timestamp'],
    #         y = np.round(weather_df['wind_speed_100m'],2),
    #         name = 'WS 100m',
    #         mode = 'lines',
    #         line_color = 'rgb(15,32,75)',
    #         ),
    #     )
    weather_figure.add_trace(
        go.Scatter(
            x = weather_df['timestamp'],
            y = weather_df['precipitation'].multiply(100), # to get %
            name = 'Nd.',
            mode = 'lines',
            line_color = 'rgb(152,143,134)',
            line_dash = 'dot',
            ),
        secondary_y = True, # y axis on right side
        )
    weather_figure.update_layout(
        plot_bgcolor = 'rgba(0,0,0,0)',
        font_family = 'Arial',
        font_color = 'black',
        grid_ygap = 0.1,
        title = str(location) + ' Weather Forecast',
        title_font_size = 20,
        title_x = 0.5,
        title_xanchor = 'center',
        )
    weather_figure.update_xaxes(
        dtick = 'D1', # daily ticks on x axis
        showgrid = True,
        gridwidth = 1,
        gridcolor = 'rgb(204,204,204)',
        )
    weather_figure.update_yaxes(
        title_text = 'WS in m/s',
        showgrid = False,
        side = 'right',
        secondary_y = False # y axis on left side
        )  
    weather_figure.update_yaxes(
        range = [0,100], # 0 - 100 % precipitation
        title_text = 'Precipitaion in %',
        showgrid = False,
        side = 'left',
        secondary_y = True # y axis on right side
        )       

    # weather_figure.add_layout_image( # not working, because picture must be saved in internet, not local
    #     dict(
    #         source=app.get_asset_url('kompass.png'),
    #         x=1,
    #         y=0.2,)
    #     )
    return weather_figure

def build_weather_div(weather_figs):
    weather_div_children = []
    for fig in weather_figs:
        weather_div_children.append(dcc.Graph(figure=fig,className='six columns blue-border dcc-weather-graph'))
    return weather_div_children


app = Dash(__name__,suppress_callback_exceptions=True)

app.layout = html.Div(
    [
     # header 
     html.Div(className='row',
              children=[
                  html.Div(className='two columns',
                           children=[
                                   html.Img(className="Logo",
                                            src=app.get_asset_url('DNV_Logo.png')                                       
                                   )
                               ] 
                           ),   
                  html.Div(className='eight columns',
                           children=[
                               html.H4('LiDAR ViEW', className='app__header__title'),
                               html.P('Choose a device and a reference')
                               ],                           
                           ), 
                  html.Div(className='two columns Selector',
                           children=[
                               dcc.Dropdown(
                                   options=[],
                                   value=None,
                                   placeholder='Select LiDAR',
                                   id='RSD-selector',
                                   clearable=False,
                                   ),
                               dcc.Dropdown(
                                   options=[],
                                   value=None,
                                   id='Ref-selector',
                                   placeholder='Select Reference',
                                   clearable=False,
                                   )
                               ]
                           ),
                   dcc.Interval(
                       id='daily-interval',
                       interval = 1000*86400,
                       )
                  ]
              ),
     html.Div(className='row',
              #style={'backgroundColor': '#0F204B'},
             children=[
                    html.Div(className='six columns blue-border',
                             # style={'border': '2px solid black'},
                             children=[
                                 dcc.Tabs([
                                         dcc.Tab(
                                             label = 'Availability',
                                             children = [
                                                 html.Br(),
                                                 html.P('LiDAR Availability Assessment'),
                                                 html.Div(id='table-ava',
                                                          children = [html.P('-')]
                                                     ),
                                                 html.Br(),
                                                 html.P('Number of Data Points'),
                                                 html.Div(id='table-num-datapoints',
                                                          children = [html.P('-')]
                                                     ),
                                                 html.Br(),
                                                 html.P('Timeline Wind Speed'),
                                                 html.Div(id='graph-timeline',
                                                          children = [html.P('-')]
                                                          )
                                                 ]
                                             ),
                                         dcc.Tab(
                                             label = 'Outages',
                                             children = [
                                                 html.Br(),
                                                 html.P('LiDAR Outages'),
                                                 dcc.Checklist([], id='rsd-outage-checklist', className='Checklist', labelClassName='Checklist-label',
                                                               labelStyle = {'display': 'inline-block',}),
                                                 html.Br(),
                                                 html.P('Mast Outages'),
                                                 dcc.Checklist([], id='ref-outage-checklist', className='Checklist', labelClassName='Checklist-label',
                                                               labelStyle = {'display': 'inline-block',}),
                                                 html.Br(),
                                                 html.Button('Apply', id='outage-apply-button', n_clicks=0, className='Button')
                                                 ]
                                             )
                                     ], className = 'Tabs')
                                 ]
                                 ),
                    html.Div(className='six columns blue-border',
                             # style={'border': '2px solid black'},
                             children=[
                                 html.Div(id='tabs-ws-corr',
                                          children = [html.P('-')]
                                     ),
                                 ]
                             )
                        ]
           ),
     html.Div(className='row',
              #style={'backgroundColor': '#0F204B'},
              children=[
                html.Div(className='twelve columns blue-border',
                         # style={'border': '2px solid black'},
                         children=[
                             html.Br(),
                             html.P('IEC Completeness: WS bins [m/s]'),
                             html.Div(id='table-ws-binned',
                                      children = [html.P('-')]
                                 )
                             ]
                         ),
                    ]
           ),
     html.Div(className='row',
              children=[],
              id='weather-plots'
                      )        
     ])

@app.callback(
    Output('RSD-selector','options'),
    Output('RSD-selector','value'),
    Output('Ref-selector','options'),
    Output('weather-plots','children'),
    Input('daily-interval','n_intervals')
    )
def daily_update(n_intervals):
    print('I got here 1')
    # get new data
    type_list = ["sta", "csv", "dat"]  

    # filenames_dict = load.get_filenames_dict(type_list)    
    # inputdata_raw_dict = load.load_files(filenames_dict)
    
    inputdata_raw_dict = sftp_connect.get_sftp()

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
            
    rsd_options = [key for key in RSDs.keys()]
    rsd_value = [key for key in RSDs.keys()][0]
    
    ref_options = [key for key in masts.keys()]+[key for key in RSDs.keys()]
        
    # get weather data
    all_weather_dict = {}
    for key, val in masts.items():
        all_weather_dict.update(get_weatherdata(loc=val.sn, lat=val.coordinates['Lat'], long=val.coordinates['Long']))
        all_weather_dict[key]['wind_speed_100m'] = all_weather_dict[key]['wind_speed'].multiply(1.5)
    
    # get daily weather data
    weather_figs = []
    for key, val in all_weather_dict.items():
        # if val == None:
        #     continue     
        weather_fig = build_weather_figure(key, val)
        weather_figs.append(weather_fig)
    weather_div = build_weather_div(weather_figs)
    
    dump_data(RSDs, 'RSDs')
    dump_data(masts, 'masts')
    
    return rsd_options, rsd_value, ref_options, weather_div

@app.callback(
    Output('Ref-selector','value'),
    Input('RSD-selector','value')
    )
def auto_select_reference(rsd_sn):
    print('selecting reference:')
    print(rsd_sn)
    if rsd_sn == None:
        raise PreventUpdate
    else:
        mast_sn = find_nearest_mast(rsd_sn)
        return mast_sn
 

@app.callback(
    Output('tabs-ws-corr', 'children'),
    Output('table-ws-binned', 'children'),
    Output('table-ava', 'children'),
    Output('table-num-datapoints', 'children'),
    Output('graph-timeline', 'children'),
    Output('rsd-outage-checklist', 'options'),
    Output('ref-outage-checklist', 'options'),
    Output('outage-apply-button', 'n_clicks'),
    Input('Ref-selector','value'),
    Input('outage-apply-button', 'n_clicks'),
    State('RSD-selector','value'),
    State('rsd-outage-checklist', 'value'),
    State('ref-outage-checklist', 'value')
    )
def new_verification(ref_sn, n_clicks_outages, rsd_sn, rsd_outages, ref_outages):
    print('I got here 3')
    trigger = ctx.triggered_id
    # print(trigger)
    # print(type(trigger))
    if rsd_sn == None or ref_sn == None:
        raise PreventUpdate
          
    elif trigger != 'outage-apply-button':
        RSDs = load_data('RSDs')
        masts = load_data('masts')
        
        # initialize Verification object that calculates all important stuff
        RSD = RSDs[rsd_sn]
        if ref_sn in masts:
            Ref = masts[ref_sn]
        elif ref_sn in RSDs:
            Ref = RSDs[ref_sn]
        currentVerification = Verification(RSD, Ref)
        
        # dump it to pickle ??
        with open('Verification.pickle', 'wb') as handle:
            pickle.dump(currentVerification, handle, protocol=pickle.HIGHEST_PROTOCOL)
        # or just create all plots from here?
        # what kind of filtering can be applied to verification?
        #   - defining power outages
        #   - apply shear correction manually
        #   - re-define startdate and enddate
        #   - sync datasets manually when auto sync fails
        #   - re-select verification heights?
        #   - choice of icing filter
        
        
    elif (trigger == 'outage-apply-button') & (n_clicks_outages > 0):
        with open('Verification.pickle', 'rb') as handle:
            currentVerification = pickle.load(handle)
            
        currentVerification.apply_outages(rsd_outages, ref_outages)
        
    else:
        raise PreventUpdate
 
    ws_corr_tabs = build_WS_tabs(currentVerification)
    ws_binned_table = currentVerification.ws_ref_binned_table
    
    ava_table = currentVerification.table_availabilities
    num_datapoints_table = currentVerification.table_num_datapoints
    timeline_ws = currentVerification.timeline_plot
    
    rsd_total_outages = currentVerification.rsd_outages
    ref_total_outages = currentVerification.ref_outages
    
    n_clicks_outages = 0
    
    return ws_corr_tabs, ws_binned_table, ava_table, num_datapoints_table, timeline_ws, rsd_total_outages, ref_total_outages, n_clicks_outages
    
    
    



if __name__ == "__main__":
    app.run_server(debug=True)