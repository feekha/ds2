# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 11:50:27 2023

@author: HENROG
"""
from datetime import date, timedelta
import requests
import pandas as pd
import numpy as np
from scipy.stats import circmean
import pickle

# def get_locs_lats_lons(df, all_WT=False):
#     '''
#     Get locations, latitudes and longitudes only when all parameters exist.
#     '''
#     locs = df['Messort']
#     lats = df['Breitengrad']
#     lons = df['Längengrad']
    
#     for i in [locs, lats, lons]:
#         idx = i[i.isna()].index
#         locs = locs.drop(idx, errors='ignore')
#         lats = lats.drop(idx, errors='ignore')
#         lons = lons.drop(idx, errors='ignore')
    
#     # drop duplicate locations
#     if all_WT == False:        
#         locs = locs.drop_duplicates(keep='first')
#         lats = lats[locs.index]
#         lons = lons[locs.index]
#     return locs, lats, lons

def get_weatherdata(loc, lat, long):
    '''
    Main function to get weather data from Brightsky API.
    '''
    
    first_date = date.today() # define start date
    last_date = date.today() + timedelta(days=10) # define forecast
    all_weather_dict = {}

    # check for lat and lon above 90 / 180
    # maybe insert every turbine soon
    
    # request data
    url = 'https://api.brightsky.dev/weather?date={}&last_date={}&lat={}&lon={}'.format(first_date, last_date, float(lat), float(long))
    data = requests.get(url)
    J = data.json()
    all_data = []
    if 'weather' not in J: # check for successful data request
        return None
    columns = list(J['weather'][0].keys())
    if 'fallback_source_ids' in columns:
        columns.remove('fallback_source_ids') # remove column that often gives error
        
    for t in J['weather']:
        if 'fallback_source_ids' in t:
            del t['fallback_source_ids'] # remove data of removed column 
        all_data.append(list(t.values()))        
    all_data_df = pd.DataFrame(all_data,columns=columns) # create dataframe
    
    # get distance to source and measurement heights (geographical) additionally
    meas_heights = []
    distance_to_source = []
    for sid in all_data_df['source_id']:
        for source in J['sources']:
            if source['id'] == sid:
                distance = source['distance']
                height = source['height']
                break
        meas_heights.append(height)
        distance_to_source.append(distance)
    
    # modify columns in dataframe
    all_data_df['wind_speed'] = all_data_df['wind_speed'].div(3.6)    
    all_data_df['precipitation'] = all_data_df['precipitation'].clip(upper=1)        
    all_data_df['distance_to_source'] = distance_to_source
    all_data_df['height_of_source'] = meas_heights
    all_data_df['timestamp'] = pd.to_datetime(all_data_df['timestamp'])
    all_data_df['day'] = all_data_df['timestamp'].dt.date
    all_data_df = all_data_df.drop(all_data_df.index[-1])
    all_weather_dict[loc] = all_data_df
    return all_weather_dict

def daily_mean(df, *args):
    '''
    Create a daily dataframe with given parameters.
    '''
    daily_weather = df.groupby('day', as_index=False)[[*args]].mean()
    if 'wind_direction' in args:
        daily_wd = []
        for day in df['day'].unique():
            df_day=df[df['day'] == day]
            daily_wd.append(np.rad2deg(circmean(np.deg2rad(df_day['wind_direction']))))
        daily_weather['wind_direction'] = daily_wd
    return daily_weather

if __name__ == "__main__":

    # all_weather_dict = get_weatherdata('Mast Janneby', 54.637616, 9.316778 )
    all_weather_dict1 = get_weatherdata('Mast HAW', 53.47098, 10.20683 )
    
    # with open('weather_data.pickle', 'wb') as handle:
    #     pickle.dump(all_weather_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
        
    

          

