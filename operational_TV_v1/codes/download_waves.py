import scipy.io as sio
import numpy as np
import datetime
import pandas as pd

    
def cm2inch(*tupl):
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple()

def read_data_into_var(st_time,end_time,dt,var_name,mat_data):
    def time_it(start, finish, time_delta):
        while finish >= start:
            yield start
            start = start + datetime.timedelta(hours=time_delta)
    
    tt = []
    data_out = []
    for hour in time_it(st_time, end_time, dt):
        tt.append(hour)
        swan_hour = hour.strftime('%Y%m%d_%H%M%S')
        full_var_name = var_name + '_' + swan_hour
        if hour == st_time:
            data_out = mat_data[full_var_name]
        else:
            z = mat_data[full_var_name]
            data_out = np.vstack((data_out,z))
            
    return(data_out,tt)

def extract_point_ts(data_array,target_lon,target_lat,lon_array,lat_array):
    def find_coordinate_ix(target_lon,target_lat,lon_array,lat_array):
        import geopy.distance
        coords_1 = (target_lat, target_lon)
        dist = geopy.distance.distance(coords_1, (lat_array[0],lon_array[0])).km
        for i in range(1,len(lon_array)):
            dist = np.append(dist,geopy.distance.distance(coords_1, (lat_array[i],lon_array[i])).km)
        idx = dist.argmin()
        return(idx,lon_array[idx],lat_array[idx])
    
    idx,model_lon,model_lat = find_coordinate_ix(target_lon,target_lat,lon_array,lat_array)
    out_data = data_array[:,idx]
    return(out_data,idx,model_lon,model_lat)

def timeConverter(first_date):
    dt = first_date.split("_")
    f_year = dt[1][:4]
    f_month = dt[1][4:6]
    f_day = dt[1][6:8]
    f_hr = dt[2][:2]
    st_time = datetime.datetime(int(f_year),int(f_month),int(f_day),int(f_hr))
    return st_time

def getDates(myarr):
    my_arr = []
    for x in myarr:
        if "Hsig" in x[0]:
            my_arr.append(x[0])
    first_date =my_arr[0]
    last_date = my_arr[-1]
    return timeConverter(first_date), timeConverter(last_date)


def generate_waves(target_lon, target_lat,rows_to_remove,matfile):
    #PARAMS
    fl_name = matfile
    var_name = 'Hsig'
    var_name2 = 'Dir'
    dt = 1 # time step in hours

    #MAIN
    mat_data = sio.loadmat(fl_name,struct_as_record=False, squeeze_me=True)
    lon = mat_data['Xp']
    lat = mat_data['Yp']
    myarr = sio.whosmat(fl_name)
    #"""
    start_date, end_date = getDates(myarr)

    #GET HS
    Hs, tt = read_data_into_var(start_date, end_date, dt, var_name,mat_data)
    Hs_ts, idx, model_lon,model_lat = extract_point_ts(Hs,target_lon,target_lat,lon,lat)

    #GET DIR
    Dir, tt = read_data_into_var(start_date, end_date, dt, var_name2,mat_data)
    Dir_ts, idx, model_lon,model_lat = extract_point_ts(Dir,target_lon,target_lat,lon,lat)

    #Wind X
    windx, tt = read_data_into_var(start_date, end_date, dt, 'Windv_x',mat_data)
    windxS, idx, model_lon,model_lat = extract_point_ts(windx,target_lon,target_lat,lon,lat)

    #Wind Y
    windy, tt = read_data_into_var(start_date, end_date, dt, 'Windv_y',mat_data)
    windyS, idx, model_lon,model_lat = extract_point_ts(windy,target_lon,target_lat,lon,lat)

    #peakdir
    peakdir, tt = read_data_into_var(start_date, end_date, dt, 'PkDir',mat_data)
    preakdirS, idx, model_lon,model_lat = extract_point_ts(peakdir,target_lon,target_lat,lon,lat)


    #peak period
    RTpeak, tt = read_data_into_var(start_date, end_date, dt, 'RTpeak',mat_data)
    RTpeakS, idx, model_lon,model_lat = extract_point_ts(RTpeak,target_lon,target_lat,lon,lat)

    #MEAD period
    meanPeriod, tt = read_data_into_var(start_date, end_date, dt, 'Tm02',mat_data)
    meanPeriodS, idx, model_lon,model_lat = extract_point_ts(meanPeriod,target_lon,target_lat,lon,lat)


    df = pd.DataFrame({'DateTime':tt, 'Hs':Hs_ts, 'Dir':Dir_ts, 'Windx':windxS, 'Windy':windyS, 'peakdir':preakdirS, 'peakperiod':RTpeakS, 'meanperiod':meanPeriodS})

    df2 = df.iloc[rows_to_remove:]
    df2.to_csv('dataset.csv', index=None, header=True)
    return None