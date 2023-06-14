from download_waves import generate_waves
from download_gfs import generate_gfs, get_start_date, interpolate_gfs
from download_tides import generate_tides
from msl import msl_plots
from plotter import produce_report
import pandas as pd
from maileroffice import send_email
import os
import json

#CONFIG
#now = '2023060200'
import datetime as dt
#now = dt.datetime(2023,6,8,18)

def hallsTailoredForecast(now, configfile):
    #READ Configurations
    f = open(configfile)
    data = json.load(f)
    output_directory = data[0]['Tailored_Forecast']['output_directory']
    location_of_small_map_png = data[0]['Tailored_Forecast']['location_of_small_map_png']
    mail_subject = data[0]['Tailored_Forecast']['mail_subject']
    mail_body = data[0]['Tailored_Forecast']['mail_body']
    receipients = data[0]['Tailored_Forecast']['receipients']

    matfile = '../runs/'+now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H") +'/output.mat'
    tideFolder = '../tmp/'
    gfsDataDir = '../Regional_tmp/out'
    location_pts_csv = data[0]['Tailored_Forecast']['forecast_locations_path']

    out_folder = output_directory+'/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H") 
    try:
        os.mkdir(output_directory)
    except OSError as error:
        pass

    try:
        os.mkdir(out_folder)
    except OSError as error:
        print(error)
    
    #GET Locations
    df = pd.read_csv(location_pts_csv)

    #Remove files 
    try:
        print('removing files')
        os.remove('dataset.csv')
        os.remove('gfs_3hourly_interpolated.csv')
        os.remove('gfs_3hourly.csv')
        os.remove('tide_data.csv')
        for file in os.listdir('.'):
            if file.endswith('.pdf'):
                os.remove(file) 
    except:
        pass

    document_arr = []
    for index, row in df.iterrows():
    
        tide_name = row['atoll']
        out_file = str(tide_name)+'_'+now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")+'.pdf'
        document_arr.append(out_folder+'/'+out_file)
        lon = row['lon']
        lat = row['lat']
        proper_name = tide_name.replace("_", " " )

        location_image = location_of_small_map_png+str(tide_name)+'.png'
        #DOWNLOAD WAVES
        rows_to_remove = 48
        generate_waves(lon, lat, rows_to_remove,matfile)

        #DOWNLOAD GFS
        generate_gfs(lon, lat,gfsDataDir)
        interpolate_gfs()
        start_date = get_start_date()

        ##DOWNLOAD TIDES
        generate_tides(tide_name, start_date,tideFolder)

        ##CREATE MAPS
        factor = 1 if start_date.hour == 12 else 0
        msl_plots(1,gfsDataDir)

        ####PLOTTER#####
        produce_report(tideFolder,out_file ,proper_name,lat,lon,tide_name,location_image,out_folder)
        print('Report Generated Successfully!')

    #SEND EMAIL
    send_email(document_arr, receipients,mail_subject,mail_body)
    print('Email Send Successfully!')

