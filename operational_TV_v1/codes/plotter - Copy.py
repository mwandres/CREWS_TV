import numpy as np
from matplotlib import pyplot as plt
import datetime
import math
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
import matplotlib as mpl
from datetime import timedelta
from download_tides import generate_tides
from download_gfs import generate_gfs, get_start_date
import locale
from datetime import timezone
locale.setlocale(locale.LC_ALL,'en_US')
plt.rcParams['axes.xmargin'] = 0

#FUNCTIONS
def unique(list1):
    # initialize a null list
    unique_list = []
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
 
def getvalue(num):
    value = 0
    if num == 0:
        value = 1000
    elif num == 1:
        value = 4000
    elif num == 2:
        value = 10000
    elif num == 3:
        value = 20000
    else:
        value = 1000000
    return value

def getColor(num):
    color = 'red'
    if num == 1000:
        color = '#BB3F3F'
    elif num == 4000:
        color = '#DB804E'
    elif num == 10000:
        color = '#EEDC5B'
    elif num == 20000:
        color = '#49759C'
    else:
        color = '#74A662'
    return color

def inter_from_256(x):
    return np.interp(x=x,xp=[0,255],fp=[0,1])

def setupAxis(ax):
    xax = ax.get_xaxis()
    months = mdates.HourLocator(byhour=[0, 12])
    months_fmt = mdates.DateFormatter('%a %d\n %H')
    xax.set_major_locator(months)
    xax.set_major_formatter(months_fmt)
    ax.xaxis.grid(True, which='major',linestyle='--')
    ax.yaxis.grid(True, which='major',linestyle='dashdot', color='#dfdfdf')
    ax.set_axisbelow(True)
    ax.tick_params(axis='both', which='both', labelsize=8)

def is_multiple_of_three(n):
    if n % 6 == 0:
        return True
    else:
        return False
   
def cm2inch(*tupl):
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple()

def get_color_from_jet(speed, min, max):
    # Normalize the speed values to the range [0, 1]
    speed_norm = (speed - min) / (max - min)

    # Get the RGBA color from the jet colormap based on the normalized speed value
    cmap = plt.get_cmap('jet')
    rgba = cmap(speed_norm, alpha=0.4)

    return rgba

def windcolor(speed):
    color = 'white'
    if float(speed) <= float(12):
        color = 'lightgreen'
    elif float(speed) > float(12) and float(speed) <= float(18):
        color = 'yellow'
    else:
        color = '#f04348'
    return color

def wavecolor(speed):
    color = 'white'
    if speed < 4:
        color = 'lightgreen'
    else:
        color = 'lightred'
    return color

def lengthcolor(speed):
    color = 'white'
    if speed < 19:
        color = 'lightgreen'
    else:
        color = 'lightred'
    return color

def getDataAdaptor():
    ##DATA ADAPTOR 1##
    df = pd.read_csv('dataset.csv')
    df['DateTime'] = pd.to_datetime(df['DateTime'])
    time = df['DateTime']
    hs = df['Hs']
    dir = df['Dir']
    meanperiod = df['meanperiod']
    windx = df['Windx']
    windy = df['Windy']
    peakdir = df['peakdir']
    peakperiod = df['peakperiod']
    ##DATA ADAPTOR 1##

    ##DATA ADAPTOR 2##
    df_vis = pd.read_csv('gfs_3hourly_interpolated.csv')
    df_vis['datetime'] = pd.to_datetime(df_vis['datetime'])
    df_vis['wind_gust'] = df_vis['wind_gust'].apply(lambda x: x*1.9438452)
    #df_vis['wind_gust'] = df_vis['wind_gust'].apply(lambda x: x*(10/0.1)**(1/7))
    vis_df = df_vis['vis']
    sst_df = df_vis['sst']
    gust_df = df_vis['wind_gust']
    vis_tt = df_vis['datetime']
    ##DATA ADAPTOR 2##

    ##DATA ADAPTOR 3##
    df_tide = pd.read_csv('tide_data.csv')
    df_tide['datetime'] = pd.to_datetime(df_tide['datetime'])
    tidal_time = df_tide['datetime']
    tide_height = df_tide['height']
    ##DATA ADAPTOR 3##

    ##CALCULATIONS
    WDIR= (270-np.rad2deg(np.arctan2(windy,windx)))%360
    WSPD = np.sqrt(np.square(windy)+np.square(windx))
    #gust_df.multiply(0.54)
    #gust_df.multiply(1.9438452)
    dfwind = pd.DataFrame({'WSPD':WSPD})
    dfwind['WSPD'] = dfwind['WSPD'].apply(lambda x: x*1.9438452)
    #WSPD.multiply(1.9438452)
    WSPDx = dfwind['WSPD']
    #new_list = [i * 10 for i in WSPD]
    #print(type(WSPD))
    dfgustcheck = pd.DataFrame({'datetime':vis_tt, 'gust':gust_df, 'windspeed':WSPDx})

    for index, row in dfgustcheck.iterrows():
        gst = row['gust']
        spd = row['windspeed']
        if gst < spd:
            dfgustcheck.at[index,'gust'] = spd
    gust_df_corrected = dfgustcheck['gust']

    return time, hs, dir, meanperiod, peakdir, peakperiod, vis_df, sst_df, gust_df_corrected, vis_tt, tidal_time, tide_height, WDIR, WSPDx.values
    

figsize = cm2inch((21,29.7))
dateLabel = datetime.datetime.now(timezone.utc).strftime("%d/%m/%Y")

##CALCULAIONS##
def first_page(folder_tmp,location_name,target_lat,target_lon,tide_name):
    time, hs, dir, meanperiod, peakdir, peakperiod, vis_df, sst_df, gust_df, vis_tt, tidal_time, tide_height, WDIR, WSPD = getDataAdaptor()
    
    degree_sign = u'\N{DEGREE SIGN}'
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    ax.text(-0.05, -0.060,"Date Printed:"+dateLabel, transform=ax.transAxes,fontsize=8, verticalalignment='top')
    ax.text(-0.05, -0.080,"All date and time are in UTC. For more information contact: tuvmet@gmail.com", transform=ax.transAxes,fontsize=6, verticalalignment='top',style='italic')
    ax.text(0.76, -0.060,"Copyright 2023 Tuvalu Met Service", transform=ax.transAxes,fontsize=8, verticalalignment='top')
    ax.text(0.01, 1.12,"Tuvalu Ocean Forecast", transform=ax.transAxes,fontsize=16, verticalalignment='top', color='#2F5F8A')
    ax.text(0.01, 1.09,"Tuvalu Meteorological Service", transform=ax.transAxes,fontsize=8, verticalalignment='top', color='#2F5F8A')
    ax.text(0.01, 1.07,"Contact: tuvmet@gmail.com", transform=ax.transAxes,fontsize=6, verticalalignment='top', color='#2F5F8A')
    ax.axis('off')
    axline = fig.add_axes([0.001, 0.90, 0.99, 0.02])
    axline.axhline(y=0.5, color='#2F5F8A', linestyle='-',linewidth = 15)
    axline.axis('off')
    issue_time = datetime.datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")
    axline.text(0.05, 1.043,"Forecast Report for "+str(round(target_lat,1))+"° "+str(round(target_lon,1))+"° "+location_name+", issued on: "+issue_time+" UTC", fontweight='bold', transform=ax.transAxes,fontsize=8, verticalalignment='top', color='white',)

    ###PLOT 1 WIND SPEED AND DIRECTION ###
    ax2 = fig.add_axes([0.06, 0.76, 0.90, 0.10])
    i = 0
    speed = 1
    for x in time:
        deg = WDIR[i]
        if deg <=0:
            deg+360
        if i == 0:
            i+=1
            continue
        if deg > 270:
            angle = math.radians(270-deg)
        else:
            angle = math.radians(abs(deg-270))   # Remember to convert to radians!

        change = [speed * math.cos(angle), speed * math.sin(angle)]
        if is_multiple_of_three(i):
            ax2.quiver(x, WSPD[i], change[0], change[1],zorder=1, color=windcolor(WSPD[i]),scale=25,width=0.01, edgecolor='black', linewidth=0.7,headwidth=2.5, headlength=2.5, headaxislength=2.4)
           # ax2.quiver(x, windspeedx[i], -1, 1,color=windcolor(windspeedx[i]),width=0.01, edgecolor='black', linewidth=0.7,headwidth=2.5, headlength=2.5, headaxislength=2.4)

        i+=1
    ax2.plot(time.values, WSPD,label='Wind Speed',zorder=3,color='black')
    ax2.plot(time.values, gust_df.values,label='Gusts',linestyle='dashed', zorder=2,color='#dcdad8')
    ax2.legend(prop={'size': 6})
    setupAxis(ax2)
    ax2.set_ylim([WSPD.min()-10,WSPD.max()+20])
    ax2.set_ylabel("Wind Speed [kts]",fontsize=7)
    ax2.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax2.tick_params(
    axis='x',          # changes apply to the x-axis
    which='minor',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
    ax2.margins(0)
    ###PLOT 1 WIND SPEED AND DIRECTION ###

    ###PLOT 2 SIGNIFICANT WAVE HEIGHT ###
    ax11 = fig.add_axes([0.06, 0.65, 0.90, 0.10],sharex=ax2)
    ax11.plot(time.values, hs.values, color='black', label='Sig Wave Height', linewidth=1.5)
    ax11.set_ylabel("Wave Height [m]",fontsize=7)
    ax11.fill_between(time.values, hs.values,color='#dcdad8')
    ax11.set_ylim([0,hs.max()+1])
    setupAxis(ax11)
    ax11.legend(prop={'size': 6})
    ax11.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
    ax11.margins(0)
    ###PLOT 2 SIGNIFICANT WAVE HEIGHT ###

    ###PLOT 3 PEAK WAVE PERIOD ###
    ax33 = fig.add_axes([0.06, 0.54, 0.90, 0.10],sharex=ax2)
    ax33.grid(True)
    ax33.plot(time.values,peakperiod.values, color='black', label='Peak Period',zorder=2)

    i = 0
    for x in time:
        if i == 0:
            i+=1
            continue
        deg = peakdir[i]
        if deg > 270:
            angle = math.radians(270-deg)
        else:
            angle = math.radians(abs(deg-270))   # Remember to convert to radians!

        change = [speed * math.cos(angle), speed * math.sin(angle)]
        if is_multiple_of_three(i):
            ax33.quiver(x, peakperiod.values[i], change[0], change[1], zorder=1, color='#dcdad8',scale=25,width=0.01, edgecolor='black', linewidth=0.7,headwidth=2.5, headlength=2.5, headaxislength=2.4)
        
        i+=1

    ax33.set_ylim([peakperiod.values.min()-5,peakperiod.values.max()+10])
    ax33.set_ylabel("Peak Period [s]",fontsize=7)
    setupAxis(ax33)
    ax33.legend(prop={'size': 6})
    ax33.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
    ax33.margins(0)
    ###PLOT 3 PEAK WAVE PERIOD ###

    ###PLOT 3 MEAN WAVE PERIOD ###
    ax66 = fig.add_axes([0.06, 0.43, 0.90, 0.10],sharex=ax2)
    ax66.grid(True)
    ax66.plot(time.values,meanperiod.values, color='black', label='Mean Period',zorder=2)

    i = 0
    for x in time:
        if i == 0:
            i+=1
            continue
        deg = dir[i]
        if deg > 270:
            angle = math.radians(270-deg)
        else:
            angle = math.radians(abs(deg-270))   # Remember to convert to radians!

        change = [speed * math.cos(angle), speed * math.sin(angle)]
        if is_multiple_of_three(i):
            ax66.quiver(x, meanperiod.values[i], change[0], change[1],zorder=1, color='#dcdad8',scale=25,width=0.01, edgecolor='black', linewidth=0.7,headwidth=2.5, headlength=2.5, headaxislength=2.4)
        
        i+=1

    ax66.set_ylim([meanperiod.values.min()-5,meanperiod.values.max()+10])
    ax66.set_ylabel("Mean Period [s]",fontsize=7)
    setupAxis(ax66)
    ax66.legend(prop={'size': 6})
    ax66.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
    ax66.margins(0)
    ###PLOT 3 MEan WAVE PERIOD ###

    ###PLOT 4 TIDES ###
    start_date = get_start_date()
    time_tide_min, tide_min, htindex, ltindex = generate_tides(tide_name,start_date,folder_tmp)

    ax4 = fig.add_axes([0.06, 0.32, 0.90, 0.10],sharex=ax2)
    setupAxis(ax4)
    ax4.fill_between(tidal_time, tide_height, color='#46769B')
    ax4.set_ylabel("Tide [m]",fontsize=7)
    ax4.set_ylim([tide_height.min()-0.5,tide_height.max()+0.5])

    for h in range(4,len(htindex)):
        text=time_tide_min[htindex[h]].strftime("%H:%M")
        ax4.text(time_tide_min[htindex[h]]-timedelta(hours=3),tide_min[htindex[h]]+0.08,text,color= 'gray',size=6) 
    
    for h in range(4,len(ltindex)):
        text=time_tide_min[ltindex[h]].strftime("%H:%M")
        ax4.text(time_tide_min[ltindex[h]]-timedelta(hours=3),tide_min[ltindex[h]]-0.20,text,color= 'gray',size=6) 
    
    ax4.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
    ax4.margins(0)
    
    ###PLOT 4 TIDES ###

    ###PLOT 5 VISIBILITY #########
    ##DUMMY
    ax55 = fig.add_axes([0.06001, 0.25, 0.90, 0.05],sharex=ax4)
    ax55.plot(vis_tt.values, vis_df.values,color='white')
    ax55.yaxis.set_ticks([])
    setupAxis(ax55)
    #DUMMY
    ax5 = fig.add_axes([0.06, 0.25, 0.90, 0.05])
    ax5.margins(x=0)
    bining_values = [1000,4000,10000,20000,1000000]
    bins = np.array([1000,4000,10000,20000,1000000])
    inds_ = np.digitize(vis_df.values, bins)
    inds = np.stack((inds_, inds_), axis=1)
    y_ =np.array([0,1])
    y, x = np.meshgrid(y_,vis_tt.values)
    cmap = mpl.colors.ListedColormap(['#BB3F3F', '#DB804E','#EEDC5B', '#49759C', '#5B7C4E'])
    ##NEW
    uniquevalues = np.unique(inds_)
    sortedvalues = np.sort(uniquevalues)
    bin_arr = [0,1,2,3,4]

    for i in sortedvalues:
        if i in bin_arr:
            bin_arr.remove(i)
            
    valuesToRemove = []
    for x in bin_arr:
        val = getvalue(x)
        valuesToRemove.append(val)

    for i in valuesToRemove:
        if i in bining_values:
            bining_values.remove(i)


    colorarr = []
    for x in bining_values:
        colorarr.append(getColor(x))
    ##END NEW
    bins = np.array(bining_values)
    inds_ = np.digitize(vis_df.values, bins)
    inds = np.stack((inds_, inds_), axis=1)
    y_ =np.array([0,1])
    y, x = np.meshgrid(y_,vis_tt.values)
    cmap = mpl.colors.ListedColormap(colorarr)

    ax5.pcolor(x,y, inds, cmap= cmap,linewidth=0,rasterized=True)
    ax5.set_ylabel("Vis",fontsize=7)
    ax5.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
    ax5.tick_params(axis='both', which='both', labelsize=6)
    ax5.yaxis.set_ticks([])
    ###PLOT 5 VISIBILITY #########
    
    ##ADDING LEGEND
    rect = mpl.patches.Rectangle(
    (0, 1.07), width=0.05, height=0.08, color="#BB3F3F", transform=ax55.transAxes,
    clip_on=False
    )
    rect2 = mpl.patches.Rectangle(
        (0.11, 1.07), width=0.05, height=0.08, color="#DB804E", transform=ax55.transAxes,
        clip_on=False
    )
    rect3 = mpl.patches.Rectangle(
        (0.21, 1.07), width=0.05, height=0.08, color="#EEDC5B", transform=ax55.transAxes,
        clip_on=False
    )
    rect4 = mpl.patches.Rectangle(
    (0.35, 1.07), width=0.05, height=0.08, color="#49759C", transform=ax55.transAxes,
    clip_on=False
    )
    rect5 = mpl.patches.Rectangle(
        (0.48, 1.07), width=0.05, height=0.08, color="#5B7C4E", transform=ax55.transAxes,
        clip_on=False
    )
    
    ax55.add_patch(rect)
    ax55.text(0.06, 1.17, 'Fog',fontsize=6, color='black',transform=ax55.transAxes, verticalalignment='top')
    ax55.text(0.17, 1.17, 'Poor',fontsize=6,color='black',transform=ax55.transAxes, verticalalignment='top')
    ax55.text(0.28, 1.17, 'Moderate',fontsize=6,color='black',transform=ax55.transAxes, verticalalignment='top')
    ax55.text(0.41, 1.17, 'Good',fontsize=6,color='black',transform=ax55.transAxes, verticalalignment='top')
    ax55.text(0.54, 1.17, 'Unlimited',fontsize=6,color='black',transform=ax55.transAxes, verticalalignment='top')
    ax55.add_patch(rect2)
    ax55.add_patch(rect3)
    ax55.add_patch(rect4)
    ax55.add_patch(rect5)
    
    ###PLOT 6 VISIBILITY ###
    data = [['Hs', 'Significant Wave Height (m)'],
        ['Tp', 'Peak Period (s)'],
        ['PkDir','Peak Direction (degrees)'],
        ['Wsp', 'Wind Speed (kts)'],
        ['Tm', 'Mean Period (s)']]

    ax6 = fig.add_axes([0.33, 0.10, 0.90, 0.05])
    ax6.text(-0.25, 1.19, 'Abbreviation:',fontsize=7,fontweight='bold',color='black',transform=ax6.transAxes, verticalalignment='top')
    table = ax6.table( cellText=data, loc='left', cellLoc='left', rowLoc='left', colLoc='left')
    table.auto_set_column_width(col=list(range(len(data))))

    # Format the table
    table.auto_set_font_size(False)
    table.set_fontsize(6)
    table.scale(1, 0.7)
    ax6.axis('off')

    data2 = [
        ['WD','Wind Direction (degrees)'],
        ['Vis', 'Visibility (km)'],
        ['SST','Sea Surface Temperature ('+degree_sign+'Celsius)'],
        ['Gst','Typical Gust Speed (kts)'],
        ['Mwd','Mean Wave Direction (degrees)']]

    ax7 = fig.add_axes([0.568, 0.10, 0.90, 0.05])
    table = ax7.table( cellText=data2, loc='left', cellLoc='left', rowLoc='left', colLoc='left')
    table.auto_set_column_width(col=list(range(len(data))))

    # Format the table
    table.auto_set_font_size(False)
    table.set_fontsize(6)
    table.scale(1, 0.7)
    ax7.axis('off')
    
    #LOGO
    im = plt.imread('Flag-Tuvalu.jpg')
    newax = fig.add_axes([0.83, 0.92, 0.15, 0.06], anchor='NE', zorder=-1)
    newax.imshow(im)
    newax.axis('off')

    #LOGO
    im = plt.imread('tuvalumet.jpg')
    newax2 = fig.add_axes([0.58, 0.92, 0.17, 0.065], anchor='NE', zorder=-1)
    newax2.imshow(im)
    newax2.axis('off')

    im = plt.imread('Coat.png')
    newax3 = fig.add_axes([0.67, 0.93, 0.15, 0.049], anchor='NE', zorder=-1)
    newax3.imshow(im)
    newax3.axis('off')
    return fig

def create_table():
    time, hs, dir, meanperiod, peakdir, peakperiod, vis_df, sst_df, gust_df, vis_tt, tidal_time, tide_height, WDIR, WSPD = getDataAdaptor()

    winDirectionDegrees = []
    for x in WDIR:
        if x <= 0:
            x+360
        winDirectionDegrees.append(x)
    
    df = pd.DataFrame({'DateTime':time,'Hs':hs, 'Tp':peakperiod,'PkDir':peakdir,'Tm':meanperiod,'Mwd':dir,'Wsp': WSPD,'Gst':gust_df, 'WD':winDirectionDegrees, 'Vis':vis_df, 'SST':sst_df})
    decimals = 1  
    df['Hs'] = df['Hs'].apply(lambda x: round(x, decimals))
    df['Tp'] = df['Tp'].apply(lambda x: round(x, decimals))
    df['PkDir'] = df['PkDir'].apply(lambda x: round(x, 0))
	
    df['Tm'] = df['Tm'].apply(lambda x: round(x, decimals))
    df['Mwd'] = df['Mwd'].apply(lambda x: round(x, 0))

    df['Wsp'] = df['Wsp'].apply(lambda x: round(x, decimals))
    df['WD'] = df['WD'].apply(lambda x: round(x, 0))
    df['Gst'] = df['Gst'].apply(lambda x: round(x, decimals))


    df['Vis'] = df['Vis'].apply(lambda x: round(x, 0))
    df['SST'] = df['SST'].apply(lambda x: round(x, decimals))

    df['Vis'] = df['Vis'].astype('int')
    df['PkDir'] = df['PkDir'].astype('int')
    df['Mwd'] = df['Mwd'].astype('int')
    df['WD'] = df['WD'].astype('int')

    
    return df

def mslp_figure():
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    ax.text(-0.05, 1.09,"Regional MSLP and surface winds Map", transform=ax.transAxes,fontsize=12, verticalalignment='top', color='#2F5F8A')
    ax.axis('off')

    ###PLOT 1 WIND SPEED AND DIRECTION ###
    im = plt.imread('msl.png')
    ax2 = fig.add_axes([0.03, 0.03, 0.95, 0.90])
    ax2.imshow(im,aspect='auto')
    ax2.axis('off')

    ax2.text(-0.05, -0.060,"Date Printed:"+dateLabel, transform=ax.transAxes,fontsize=8, verticalalignment='top')
    ax2.text(-0.05, -0.080,"All date and time are in UTC. For more information contact: tuvmet@gmail.com", transform=ax.transAxes,fontsize=6, verticalalignment='top',style='italic')
    ax2.text(0.76, -0.060,"Copyright 2023 Tuvalu Met Service", transform=ax.transAxes,fontsize=8, verticalalignment='top')
    return fig

def produce_report(folder_tmp, report_name,location_name,target_lat,target_lon,tide_name):
    with PdfPages(report_name) as pdf:
        #plotting First Page
        
        fig = first_page(folder_tmp,location_name,target_lat,target_lon,tide_name)
        pdf.savefig()
        df = create_table()
        df_list = np.array_split(df, 3)
        for table in df_list:
            
            colors = []
            for _, row in table.iterrows():
                colors_in_column = ["white", "white", "white","white","white","white","white","white","white","white","white"]
                #colors_in_column[1] =wavecolor(row["Hs"])
                #colors_in_column[2] =lengthcolor(row["Tp"])
                colors_in_column[1] =get_color_from_jet(row["Hs"],0,3)
                colors_in_column[2] =get_color_from_jet(row["Tp"],5,20)
                colors_in_column[4] =get_color_from_jet(row["Tm"],5,15)
                colors_in_column[6] =windcolor(row["Wsp"])  
                colors.append(colors_in_column)
            
            fig, ax = plt.subplots(figsize=figsize, dpi=150)
            tab2 = ax.table(cellText=table.values, colLabels=table.columns, loc='center',
                    colColours=["#c0c0c0"] * 11,cellColours=colors)
            tab2.auto_set_column_width(col=list(range(len(df.columns))))
            ax.axis('off')
            ax.text(-0.05, -0.090,"Date Printed:"+dateLabel, transform=ax.transAxes,fontsize=8, verticalalignment='top')
            ax.text(-0.05, -0.11,"All date and time are in UTC. For more information contact: tuvmet@gmail.com", transform=ax.transAxes,fontsize=6, verticalalignment='top',style='italic')
            ax.text(0.76, -0.090,"Copyright 2023 Tuvalu Met Service", transform=ax.transAxes,fontsize=8, verticalalignment='top')
            pdf.savefig()
        
        fig=mslp_figure()
        pdf.savefig()
        plt.close('all')
        plt.figure().clear()
    return None
