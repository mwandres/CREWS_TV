# -*- coding: utf-8 -*-
"""
Created on Thu May 20 17:51:04 2021

@author: antonioh
"""
import datetime as dt
from matplotlib import path
import os, shutil
from adcircpy import AdcircMesh
import matplotlib.pyplot as plt
# from AdcircPy import read_mesh
import numpy as np
from subprocess import Popen, PIPE

def execute(cmd):
    popen = Popen(cmd, stdout=PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return()

def startenddate(start_str,end_str):
    fmt = '%Y%m%d%H'
    dt1 = dt.datetime.strptime(start_str, fmt)
    dt2 = dt.datetime.strptime(end_str, fmt)
    dt1f = dt1.strftime('%Y%m%d.%H%M%S')
    dt2f = dt2.strftime('%Y%m%d.%H%M%S')
    return(dt1f,dt2f)

def parallelization_fort26(folder_name,t1,t2):
    # directory = os.getcwd()
    # fpath = os.path.join(directory,folder_name)
    flist = os.listdir(folder_name)
    
    pmesh = AdcircMesh.open(folder_name+"/fort.14",crs=4326)
    #pmesh = read_mesh(folder_name+"/fort.14")
    px = pmesh.x
    py = pmesh.y
    pbnd = pmesh.open_boundaries
    pbnd_val = pbnd[0]["indexes"]
    pbnd_x = px[pbnd_val]
    pbnd_y = py[pbnd_val]

    for fl in flist:
        if fl.startswith("PE"):
            subpath = os.path.join(folder_name,fl)
            # sublist = os.listdir(subpath)
            # for fsub in sublist:
            submesh = AdcircMesh.open(subpath+"/fort.14",crs=4326)
            subbnd = submesh.open_boundaries
            
            if  len(subbnd)>0:
                subx = submesh.x
                suby = submesh.y
                subbnd_val = subbnd[0]["indexes"]
                subbnd_val2 = list(set(subbnd_val))        
                subbnd_x = subx[subbnd_val2]
                subbnd_y = suby[subbnd_val2]
                
                pid = np.empty(0, int)
                for sbx,sby in zip(subbnd_x,subbnd_y):
                    pos = np.sqrt((pbnd_x-sbx)**2+(pbnd_y-sby)**2).argmin()
                    # pos = np.sqrt((pbnd_x-subbnd_x[0])**2+(pbnd_y-subbnd_y[0])**2).argmin()
                    pid_val = pbnd_val[pos]
                    pid = np.append(pid,pid_val)
                    
                pidf = pid+1
                subidf = np.array(subbnd_val)+1
                bnd_txt = []
                for pf14id,subid in zip(pidf,subidf):
                    s = ("BOUN SEGMENT IJ %s VAR FILE 0 'Pto_%s.sp2' 1, &\n" % (subid,pf14id))
                    bnd_txt.append(s)
                bnd_txt[-1] = bnd_txt[-1][:-4].strip()
                bnd_str = ''.join(bnd_txt)
        

                # Read in the file
                with open(subpath+"/fort.26", 'r') as f:
                    filedata = f.read()

                # Replace the target string
                filedata = filedata.replace('$%%BOUND_COMMAND%%', bnd_str)
                filedata = filedata.replace('%%BOUND_STARTDATE%%', t1)
                filedata = filedata.replace('%%BOUND_ENDDATE%%', t2)

                # Write the file out again
                with open(subpath+"/fort.26", 'w') as fout:
                    fout.write(filedata)
                    
                   
            else:
                # Read in the file
                with open(subpath+"/fort.26", 'r') as f:
                    filedata = f.read()

                # Replace the target string
                filedata = filedata.replace('%%BOUND_STARTDATE%%', t1)
                filedata = filedata.replace('%%BOUND_ENDDATE%%', t2)

                # Write the file out again
                with open(subpath+"/fort.26", 'w') as fout:
                    fout.write(filedata)
    # return(print("fort26 edited for parallelization"))

def points_subdomains_fort26(folder_name,pointlist,t1):
    flist = os.listdir(folder_name)

    for fl in flist:
        if fl.startswith("PE"):
            subpath = os.path.join(folder_name,fl)
            submesh = AdcircMesh.open(subpath+"/fort.14",crs=4326)           
            kkout=submesh.outer_ring_collection[0][:,0]# outer boundary
            kkinner=submesh.inner_ring_collection[0]

            subx = submesh.x
            suby = submesh.y
            
            ## check whether the points are inside each subdomain
            p = path.Path(np.stack((subx[kkout],suby[kkout]),axis=-1))
            
            # plt.plot(subx,suby,'.')
            # plt.plot(subx[kkout],suby[kkout],'r')
            # # plt.plot(subx[kkin],suby[kkin],'g')
            # plt.plot(pointlist[:,0],pointlist[:,1],'*k')
            
            inside = p.contains_points(pointlist)
                   
            points_txt = []
            table_txt = []
            spect_txt = []
            
       
            for npp in range(len(pointlist[:,1])):
                ninner=np.shape(kkinner)
                ninner=ninner[0]                
                dentro= False
                for ni in range(ninner): 
                #for ni in range(np.size(submesh.inner_ring_collection)): 
                    ## check whether the points are inside a smaller subdomain that is contained in a larger subdamain,
                    ## otherwise there would be several subdomains archieving the same points
                    kkin=submesh.inner_ring_collection[0][ni][:,1]                    
                    p = path.Path(np.stack((subx[kkin],suby[kkin]),axis=-1))
                    inside2 = p.contains_points(pointlist)
                    
                    if inside2[npp]== True:
                        dentro = True
                        
                    
                    
                if inside[npp] == True and dentro == False:
                   
                       po = ("POINTS 'Pto_%s' %s %s \n" % (npp+1,pointlist[npp,0],pointlist[npp,1]))
                       points_txt.append(po)
                       
                       tb=("TABLE  'Pto_%s' HEAD   '../Pto_%s.tab' TIME DEP HS HSWELL RTP PER DIR WIND OUT %s 1 HR \n" % (npp+1,npp+1,t1))
                       table_txt.append(tb)
                       
                       sp=("SPECOUT 'Pto_%s' SPEC2D ABS '../Pto_%s.spec' OUT %s 1 HR \n" % (npp+1,npp+1,t1))
                       spect_txt.append(sp)
                   
                
                
            # Read in the file
            with open(subpath+"/fort.26", 'r') as f:
                filedata = f.read()
            # Replace the target string
            filedata = filedata.replace('$%%POINTS%%', ''.join(points_txt))
            filedata = filedata.replace('$%%TABLE%%', ''.join(table_txt))
            filedata = filedata.replace('$%%SPECT%%', ''.join(spect_txt))

            # Write the file out again
            with open(subpath+"/fort.26", 'w') as fout:
                 fout.write(filedata)

############################################################################################################################################
def par_run(now):
  
    out_name='../runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") + now.strftime("%H")  +'/'
    cygfolder = '/cygdrive/d/CREWS_TV/operational_TV_v1/runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H")  +'/'
    ini = now - dt.timedelta(days=2)
    start_str = ini.strftime("%Y") + ini.strftime("%m") + ini.strftime("%d") + ini.strftime("%H") 
    end = ini + dt.timedelta(days=9.5)
    end_str = end.strftime("%Y") + end.strftime("%m") + end.strftime("%d") + end.strftime("%H") 
    
    
    t1,t2 = startenddate(start_str,end_str)
    
    source_dir ='../common'
    file_names = os.listdir(source_dir)
    for file_name in file_names:
        shutil.copy(os.path.join(source_dir, file_name), out_name)
        
        
    os.system('C:/cygwin64/bin/bash --login -c "cd %s;  /home/executables/adcprep --np 16 --partmesh"' % cygfolder)
    os.system('C:/cygwin64/bin/bash --login -c "cd %s;  /home/executables/adcprep --np 16 --prepall"' % cygfolder)
    parallelization_fort26(out_name,t1,t2)
    
    pointlist=np.array([[179.2139623911262, -8.522669489113921],[179.0812094845414, -8.652164074276669],[179.0000481854969, -8.513060353962492],[179.1357025559911, -8.411012345526808]])
    #pointlist=np.array([[179.2139623911262, -8.522669489113921],[179.0812094845414, -8.652164074276669],[179.01064,-8.4887],[179.1357025559911, -8.411012345526808]])
    points_subdomains_fort26(out_name,pointlist,t1)
    
    os.system('C:/cygwin64/bin/bash --login -c "cd %s; mpirun --use-hwthread-cpus -np 16 /home/executables/punswan"' % cygfolder )
    
    try: 
        os.system('taskkill /F /IM punswan_4131.exe')
    except:
        print('SWAN run finished succesfully')


#now = dt.datetime(2021,12,21)
#par_run(now)