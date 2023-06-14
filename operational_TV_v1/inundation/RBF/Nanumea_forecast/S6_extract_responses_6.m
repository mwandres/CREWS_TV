clear all
close all
%%
addpath('F:\XBeach\Tuvalu\1D_sims\codes')

path_transect = 'F:\XBeach\Tuvalu\1D_sims\Nanumea_forecast\Nanumea\Nanumea_Merged_Profilev3.csv';
path_eva = 'F:\XBeach\Tuvalu\1D_sims\Nanumea\';
path_tcs = [path_eva 'TC_Nanumea_synth.mat'];
base_path = 'F:\XBeach\Tuvalu\1D_sims\Nanumea_forecast';

%%
transects_csv = path_transect;
transects = readtable(transects_csv);
transect_labels = unique(transects.LINE_ID);
% load data from transects
for i = 1:length(transect_labels)
    t{i}.label = transect_labels(i);
    ix = find(transects.LINE_ID == transect_labels(i));
    t{i}.x = transects.X(ix);
    t{i}.y = transects.Y(ix);
    [t{i}.lat,t{i}.lon] = utm2ll(t{i}.x,t{i}.y,-60);
    t{i}.z = transects.Z(ix);
    t{i}.angle = transects.Angle(ix);
    t{i}.dist = transects.DIST(ix);
    t{i}.dist_surf = transects.DIST_SURF(ix);
end
number_of_transects = length(t);
figure()
hold on
for i = 1:length(transect_labels)
plot(t{i}.lon,t{i}.lat)
text(t{i}.lon(end),t{i}.lat(end),num2str(transect_labels(i)))
end
fileprint=['Transects.png'];
print('-dpng','-r200',fileprint)
close
%%

%%
fff = 0;
for i = 268:319%length(transect_labels)
    dirName1 = [base_path '/t_' num2str(t{i}.label) '_base'];
    load([dirName1 '/centroid.mat']);
    dist = t{i}.dist';
    dist = -dist';
    Sx.max_twl_nearshore = nan(length(subset),1);
    Sx.inundation_extent = nan(length(subset),1);
    for j = 1:length(subset)
        dirName2 = [dirName1 '/Run_' num2str(j)];
        ncpath = [dirName2 '/xboutput.nc'];
        disp(ncpath)
        try
            zs_max = ncread(ncpath,'zs_max');
            disp('zs_max loaded')
            if mean(zs_max) > 100
                breaker = ncread(ncpath,'breaker'); %use this to break the try statement
            end
            zb = squeeze(ncread(ncpath,'zb'));
            zb = zb(:,1);
            disp('zb loaded')
            x_xb = ncread(ncpath,'globalx');
            disp('x_xb loaded')
            u = ncread(ncpath,'u');
            utest = squeeze(u(:,1,300));
            disp('u loaded')
            zs = ncread(ncpath,'zs');
            disp('zs loaded')
            for k = 1:length(u(:,1,1))
                umax(k,1) = max(u(k,1,:));
            end
            aaa = find(zb<0.1);
            zb(aaa) = nan;
            inundation(:,1) = zs_max - zb;
            bbb = find(inundation <= 0.05);
            inundation(bbb) = nan;
            ccc = find(umax <= 0.05);
            inundation(ccc) = nan;
            %     ddd = find(vmax <= 0.0001);
            %     inundation(ddd) = nan;
            %     inundation(:,:) = inundation(:,:);
            aaa = find(~isnan(inundation));
        
            inundationX(:,1) = interp1(x_xb,inundation,dist);
            zs_interp = interp1(x_xb,squeeze(zs),dist);
            
            ix1  = find(~isnan(inundationX),1,'first');
            ix2  = find(~isnan(inundationX),1,'last');
            if isempty(ix1)
                Sx.max_twl_nearshore(j) = 0;
                Sx.inundation_extent(j) = 0;
                ix2=ix3;
            else
                %S.max_twl_nearshore(j) = inundationX(ix2);
                Sx.max_twl_nearshore(j) = prctile(zs_interp(ix2+10,:),98);
                Sx.inundation_extent(j) = dist(ix1) - dist(ix2);
                ix3=ix2;
            end
            %Sx.inundationX(isnan(inundationX)) = 0;
            %Sx.inundation_transect(j,:) = inundationX;
            disp('Clearing variables')
            clear aaa bbb ccc inundationX inundation umax zs x_xb u zb zs_interp zs_max;
            disp('Cleared')
        catch
            Sx.max_twl_nearshore(j) = nan;
            Sx.inundation_extent(j) = nan;
            clear aaa bbb ccc inundationX inundation umax zs x_xb u zb zs_interp zs_max;
            fff = fff+1;
            disp(['Problem in ' ncpath])
            params_path = [dirName2 '/params.txt'];
            disp(['Changing CFL in ' params_path])
            %S = read_write_entire_textfile(params_path);
            fid  = fopen(params_path,'r');
            f=fread(fid,'*char')';
            fclose(fid);
            textobnd = 'CFL          = 0.100000';
            f =  regexprep(f,'CFL          = 0.700000',textobnd);
            f =  regexprep(f,'CFL          = 0.600000',textobnd);
            f =  regexprep(f,'CFL          = 0.500000',textobnd);
            f =  regexprep(f,'CFL          = 0.400000',textobnd);
            f =  regexprep(f,'CFL          = 0.300000',textobnd);
            f =  regexprep(f,'CFL          = 0.200000',textobnd);
            %f =  regexprep(f,'dx           = 0.500000',textobnd2);
            fid  = fopen(params_path,'w');
            fprintf(fid,'%s',f);
            fclose(fid);
            
            out_failed = 'Failed_runs_6/';
            out_path_failed = [out_failed 't_' num2str(t{i}.label) '_base/Run_' num2str(j)  '/'];
            if exist(out_path_failed,'dir')
            else
                mkdir(out_path_failed);
            end
            filenamelist=dir(dirName2);
            for ccc=3:length(filenamelist)
                copyfile([dirName2 '/' filenamelist(ccc).name],out_path_failed)
            end
            binDir = ['D:\XBeach\Tuvalu\1D_sims\xbeach_5848_x64_netcdf\'];
            % Make batchfile to call XBeach executable
            fid = fopen([out_path_failed 'run_xbeach.bat'],'w');
            fprintf(fid,'%s',['call "' binDir  'xbeach.exe"']);
            fclose(fid);
            
            run_path = '/scale_wlg_nobackup/filesets/nobackup/spc03223/TCAP/Forecast/Nanumea_forecast/Failed_runs_6/';
            fid  = fopen(['Failed_runs_6/failed_job_' num2str(fff) '.sl'],'w');
            fprintf(fid,'%s',['#!/bin/bash -e']);
            fprintf(fid,'\n');
            fprintf(fid,'%s',['#SBATCH --job-name=XBeach_Nanumea # job name (shows up in the queue)']);
            fprintf(fid,'\n');
            fprintf(fid,'%s',['#SBATCH --time=60:00:00 # Walltime (HH:MM:SS)']);
            fprintf(fid,'\n');
            fprintf(fid,'%s',['#SBATCH --mem-per-cpu=2048 # number of tasks (e.g. MPI)']);
            fprintf(fid,'\n');
            fl_path = [run_path  't_' num2str(t{i}.label) '_base/Run_' num2str(j)  '/'];
            fprintf(fid,'%s',['cd ' fl_path]);
            fprintf(fid,'\n');
            fprintf(fid,'%s','srun /scale_wlg_persistent/filesets/project/spc03223/Killo/XBeach/trunk/executables_serial/bin/xbeach');
            fprintf(fid,'\n');
			fclose all;
        end
    end
    S.transect_angle = t{i}.angle;
    [lat, lon]=utm2ll(t{i}.x,t{i}.y,-60);
    S.transect_lon = lon;
    S.transect_lat = lat;
    S.transect_dist = dist;
    posok = find(~isnan(Sx.max_twl_nearshore));
    
    S.max_twl_nearshore = Sx.max_twl_nearshore(posok);
    S.inundation_extent = Sx.inundation_extent(posok);
    %S.inundation_transect = Sx.inundation_transect(posok,:);
    S.forcing = subset(posok,:);
    
    outname = ['t_' num2str(t{i}.label) '_results.mat'];
    save(outname,'S');
    A = [S.forcing, S.max_twl_nearshore, S.inundation_extent];
    outname = ['t_' num2str(t{i}.label) '_results.csv'];
    csvwrite(outname,A)
    B = [t{i}.lon(end) t{i}.lat(end) S.transect_lon(ix2+10) S.transect_lat(ix2+10)];
    outname_2 = ['t_' num2str(t{i}.label) '_coordinates.csv'];
    dlmwrite(outname_2,B, 'delimiter', ',', 'precision', 9);
    clear S A Sx;
    disp(i)
end
fid  = fopen('Failed_runs_6/run_all_models.sh','w');
for j = 1:fff
    fprintf(fid,'%s',['sbatch failed_job_' num2str(j) '.sl']);
    fprintf(fid,'\n');
end
fclose(fid);
%%




