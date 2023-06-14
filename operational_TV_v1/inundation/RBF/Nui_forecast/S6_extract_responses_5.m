clear all
close all
%%
addpath('F:\XBeach\Tuvalu\1D_sims\codes')

path_transect = 'F:\XBeach\Tuvalu\1D_sims\Nui_forecast\Nui\Nui_Oceanside_Profilesv2.csv';
path_eva = 'F:\XBeach\Tuvalu\1D_sims\Nui\';
path_tcs = [path_eva 'TC_Nui_synth.mat'];
base_path = 'F:\XBeach\Tuvalu\1D_sims\Nui_forecast';

%%
transects_csv = path_transect;
transects = readtable(transects_csv);
transects.LINE_ID = str2double(transects.LINE_ID);
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
close all
%%
% rerun_counter = 0;
% for i = 1:length(transect_labels)
%     dirName1 = [base_path '/t_' num2str(t{i}.label) '_base'];
%     load([dirName1 '/centroid.mat']);
%     %dist = t{i}.dist';
%     %dist = -dist';
%     %Sx.max_twl_nearshore = nan(length(subset),1);
%     %Sx.inundation_extent = nan(length(subset),1);
%     for j = 1:length(subset)
%         dirName2 = [dirName1 '/Run_' num2str(j)];
%         ncpath = [dirName2 '/xboutput.nc'];
%         if ~exist(ncpath)
%             disp([ncpath ' does not exist'])
%         else
%             s=dir(ncpath);
%             fl_size = s.bytes;
%             if fl_size < 300000
%                 disp([ncpath ' is very small. Probably crashed.'])
%                 rerun_counter = rerun_counter+1;
%                 params_path = [dirName2 '/params.txt'];
%                 %S = read_write_entire_textfile(params_path);
%                 fid  = fopen(params_path,'r');
%                 f=fread(fid,'*char')';
%                 fclose(fid);
%                 textobnd = 'CFL          = 0.100000';
%                 f =  regexprep(f,'CFL          = 0.700000',textobnd);
%                 f =  regexprep(f,'CFL          = 0.600000',textobnd);
%                 f =  regexprep(f,'CFL          = 0.500000',textobnd);
%                 f =  regexprep(f,'CFL          = 0.400000',textobnd);
%                 f =  regexprep(f,'CFL          = 0.300000',textobnd);
%                 f =  regexprep(f,'CFL          = 0.200000',textobnd);
%                 fid  = fopen(params_path,'w');
%                 fprintf(fid,'%s',f);
%                 fclose(fid);
%                 binDir = ['D:\XBeach\Tuvalu\1D_sims\xbeach_5848_x64_netcdf\'];
%                 % Make batchfile to call XBeach executable
%                 fid = fopen([dirName2 '\run_xbeach.bat'],'w');
%                 fprintf(fid,'%s',['call "' binDir  'xbeach.exe"']);
%                 fclose(fid);
%                 rerun_dir{rerun_counter} = dirName2;
%             end
%         end
%     end
% end
% fid = fopen(['rerun_failed_xbeach.bat'],'w');
% 
% for i = 1:rerun_counter
%     % Make batchfile to call XBeach executable
%     
%     fprintf(fid,'%s\n',['cd ' rerun_dir{i}]);
%     %    fprintf(fid,'%s\n','START /wait CMD /C run_xbeach.bat');
%     fprintf(fid,'%s\n','START CMD /C run_xbeach.bat');
%     % fprintf(fid,'%s\n',['run_xbeach.bat']);
%     
% end
% fclose(fid);
%%
fff = 0;
for i = 5:6:length(transect_labels)
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
                try
                    Sx.max_twl_nearshore(j) = prctile(zs_interp(ix2+10,:),98);
                    Sx.inundation_extent(j) = dist(ix1) - dist(ix2);
                    ix3=ix2;
                catch
                    Sx.max_twl_nearshore(j) = prctile(zs_interp(ix2,:),98);
                    Sx.inundation_extent(j) = dist(ix1) - dist(ix2);
                    ix3=ix2;
                end
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
            
            out_failed = 'Failed_runs_5/';
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
            
            run_path = '/scale_wlg_nobackup/filesets/nobackup/spc03223/TCAP/Forecast/Nui_forecast/Failed_runs_5/';
            fid  = fopen(['Failed_runs_5/failed_job_' num2str(fff) '.sl'],'w');
            fprintf(fid,'%s',['#!/bin/bash -e']);
            fprintf(fid,'\n');
            fprintf(fid,'%s',['#SBATCH --job-name=XBeach_Nui # job name (shows up in the queue)']);
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
fid  = fopen('Failed_runs_5/run_all_models.sh','w');
for j = 1:fff
    fprintf(fid,'%s',['sbatch failed_job_' num2str(j) '.sl']);
    fprintf(fid,'\n');
end
fclose(fid);
%%
load('t_1_results.mat')
figure()
subplot(6,1,1)
plot(S.forcing(:,1))
ylabel('H_s (m)')
subplot(6,1,2)
plot(S.forcing(:,2))
ylabel('T_p (s)')
subplot(6,1,3)
plot(S.forcing(:,3))
ylabel('\theta_p (\circ)')
ylim([0 360])
subplot(6,1,4)
plot(S.forcing(:,4))
ylabel('MSLA (m)')
subplot(6,1,5)
plot(S.max_twl_nearshore)
ylabel('TWL nearshore (m)')
subplot(6,1,6)
plot(S.inundation_extent)
ylabel('Inundation extent (m)')
fileprint=['t_1_results.png'];
print('-dpng','-r200',fileprint)
%%
% ARI_name = {'5-year ARI';'10-year ARI';'25-year ARI';'50-year ARI';'100-year ARI';'250-year ARI'};
% ARI_name_all = [4.9,5,5.1;9.5,10,10.5;24,25,26;49,50,51;99,100,101;240,250,260];
% ARI = [4.9;5;5.1;9.5;10;10.5;24;25;26;49;50;51;99;100;101;240;250;260];
% 
% failed_job_counter = 0;
% for nnn = 1:length(ARI)
%     for i = 1:length(t)
%         dist = t{i}.dist';
%         dist = -dist';
%         %         inundationX = nan(length(dist),3);
%         %         for k = 1:3
%         dirName1 = ['t_' num2str(t{i}.label) '_base'];
%         %dirName2 = [dirName1 '\t_' num2str(t{i}.label) '_ARI_' num2str(ARI_name_all(nnn,k)) '\'];
%         dirName2 = [dirName1 '\t_' num2str(t{i}.label) '_ARI_' num2str(ARI(nnn)) '/'];
%         dirName2 = strrep(dirName2,'\','/');
%         %dirName2 = [dirName1 '\' 't_' num2str(t{i}.label) '_' ARI_name{nnn} '\'];
%         ncpath = [dirName2 'xboutput.nc'];
%         zs_max = ncread(ncpath,'zs_max');
%         zb = squeeze(ncread(ncpath,'zb'));
%         zb = zb(:,1);
%         x_xb = ncread(ncpath,'globalx');
%         u = ncread(ncpath,'u');
%         for j = 1:length(u(:,1,1))
%             umax(j,1) = max(u(j,1,:));
%         end
%         aaa = find(zb<0.1);
%         zb(aaa) = nan;
%         inundation(:,1) = zs_max - zb;
%         bbb = find(inundation <= 0.05);
%         inundation(bbb) = nan;
%         ccc = find(umax <= 0.05);
%         inundation(ccc) = nan;
%         %     ddd = find(vmax <= 0.0001);
%         %     inundation(ddd) = nan;
%         %     inundation(:,:) = inundation(:,:);
%         aaa = find(~isnan(inundation));
%         
%         
%         inundationX(:,1) = interp1(x_xb,inundation,dist);
%         
%         [lat, lon]=utm2ll(t{i}.x,t{i}.y,-60);
%         S.inundation_lon{i}{nnn} = lon;
%         S.inundation_lat{i}{nnn} = lat;
%         S.inundation_result{i}{nnn} = inundationX;
%         ix1  = find(~isnan(inundationX),1,'first');
%         ix2  = find(~isnan(inundationX),1,'last');
%         S.inundation_max_twl_nearshore{i}{nnn} = inundationX(ix2);
%         S.inundation_extent{i}{nnn} = dist(ix1) - dist(ix2);
%         S.inundation_scenario{nnn} = ARI(nnn);
%         S.t_name{i} = dirName1;
%         
%         load([dirName1 '/forcing.mat']);
%         S.forcing{i} = forcing{length(forcing)};
%         
%         clear inundation; clear umax; clear x_xb; clear u; clear zs_max;
%         clear aaa; clear bbb; clear ccc; clear inundationX; clear forcing;
%     end
%     
%     disp(nnn)
% end
% 
% figure(100)
% hold on
% for i = 1:154
%     %plot(inundation_result{1}{1})
%     scatter(S.inundation_lon{i}{1},S.inundation_lat{i}{1},15,S.inundation_result{i}{1},'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
% end
% 
% figure(200)
% hold on
% for i = 1:154
%     plot(S.inundation_result{i}{1})
% end
% 
% save('Nui_results.mat','S')




