clear all
close all
%%

files = dir('Failed_runs\');
dirFlags = [files.isdir];
subdirs = files(dirFlags);
cnt = 0;
for i = 3:length(subdirs)
    files = dir([subdirs(i).folder '/' subdirs(i).name]);
    dirFlags = [files.isdir];
    subsubdirs = files(dirFlags);
    for j = 3:length(subsubdirs)
        cnt = cnt+1;
        path{cnt} = [subdirs(i).name '/' subsubdirs(j).name];
        
    end
end
%%
run_path = '/nesi/nobackup/spc03223/TCAP/Forecast/Nui_forecast/Failed_runs';
fid  = fopen(['Failed_runs/Run_failed_jobs' '.sl'],'w');
fprintf(fid,'%s',['#!/bin/bash -e']);
fprintf(fid,'\n');
fprintf(fid,'%s',['#SBATCH --job-name=XBeach_Nui # job name (shows up in the queue)']);
fprintf(fid,'\n');
fprintf(fid,'%s',['#SBATCH --time=48:00:00 # Walltime (HH:MM:SS)']);
fprintf(fid,'\n');
fprintf(fid,'%s',['#SBATCH --mem-per-cpu=2048 # number of tasks (e.g. MPI)']);
fprintf(fid,'\n');
fprintf(fid,'%s',['#SBATCH --array=0-' num2str(length(path)-1) '%800']);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'BIN=/nesi/project/spc03223/Killo/XBeach/trunk/executables_serial/bin');
fprintf(fid,'\n');
fprintf(fid,['DATA=' run_path]);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'shopt -s dotglob');
fprintf(fid,'\n');
fprintf(fid,'shopt -s nullglob');
fprintf(fid,'\n');
fprintf(fid,'subdir_paths=($DATA/t*/Run_*)');
fprintf(fid,'\n');
fprintf(fid,['cd ' '"' '${subdir_paths[${SLURM_ARRAY_TASK_ID}]}' '"' ' && $BIN/xbeach']);
fclose all;
