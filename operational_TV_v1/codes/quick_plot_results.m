clear all
close all
%%
load('D:\CREWS_TV\operational_TV_v1\runs\2023061218\output.mat')

[fem,bnd]=read_adcirc_mesh(['D:\CREWS_TV\operational_TV_v1\common\fort.14']);

%plot_date = datestr(datenum(yr,mo,dd,HH,0,0),'yyyymmdd_HHMMSS');
        
 Xp=double(Xp);
 Yp=double(Yp);
 
 x = Windv_x_20230616_000000;
 y = Windv_y_20230616_000000;
 
 figure
 h=trisurf(fem.e,Xp,Yp,x);
 view(0,90);shading interp;
 caxis([-3 3])
 colormap('jet')
 c = colorbar;
 ylabel(c,'x wind (m/s)')
 title('x')
 
 
 figure
 h=trisurf(fem.e,Xp,Yp,y);
 view(0,90);shading interp;
 caxis([-6 6])
 colormap('jet')
 c = colorbar;
 ylabel(c,'y wind (m/s)')
 title('y')