addpath ../../matlab_leica
addpath ../../nuova_analisi_confocale
addpath ../ 

clear all
close all


% command to have the figures docked into the same window.
set(0,'DefaultFigureWindowStyle','docked')
warning('off','MATLAB:namelengthmaxexceeded')

DIR='./';
% ves='v3.lif';
% NSERIE=2;

ves='10nov14b.lif';
NSERIE=10;
zoom=1;
%% CORE OF THE CODE
Nrad=360;          % number of points used to map the contour
dpix=0;             % the real value is set from metadata of files when "flickering" is called
Temp=273.15+25;     % temperature, in K
WIDTH=20;           % width (in pixels) of the detection ring
plotting=true;
fluo=false;

movie=flickering_fluo([DIR,ves],Nrad,dpix,Temp,NSERIE,zoom);
% movie=flickering([DIR,ves],Nrad,dpix,Temp,NSERIE,zoom);

movie.video.pinhole*1e6
movie.set_center('manual')
movie.set_radius(WIDTH,'manual')
movie.analyse_contour(true,fluo,500)

%% saving the results

cc=movie.contour_fine(~isnan(movie.contour_fine(:,1)),:);
cc=cc(2:end,:);cc=bsxfun(@rdivide,cc,mean(cc,2)).*mean(cc(:,2));

fftlog=log10(abs(fft(cc')));

figure(11)
imagesc(fftlog)

movie.get_spectrum(true);
%% fit
%%saving the results
pinhole=num2str(movie.video.pinhole*1e6,'%.0f');

movie.fit_vesicles_fluctuations(true);
save([ves(1:end-4),'_pinhole',pinhole,'.mat'],'movie')