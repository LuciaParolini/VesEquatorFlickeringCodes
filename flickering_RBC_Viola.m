
%%%%%% flickering red blood cells %%%%%%%

addpath('E:\MATLABscripts')

clear all
close all

% command to have the figures docked into the same window.
set(0,'DefaultFigureWindowStyle','docked')
warning('off','MATLAB:namelengthmaxexceeded')

name = 'malaria_fast.15May2016_13.02.31.movie';

zoom=1;
% CORE OF THE CODE
Nrad=360;          % number of points used to map the contour
dpix=0.0973*10^-6;        % um 60x e grasshopper
Temp=273.15+37;     % temperature, in K
WIDTH=20;           % width (in pixels) of the detection ring
plotting=true;
fluo=false;

%try
    %movie=flickering_fluo_mod([DIR,ves],Nrad,dpix,Temp,NSERIE,zoom);
    movie=flickering_Davide_3(name,Nrad,dpix,Temp);
%catch
%    disp('There is no file')
%    return
%end
%movie.video.pinhole*1e6;
movie.set_center('manual')
movie.set_radius(WIDTH,'manual')
movie.analyse_contour(true,fluo)

%% saving the results

% cc=movie.contour_fine(~isnan(movie.contour_fine(:,1)),:);
% cc=cc(2:end,:);cc=bsxfun(@rdivide,cc,mean(cc,2)).*mean(cc(:,2));
%
% fftlog=log10(abs(fft(cc')));
%
% figure(11)
% imagesc(fftlog)
%
% movie.get_spectrum(true);

%%
close all
ttFrames = [1:round(size(movie.contour_fine,1))];
BADframes = checkCONTOUR(movie);
GoodFrames = setdiff(ttFrames,BADframes);
movie.get_spectrum(true,GoodFrames);

% fit
%%saving the results
%pinhole=num2str(movie.video.pinhole*1e6,'%.0f');

movie.fit_vesicles_fluctuations(true);

%%
%save([ves(1:end-4),'_flickering_WIDTH',num2str(WIDTH),'_ves_',num2str(NSERIE),'.mat'])

%movie.results.k
