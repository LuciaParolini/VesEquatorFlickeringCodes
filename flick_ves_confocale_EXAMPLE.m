
% Example of code to run the flickering analysis
% This is for now mainly for confocal and epifluorescence images. For
% confocal images works very well. Maybe for epifluorescence it needs to be
% double checked. For movies taken with confocal the pixel dimension is
% extracted from the metadata.
% Could work also for the .movies - needs to be checked tough - (SERIES)

clear all
close all

addpath('E:\MATLABscripts\VesEquatorFlickeringCodes') % put here the right path
addpath('E:\MATLABscripts\ReaderFunctions')           % put here the right path
addpath('E:\MATLABscripts\bfmatlab5.1.10')            % put here the right path

% command to have the figures docked into the same window.
set(0,'DefaultFigureWindowStyle','docked')
warning('off','MATLAB:namelengthmaxexceeded')

DIR='./';

ZeissFiles = dir('*.czi');
nameFILEtt = {ZeissFiles(:).name};

% few things to define for the analysis (put this into the for if the parameters change
% for the different files)

zoom = 1;
NSERIE = 1;
Nrad=360;                   % number of points used to map the contour
dpix = [];                  % set the value if the reader doesn't extract it from the metadata (in m)
Temp=273.15+25;             % temperature of the experiment, in K
WIDTH=30;                   % width (in pixels) of the detection ring
plotting=true;
fluo=false;                 % false confocal, true epifluorescence

for FILES = 1:size(nameFILEtt,2)
    
    ves = char(nameFILEtt{FILES});
    
    % INTRODUCE FOR ON SERIES IF YOU HAVE MORE THAN ONE
    
    % if nothing is specified it goes from the first to the last frame,
    % otherwise specify it in the argument as a vector with first frame as
    % first element and last frame as second [first last], put it after
    % zoom
    
    try
        movie=flickering_fluo_mod([DIR,ves],Nrad,Temp,NSERIE,zoom);
        %movie=flickering_Davide_3([DIR,ves],Nrad,dpix,Temp,NSERIE,zoom);
        %(use this for brightfield)
    catch
        disp('There is no file')
        break
    end
    
    if isempty(movie.dpix) == 1
        movie.dpix = dpix;
    end
    
    % initialization of the radius and center and contour finding
    
    movie.set_center('auto')
    movie.set_radius(WIDTH,'auto')
    movie.analyse_contour(true,fluo)
    
    % select frames without spikes and compute the spectrum
    
    close all
    ttFrames = [movie.Frames2Analyse(1):size(movie.contour_fine,1)];
    BADframes = checkCONTOUR(movie);
    GoodFrames = setdiff(ttFrames,BADframes);
    movie.get_spectrum(true,GoodFrames);
    
    % fit 
    
    movie.fit_vesicles_fluctuations(true);
    
    % save (change name if you want to)
    
    save([ves(1:end-4),'_flickering_WIDTH',num2str(WIDTH),'.mat'])
    
    % to print fig in pdf
    %{
    set(gcf,'paperunits','centimeters');
    set(gcf,'paperposition', [1 1 11 7]); % left bottom width height
    set(gca,'ticklength',[0.02 0.01],'Linewidth',0.7)
    namefig_eps = ['ves_' num2str(FILES) '.eps'];
    namefig_pdf = ['ves_' num2str(FILES) '.pdf'];
    print(namefig_eps, '-depsc');
    system(sprintf(['ps2pdf -dEPSCrop ' namefig_eps ' ' namefig_pdf]));
    %}
    
end