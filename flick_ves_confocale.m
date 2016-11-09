addpath('E:\MATLABscripts\matlab_leica')
addpath('E:\MATLABscripts')

clear all
close all


% command to have the figures docked into the same window.
set(0,'DefaultFigureWindowStyle','docked')
warning('off','MATLAB:namelengthmaxexceeded')

DIR='./';
% ves='v3.lif';
% NSERIE=2;

nameFILEtt = {'DOPC_305mMglucose_100uMGR24.lif'};

for FILES = 1
    
    ves = char(nameFILEtt{FILES});
    
    for NSERIE = 1:40
        
        zoom=1;
        % CORE OF THE CODE
        Nrad=360;          % number of points used to map the contour
        dpix=0;             % the real value is set from metadata of files when "flickering" is called
        Temp=273.15+25;     % temperature, in K
        WIDTH=20;           % width (in pixels) of the detection ring
        plotting=true;
        fluo=false;
        try
            movie=flickering_fluo_mod([DIR,ves],Nrad,dpix,Temp,NSERIE,zoom);
            %movie=flickering_Davide_3([DIR,ves],Nrad,dpix,Temp,NSERIE,zoom);
        catch
            disp('There is no file')
            break
        end
        movie.video.pinhole*1e6;
        movie.set_center('auto')
        movie.set_radius(WIDTH,'auto')
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
        save([ves(1:end-4),'_flickering_WIDTH',num2str(WIDTH),'_ves_',num2str(NSERIE),'.mat'])
        
        movie.results.k
        
    end
    
end