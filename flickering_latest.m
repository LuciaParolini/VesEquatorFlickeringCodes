%% FLICKERING CLASS: a Matlab algorithm for vesicles contour detection and analysis

classdef flickering < handle
    % Class properties (variables) that can be set by the user in the
    % Matlab prompt, instead that using the class methods
    properties
        video           % contains the VideoReader object (with video infos)
        frame           % current frame
        current_frame   % ID of the current frame
        center          % center of the vescicle
        mean_radius     % mean value of the radius
        radius_limits   % limits of the ring region that contains the membrane
        sigma           % parameter used to calculate radius_limits from mean_radius
        dt              % Video Framerate
        dpix            % pixel to m conversion, in m
    end
    
    
    % Class properties (variables) that can be set only by the class methods
    properties (SetAccess = private)
        cen
        spectrum        % the spectra obtained from the contour
        results         % results of fitting
        Nrad            % Number of angular rays that maps the contour
        Temp            % temperature, in K
        ring            % contour region
        contour_raw     % logging contours (local parabolic fit) as a 2D matrix
        contour_fine    % logging contours (global fit) as 2D matrix
        xc              % X cartesian coordinate of the contour (for plotting)
        yc              % Y cartesian coordinate of the contour (for plotting)
    end
    
    
    % Class methods (accessible by the users)
    methods
        
        function obj=flickering(filename,Nrad,dpix,Temp)
            % Class constructor. It loads and initializes a
            % video file, creating the initial data structure.
            %
            % USAGE: obj=flickering(filename, Nrad, dpix, Temperature);
            %
            %   filename: the filename of the video to analyze;
            %
            %   Nrad:     number of points used to map the contour
            %             (usually, 360).
            %
            %   dpix:     pixel to meter conversion
            %             (spatial resolution of the camera).
            %
            %   Temperature: the temperature of the sample during the
            %                measurement. Required for correct normalization of the
            %                optput parameters of the flickering analysis
            
            % checking if the matlab pool of workers is active, to disable
            % some warnings.
            if matlabpool('size')>0
                pctRunOnAll warning('off','MATLAB:nearlysingularmatrix')
            else
                warning('off','MATLAB:nearlysingularmatrix')
            end
            
            
            obj.Nrad=Nrad;
            obj.Temp=Temp;
            obj.dpix=dpix;
            
            
            %reading the video file.
            
            if strfind(filename,'.movie')
                if which('moviereader');
                    obj.video=moviereader(filename);
                else
                    disp('Cannot open .movie videos');
                end
            elseif strfind(filename,'.lif')
                if which('lifreader');
                    obj.video=lifreader(filename,0.25);
                    obj.dpix=obj.video.dpix;
                else
                    disp('Cannot open .lif videos');
                end
            else
                obj.video=VideoReader(filename);
            end
            %load first frame
            obj.frame=get_frame(obj.video,1);
            obj.current_frame=1;
            
            obj.dt=1./obj.video.FrameRate;
        end
        
        function obj=set_center(obj,varargin)
            % Function that sets the center of the vescicles.
            %
            % USAGE:
            %   obj.set_center;
            %   obj.set_center(method);
            %   obj.set_center(method,plotting);
            %
            % plotting: boolean value that activates plotting (default=true)
            %
            % method is a string that sets the center detection algorithm:
            %
            %   'manual':     center selected by a mouse click;
            %
            %   'basic_symm': (default) the center is desumed roughly as the center of
            %                 gravity of the grayscale images. The vesicles
            %                 has to be roughly at the center of the image;
            %
            %   'auto':       detected as the center of the largest object
            %                 in the automatically-thresholded image. Works
            %                 only with low-noise frames.
            
            method='basic_symm';
            plotting=true;
            
            if nargin==2;
                method=varargin{1};
                plotting=true;
            end
            
            if nargin==3
                method=varargin{1};
                plotting=varargin{2};
                if strcmp(method,'manual')
                    plotting=true;
                end
            end
            
            if plotting
                show_img(1,obj.frame)
            end
            if strcmp(method,'manual')
                title('Click on the center of the vesicle')
                obj.center=ginput(1);
                
            elseif strcmp(method,'basic_symm')
                xline=mean(obj.frame,1);
                yline=mean(obj.frame,2);
                x=1:length(xline);
                y=(1:length(yline))';
                obj.center=[sum(x.*xline)./sum(xline),sum(y.*yline)./sum(yline)];
                
            elseif strcmp(method,'auto')
                
                [gx,gy]=gradient(obj.frame);
                
                gg=sqrt(gx.^2+gy.^2);
                bckg=imfilter(gg,fspecial('disk',10));
                IM=gg-bckg;
                L=graythresh(IM);
                
                bw=IM>L;
                bw=bwmorph(bw,'open',1);
                bw=bwmorph(bw,'clean',1);
                bw=bwmorph(bw,'close',1);
                
                C=regionprops(bw,'Area','Centroid');
                [~,ind]=max([C.Area]);
                
                obj.center=C(ind).Centroid;
                
            end
            
            if plotting
                hold on
                plot(obj.center(1),obj.center(2),'r.')
            end
        end
        
        function obj=load_frame(obj,ind)
            % USAGE obj.load_frame( iFrame );
            % load the ith frame of the video, for debugging/control purposes.
            obj.current_frame=ind;
            obj.frame=get_frame(obj.video,ind);
        end
        
        function obj=set_radius(obj,delta,varargin)
            % Function that sets the radius of the vescicles.
            % USAGE:
            %   obj.set_radius(width);
            %   obj.set_radius(width,method);
            %   obj.set_radius(width,method,plotting);
            %
            % width: width in pixels of a ring, roughly centered on the vesicle's
            %        contour, that contains the contour itself (usually set
            %        to 70)
            %
            % plotting: boolean value that activates plotting (default=true)
            %
            % method is a string that sets the radius detection algorithm:
            %
            %   'manual':     radius selected by a mouse click;
            %
            %   'basic_symm': (default) the radius is desumed roughly as the distance
            %                 of the maximum of the grayscale image center
            %                 with respect to the center of gravity of the
            %                 grayscale images. The vesicles (contour in white)
            %                 has to be roughly at the center of the image;
            %
            %   'auto':       detected as the radius of the largest object
            %                 in the automatically-thresholded image. Works
            %                 only with low-noise frames.
            
            method='basic_symm';
            plotting=true;
            
            if nargin==3;
                method=varargin{1};
                plotting=true;
            end
            
            if nargin==4
                method=varargin{1};
                plotting=varargin{2};
                if strcmp(method,'manual')
                    plotting=true;
                end
            end
            
            if plotting
                show_img(1,obj.frame)
                hold on
                plot(obj.center(1),obj.center(2),'r.')
            end
            
            if strcmp(method,'manual')
                title('Select a point on the membrane')
                p=ginput(1);
                
                obj.mean_radius=sqrt(sum((p-obj.center).^2));
                
            elseif strcmp(method,'basic_symm')
                a=abs(mean(obj.frame(round(obj.center(1))+(-10:1:10),:),1));
                [~,II]=max(a);
                r1=abs(II-round(obj.center(1)));
                
                b=abs(mean(obj.frame(:,round(obj.center(2))+(-10:1:10)),2));
                [~,II]=max(b);
                r2=abs(II-round(obj.center(2)));
                obj.mean_radius=mean([r1,r2]);
                
            elseif strcmp(method,'auto')
                
                [gx,gy]=gradient(obj.frame);
                
                gg=sqrt(gx.^2+gy.^2);
                bckg=imfilter(gg,fspecial('disk',10));
                IM=gg-bckg;
                L=graythresh(IM);
                
                bw=IM>L;
                bw=bwmorph(bw,'open',1);
                bw=bwmorph(bw,'clean',1);
                bw=bwmorph(bw,'close',1);
                
                C=regionprops(bw,'Area','BoundingBox');
                [~,ind]=max([C.Area]);
                
                obj.mean_radius=mean(C(ind).BoundingBox(3:4))/2;
                
            end
            
            % calculating radius and limits
            
            obj.sigma=delta/2/obj.mean_radius;
            obj.radius_limits=obj.sigma*obj.mean_radius*[1,-1]+obj.mean_radius;
            
            if plotting
                title('')
                ang=(0:0.1:360)';
                r=obj.mean_radius*[cosd(ang),sind(ang)];
                r1=obj.radius_limits(1)*[cosd(ang),sind(ang)];
                r2=obj.radius_limits(2)*[cosd(ang),sind(ang)];
                
                plot(obj.center(1)+r(:,1),obj.center(2)+r(:,2),'b--')
                plot(obj.center(1)+r1(:,1),obj.center(2)+r1(:,2),'g-.')
                plot(obj.center(1)+r2(:,1),obj.center(2)+r2(:,2),'g-.')
            end
        end
        
        function plot_contour(obj,fig)
            % function to plot the contour stored into obj.xc and obj.yc
            % (for debugging purposes)
            %
            % USAGE obj.plot_contour(figure_number)
            
            show_img(fig,obj.frame);
            hold on
            title(['Frame # ',num2str(obj.current_frame)])
            plot(obj.center(1),obj.center(2),'r.')
            plot(obj.xc+obj.center(1),obj.yc+obj.center(2),'r-');
        end
        
        function obj=analyse_contour(obj,varargin)
            % Function that detects the contour of the vesicle in every
            % frame of the movie. The video is split into 60-frames bunches,
            % that are processed in parallel if the "matlabpool numworkers"
            % command has been issued before.
            %
            % USAGE:
            %   obj.analyse_contour;
            %   obj.analyse_contour(plotting);
            %   obj.analyse_contour(plotting,method);
            %
            % plotting: boolean variable to activate plotting (default: true)
            %
            % method:   "profile" (default) the contour is detected by
            %           correlating each radial profile with a previously
            %           determined template.
            
            NF=obj.video.NumberOfFrames;
            if NF>6000
                NF=6000;
                disp(['Analyszing ',num2str(NF),' frames of ',num2str(obj.video.NumberOfFrames)]);
            end
            
            contour_raw=zeros(NF,obj.Nrad);
            contour_fine=zeros(NF,obj.Nrad);
            cen=zeros(NF,2);
            rc=cell(1,NF);
            
            fluo=0;
            plotting=0;
            method='profile';
            if nargin==2
                plotting=varargin{1};
            elseif nargin==3
                plotting=varargin{1};
                fluo=varargin{2};
            end
            
            

            
            av30=fspecial('average',30);
            ga=fspecial('gaussian',7,3);
            if fluo
                obj.frame=medfilt2(obj.frame,[7,7]);
                obj.frame=imfilter(obj.frame,ga,'replicate');
                obj.frame=-obj.frame+max(obj.frame(:));
                obj.frame=obj.frame+imfilter(max(obj.frame(:))-obj.frame,av30,'replicate');
                

            end
            %basic initialization for checks
            obj.ring=get_ring_interpolation(obj.frame,obj.Nrad,obj.radius_limits,obj.center);
            template=corr_template(obj.ring,method);
            [~,firstrc]=find_contour(obj.ring,obj.radius_limits,template);
            cc=obj.radius_limits(1)-fit_contour(firstrc)-floor(length(template)/4);
            
            obj.update_info(obj.ring,cc);
            
            %plotting the results every 60 frames
            if plotting
                figure(2)
                clf
                IM_handle=imagesc(obj.frame);
                colormap(gray)
                axis image
                hold on;
                CEN_handle=plot(obj.center(1),obj.center(2),'r.');
                CON_handle=plot(obj.xc+obj.center(1),obj.yc+obj.center(2),'r-');
                TITLE_handle=title('');
                drawnow
            end
            
            
            disp('Wait for contour calculation');
            tic
            t0=toc;
            
            %dividing the dataset into N-frames bunches
            step=30;
            range=unique([1:step:NF,NF]);
            
            % some strange redefinement of variables, used to reduced the
            % amount of bytes passed to each processor when the algorithm
            % is used in parallel
            circleinfo.sigma=obj.sigma;
            circleinfo.radius_limits=obj.radius_limits;
            circleinfo.Nrad=obj.Nrad;
            circleinfo.center=obj.center;
            cen(1,:)=obj.center;
            
            ang=linspace(0,2*pi,obj.Nrad);
            for nn=1:length(range)-1;
                
                % removing unnecessary color info (videos are in BW) by
                % selecting the red channel
                ffs=read(obj.video,[range(nn),range(nn+1)]);
                if length(size(ffs))>3
                    ffs=squeeze(ffs(:,:,1,:));
                end
                try
                    parfor k=range(nn):range(nn+1);
                        
                        sigma=circleinfo.sigma; %#ok<*PFBNS>
                        radius_limits=circleinfo.radius_limits;
                        Nrad=circleinfo.Nrad;
                        center=circleinfo.center;
                        
                        frame=double(ffs(:,:,k-range(nn)+1));
                        if fluo
                            frame=medfilt2(frame,[7,7]);
                            frame=imfilter(frame,ga,'replicate');
                            frame=-frame+max(frame(:));
                            frame=frame+imfilter(max(frame(:))-frame,av30,'replicate');
                        end
                        
                        
                        %executing the raw fitting procedure until the center
                        %position is definitively found.
                        for j=1:3
                            contour_raw(k,:)=find_contour(get_ring_interpolation(frame,Nrad,radius_limits,center),radius_limits,template);
                            
                            [xc,yc]=pol2cart(ang,contour_raw(k,:));%#ok<*PROP>
                            mean_radius=mean(contour_raw(k,:));
                            radius_limits=sigma*mean_radius*[1,-1]+mean_radius;
                            center=center+[mean(xc),mean(yc)];
                            
                        end
                        
                        [contour_raw(k,:),rc{k}]=find_contour(get_ring_interpolation(frame,Nrad,radius_limits,center),radius_limits,template); %#ok<*PFOUS>
                        
                        contour_fine(k,:)=radius_limits(1)-fit_contour(rc{k})-floor(length(template)/4);
                        
                        [xc,yc]=pol2cart(ang,contour_fine(k,:));
                        cen(k,:)=center+[mean(xc),mean(yc)];
                    end
                catch err
                    contour_fine=contour_raw(1:range(nn)-1,:);
                    contour_raw=contour_raw(1:range(nn)-1,:);
                    cen=cen(1:range(nn)-1,:);
                    disp(' ')
                    disp('Analysis interrupted by an undesired error.')
                    disp('The results calculated so far have been saved.')
                    break
                end
                
                ndone=nn*step;
                if ndone>NF
                    ndone=NF;
                end
                str=['FPS: ',num2str(size(ffs,3)/(toc-t0),'%.1f'),'  -  ',num2str(100*ndone./NF,'%.1f'),'%% completed'];
                t0=toc;
                
                if plotting
                    cc=cen(range(nn+1),:);
                    [xc,yc]=pol2cart(ang,contour_fine(range(nn+1),:));
                    set(IM_handle,'CData', ffs(:,:,end));
                    set(CEN_handle,'XData',cc(1),'YData',cc(2));
                    set(CON_handle,'XData',xc+cc(1),'YData',yc+cc(2));
                    set(TITLE_handle,'string',str);
                    drawnow;
                else
                    if nn>1
                        fprintf(repmat(char('\b'),1,Lstr));
                    end
                    fprintf(str);
                    Lstr=length(str)-1;
                end
            end
            obj.frame=get_frame(obj.video,size(contour_fine,1));
            obj.current_frame=size(contour_fine,1);
            
            disp(' ');
            disp('Done!')
            obj.contour_fine=contour_fine;
            obj.contour_raw=contour_raw;
            obj.cen=cen;
            
            figure(11)
            imagesc(obj.contour_fine')
            axis image
            xlabel('frame #')
            ylabel('azimuth angle #')
        end
        
        function obj=analyse_contour_raw(obj,varargin)
            % Function that detects the contour of the vesicle in every
            % frame of the movie. The video is analysed frame-by-frame
            % using a fast but raw fitting procedure.
            %
            % USAGE:
            %   obj.analyse_contour_raw;
            %   obj.analyse_contour_raw(plotting);
            %   obj.analyse_contour_raw(plotting,method);
            %
            % plotting: boolean variable to activate plotting (default: true)
            %
            % method:  "profile" (default) the contour is detected by
            %           correlating each radial profile with a previously
            %           determined template.
            
            contour_raw=zeros(obj.video.NumberOfFrames,obj.Nrad);
            cen=zeros(obj.video.NumberOfFrames,2);
            
            fluo=0;
            plotting=0;
            method='profile';
            if nargin==2
                plotting=varargin{1};
            elseif nargin==3
                plotting=varargin{1};
                fluo=varargin{2};
            end
            
            
            NF=obj.video.NumberOfFrames;
            
            
            if fluo
                %                     obj.frame=imfilter(obj.frame,ones(3,3));
                obj.frame=-obj.frame+max(obj.frame(:));
            end
            %basic initialization for checks
            for fooind=1:3
                obj.ring=get_ring_interpolation(obj.frame,obj.Nrad,obj.radius_limits,obj.center);
                template=corr_template(obj.ring,method);
                [cc,~]=find_contour(obj.ring,obj.radius_limits,template);
                obj.update_info(obj.ring,cc);
            end
            
            %plotting the results every 60 frames
            if plotting
                figure(2)
                clf
                IM_handle=imagesc(obj.frame);
                colormap(gray)
                axis image
                hold on;
                CEN_handle=plot(obj.center(1),obj.center(2),'r.');
                CON_handle=plot(obj.xc+obj.center(1),obj.yc+obj.center(2),'r-');
                TITLE_handle=title('');
                xlim(round(obj.center(1)+obj.mean_radius*[-2,2]))
                ylim(round(obj.center(2)+obj.mean_radius*[-2,2]))
                drawnow
            end
            
            
            disp('Wait for contour calculation');
            tic
            t0=toc;
            
            
            for nn=1:NF;
                
                % removing unnecessary color info (videos are in BW) by
                % selecting the red channel
                obj.frame=get_frame(obj.video,nn);
                if fluo
                    %                     obj.frame=imfilter(obj.frame,ones(3,3));
                    obj.frame=-obj.frame+max(obj.frame(:));
                end
                obj.current_frame=nn;
                
                try
                    for fooind=1:3
                        obj.ring=get_ring_interpolation(obj.frame,obj.Nrad,obj.radius_limits,obj.center);
                        %                         if fooind==1
                        %                             template=corr_template(obj.ring,method);
                        %                         end
                        [contour_raw(nn,:),~]=find_contour(obj.ring,obj.radius_limits,template);
                        obj.update_info(obj.ring,contour_raw(nn,:));
                        
                        cen(nn,:)=obj.center;
                    end
                    if isnan(contour_raw(nn,:))
                        contour_raw=contour_raw(1:nn-3,:);
                        cen=cen(1:nn-3,:);
                        break
                    end
                catch err
                    contour_raw=contour_raw(1:nn-1,:);
                    cen=cen(1:nn-1,:);
                    disp(' ')
                    disp('Analysis interrupted by an undesired error.')
                    disp('The results calculated so far have been saved.')
                    break
                end
                
                str=['FPS: ',num2str(1/(toc-t0),'%.1f'),'  -  ',num2str(100*nn./NF,'%.1f'),'%% completed'];
                t0=toc;
                
                if plotting
                    set(IM_handle,'CData', obj.frame);
                    set(CEN_handle,'XData',obj.center(1),'YData',obj.center(2));
                    set(CON_handle,'XData',obj.xc+obj.center(1),'YData',obj.yc+obj.center(2));
                    set(TITLE_handle,'string',strrep(str,'%%','%'));
                    drawnow;
                else
                    if nn>1
                        fprintf(repmat(char('\b'),1,Lstr));
                    end
                    fprintf(str);
                    Lstr=length(str)-1;
                end
            end
            obj.frame=get_frame(obj.video,size(contour_raw,1));
            obj.current_frame=size(contour_raw,1);
            
            disp(' ');
            disp('Done!')
            
            rr=contour_raw;
            rr=rr-repmat(mean(rr,2),1,size(rr,2))+mean(rr(:));
            
            
            obj.contour_raw=rr;
            obj.cen=cen;
            
            figure(11)
            imagesc(obj.contour_raw')
            axis image
            xlabel('frame #')
            ylabel('azimuth angle #')
        end
        
        function obj=get_spectrum(obj,varargin)
            %function that analyse the contour_fine property and extract
            %the relevant physical parameters of the membrane, and it stores
            %them in the "results" property.
            %
            % USAGE:
            %   obj.get_spectrum;
            %   obj.get_spectrum(plotting);
            %   obj.get_spectrum(plotting,frame_to_analyze);
            %
            % plotting: boolean value to trigger plotting (default=true)
            %
            % frame_to_analyze: array containing the frames that have to be
            %                   considered into this analysis. By default, the whole video is
            %                   considered. This parameter is useful if some contour shows
            %                   defects.
            
            plotting=false;
            frame_range=1:size(obj.contour_fine,1);
            cc=obj.contour_fine;
            if nargin==2
                plotting=varargin{1};
            elseif nargin==3
                plotting=varargin{1};
                frame_range=varargin{2};
            elseif nargin==4
                plotting=varargin{1};
                frame_range=varargin{2};
                cc=varargin{3};
            end
            
            data.T=obj.Temp;
            
            
            dpix=obj.dpix;
            if dpix>0.1
                dpix=dpix*1e-6;
            end
            dt=obj.dt;
            clog=cc(frame_range,:);
            clear cc
            
            
            if numel(clog)==0
                error('Contour matrix is empty. Run the "analyse_contour" routine!');
            end
            
            %data initialization
            rr=clog*dpix;
            data.R=mean(rr(:));
            
            %perimeter calculation
            ang=linspace(0,2*pi,obj.Nrad)';
            [x,y]=pol2cart(repmat(ang,1,size(rr,1)),rr');
            dx=x-circshift(x,[-1,0]);
            dy=y-circshift(y,[-1,0]);
%             data.L=mean(sum(sqrt(dx.^2+dy.^2),1)); % real perimeter
            data.L=2*pi*data.R;                      % 2 pi R
            
            % subtracting the mean contour to all the contours.
            rr=rr-repmat(mean(rr,2),1,size(rr,2));
            rr=rr-repmat(mean(rr,1),size(rr,1),1);
            
            % STATIC power spectrum, averaged over frames.
            FFTq=fft(rr,[],2);
            normq=size(rr,2).^2;
            FFTqm=mean( FFTq.*conj(FFTq)./normq ,1);
            data.static.ps = FFTqm(1:end/2)'*4; %il per 4 viene da un test su un sin
            %             data.static.q  = 1/(data.R)*(0:length(data.static.ps)-1)';
            data.static.q = 2*pi*obj.Nrad/(data.L)*linspace(0,1,obj.Nrad)';
            data.static.q = data.static.q(1:end/2);
            
            data.static.range=data.static.q>data.static.q(5)&data.static.q<data.static.q(22);
            data.static.err=std( FFTq.*conj(FFTq)./normq,0,1);
            data.static.err=data.static.err(1:end/2);
            
            % DYNAMIC power spectrum, averaged over modes.
            FFTw=fft(rr,[],1);
            normw=size(rr,1).*(2*pi/dt);
            FFTwm=mean( FFTw.*conj(FFTw)./normw ,2);
            data.dynamic.ps=FFTwm(1:floor(length(FFTwm)/2));
            data.dynamic.omega=2*pi/dt*(1:length(data.dynamic.ps))'/length(data.dynamic.ps);
            
            data.dynamic.range=data.dynamic.omega>data.dynamic.omega(2)&data.dynamic.omega<data.dynamic.omega(3)*30;
            data.dynamic.err=std( FFTw.*conj(FFTw)./normw ,0,2);
            data.dynamic.err=data.dynamic.err(1:floor(length(FFTwm)/2));
            
            % Autocorrelation of the temporal evolution of the static spectrum's modes
            plot_autocor=false; %true for debugging
            tocorr=FFTq.*conj(FFTq)./normq;
            data.corr.g2=autocor( tocorr',dt,plot_autocor);
            data.corr.g2=data.corr.g2(2:end,:); % the first point is crappy.
            
            %fitting with a simple exponential, and obtaining the relaxation time of
            %the modes.
            Fexp=fittype('a+b*exp(-x./tau)');
            options=fitoptions( 'Method','nonLinearLeastSquares',...
                'startpoint',[1,1,0.1],'lower',[0.9,0.5,1e-3],'upper',[2,2,10]);
            
            qrange=(2:30);
            j=0;
            rt=data.corr.g2(:,1)<1000*dt;
            data.corr.tau=zeros(length(qrange),1);
            data.corr.a=zeros(length(qrange),1);
            data.corr.b=zeros(length(qrange),1);
            data.corr.fit_g2=cell(length(qrange),1);
            
            for ii=qrange
                j=j+1;
                
                data.corr.fit_g2{j}=fit(data.corr.g2(rt,1),data.corr.g2(rt,ii+1),Fexp,options);
                data.corr.tau(j)=data.corr.fit_g2{j}.tau;
                data.corr.a(j)=data.corr.fit_g2{j}.a;
                data.corr.b(j)=data.corr.fit_g2{j}.b;
            end
            
            data.corr.err=linspace(0.01,0.1,length(qrange))'.*data.corr.tau;
            data.corr.q=data.static.q(qrange);
            data.corr.range=(1:5);
            
            obj.spectrum=data;
            
            if plotting
                warning('OFF','MATLAB:Axes:NegativeDataInLogAxis');
                
                figure(3)
                clf
                errorbar(data.static.q(2:end),data.static.ps(2:end),data.static.err(2:end),'r.')
                set(gca,'yscale','log','xscale','log')
                xlabel('q (m^{-1})')
                ylabel('<|u(q_x)|^2> (m^2)')
                
                figure(4)
                clf
                errorbar(data.dynamic.omega(3:end),data.dynamic.ps(3:end),data.dynamic.err(3:end),'r.')
                set(gca,'yscale','log','xscale','log')
                xlabel('\omega (rad/s)')
                ylabel('PSD (m^2/s)')
                
                figure(5)
                clf
                loglog(data.corr.q,data.corr.tau,'ro')
                xlabel('q (m^{-1})')
                ylabel('\tau (s)')
                
                drawnow
            end
        end
        
        function obj=fit_vesicles_fluctuations(obj,varargin)
            plotting=false;
            if nargin==2
                plotting=varargin{1};
            end
            data=obj.spectrum;
            
            % fitting the static spectrum
            start=[-8,-19];
            lower=[-13,-25];
            upper=[-4,-16];
            modelPSDq=fittype(@(s,k,x) log10(PSDq(s,k,0,data.T,data.L,x)));
            options=fitoptions('Method','nonlinearleastsquares',...
                'StartPoint',start,'lower',lower,'upper',upper);%,...
            %                 'Weight',1./(data.static.err(data.static.range)'./data.static.ps(data.static.range).*log(10)));
            [fitPSDq,Gs]=fit(data.static.q(data.static.range),...
                log10(data.static.ps(data.static.range)),...
                modelPSDq,options);
            confidence_intervals=confint(fitPSDq,0.68);
            
            data.k=10.^fitPSDq.k;
            data.dk=mean(abs(10.^confidence_intervals(:,2)-data.k));
            data.s=10.^fitPSDq.s;
            data.ds=mean(abs(10.^confidence_intervals(:,1)-data.s));
            data.gamma=0;
            
            %             % fitting the dynamic spectrum
            %             start=[fitPSDq.s,fitPSDq.k,-2];
            %             lower=[fitPSDq.s-0.5,fitPSDq.k-0.5,-3];
            %             upper=[fitPSDq.s+0.5,fitPSDq.k+0.5,-1];
            %             modelPSDw=fittype(@(s,k,eta,x) log10(PSDw(s,k,eta,data.T,data.R,x)));
            %             options=fitoptions('Method','nonlinearleastsquares',...
            %                                'StartPoint',start,'lower',lower,'upper',upper,...
            %                                'Weights',data.dynamic.ps(data.dynamic.range)./(data.dynamic.err(data.dynamic.range)));
            %             [fitPSDw,Gw]=fit(data.dynamic.omega(data.dynamic.range),log10(data.dynamic.ps(data.dynamic.range)),...
            %                       modelPSDw,options);
            %             kd=10.^fitPSDw.k;
            %             sd=10.^fitPSDw.s;
            %             data.eta=10.^fitPSDw.eta;
            
            
            start=[-2];
            lower=[-5];
            upper=[-1];
            modeltauq=fittype(@(eta,x) tauq(eta,log10(data.s),log10(data.k),0,data.R,x));
            options=fitoptions('Method','nonlinearleastsquares',...
                'StartPoint',start,'lower',lower,'upper',upper);
            [fittau,Gt]=fit(data.corr.q(data.corr.range),data.corr.tau(data.corr.range),...
                modeltauq,options);
            data.eta=10.^fittau.eta;
            data.fit_goodness=Gs.sse+Gt.sse;
            
            
            
            ms=PSDq(log10(data.s),log10(data.k),data.gamma,data.T,data.R,data.static.q);
            data.dev_from_model=sum(abs(data.static.ps(2:end)-ms(2:end))./ms(2:end))./length(ms);
            
            
            obj.results=data;
            
            % PLOTTING THE RESULTS % ######################################
            if plotting
                warning('OFF','MATLAB:Axes:NegativeDataInLogAxis');
                r=find(data.static.range);
                qs=logspace(log10(data.static.q(r(1))),log10(data.static.q(r(end))),1000)';
                modelstatic=PSDq(log10(data.s),log10(data.k),data.gamma,data.T,data.L,qs);
                
                om=logspace(log10(data.dynamic.omega(1)),log10(data.dynamic.omega(end)),1000)';
                modeldynamic=PSDw(log10(data.s),log10(data.k),log10(data.eta),data.T,data.R,om);
                
                qt=logspace(log10(data.corr.q(1)/2),log10(data.corr.q(end)*2),1000)';
                modelcorr=tauq(log10(data.eta),log10(data.s),log10(data.k),data.gamma,data.R,qt);
                
                figure(3)
                clf
                errorbar(data.static.q(2:end),data.static.ps(2:end),data.static.err(2:end),'r.')
                hold on
                plot(qs,modelstatic,'k-','linewidth',2)
                set(gca,'yscale','log','xscale','log')
                xlabel('q (m^{-1})')
                ylabel('<|u(q_x)|^2> (m^2)')
                
                %                 figure(4)
                %                 clf
                %                 errorbar(data.dynamic.omega(3:end),data.dynamic.ps(3:end),data.dynamic.err(3:end),'r.')
                %                 hold on
                %                 plot(om,modeldynamic,'k-','linewidth',2)
                %                 set(gca,'yscale','log','xscale','log')
                %                 xlabel('\omega (rad/s)')
                %                 ylabel('PSD (m^2/s)')
                %
                %                 figure(5)
                %                 clf
                %                 loglog(data.corr.q,data.corr.tau,'ro')
                %                 hold on
                %                 loglog(qt,modelcorr,'k-','linewidth',2)
                %                 xlabel('q (m^{-1})')
                %                 ylabel('\tau (s)')
                %
                drawnow
            end
        end
        
        
    end
    
    methods (Hidden)
        % Hidden function that uptades the structure after each
        
        function obj=update_info(obj,ring,contour)
            obj.ring=ring;
            ang=linspace(0,2*pi,length(contour));
            [obj.xc,obj.yc]=pol2cart(ang,contour);
            
            obj.mean_radius=mean(contour);
            obj.radius_limits=obj.sigma*obj.mean_radius*[1,-1]+obj.mean_radius;
            obj.center=obj.center+[mean(obj.xc),mean(obj.yc)];
            
        end
    end
    
    
end

%% BASIC CONTOUR FITTING (PARABOLIC LOCAL FIT)
function frame=get_frame(movie,index)
%loading a frame from the movie object
frame=double(read(movie,index));
if length(size(frame))>2
    frame=squeeze(double(frame(:,:,1)));
end
end

function show_img(fig,img)
% a basic function to plot the image with the correct map and axis
if fig>0
    figure(fig)
    clf
    
end
imagesc(img)
colormap gray
axis image off
end

function II=clean_cont(II,lim,lmax)
% NB=smooth(abs(gradient(II)),3)'>2;
% II(NB)=round(mean(II(~NB)));
% 
% 
% PB_low=II>=lim+2;
% PB_high=II<=(lmax-(lim+2));
% 
% II(~PB_low)=lim+2;
% II(~PB_high)=lmax-(lim+2);

NB=(abs(gradient(II))>2);
NB=NB+(abs(gradient(gradient(II)))>1);
NB=NB+(abs(gradient(gradient(gradient(II))))>0.5);
NB=NB>1;

II(NB)=ceil(nanmean(II(~NB)));

end

function [vx,vy,a]=lsq_parab_vertex(x,y)
%least square parabolic fitting
s0 = length(x); s1 = sum(x); s2 = sum(x.^2); s3 = sum(x.^3); s4 = sum(x.^4);
A = [s4,s3,s2;s3,s2,s1;s2,s1,s0];
d = [sum(x.^2.*y);sum(x.*y);sum(y)];
a = A\d;

vx=-a(2)./2./a(1);
vy = a(1)*vx.^2+a(2)*vx+a(3);
end


function ring=get_ring_interpolation(frame,Nrad,radius_limits,center)
%generating a coordinate array of the ring that contains the
%membrane
ang=linspace(0,2*pi,Nrad);
Np=round(abs(diff(radius_limits)));
rad=linspace(radius_limits(1),radius_limits(2),Np);
[theta,rho]=meshgrid(ang,rad);
[Xp,Yp]=pol2cart(theta,rho);

%generating an interpolant function of the whole frame
[X,Y]=meshgrid(1:size(frame,2),1:size(frame,1));
F=griddedInterpolant(X',Y',frame','linear');

%interpolating (30 times faster than multiple line interpolation);
ring=F(Xp+center(1),Yp+center(2));

end

function template=corr_template(r1,method)


NN=8;
%%gradient template
if strcmp(method,'gradient')
    template=(NN:-1:-NN)';
elseif strcmp(method,'profile')
    r1(isnan(r1))=nanmean(r1(:));
    
    filt=(max(fspecial('gaussian',30,2),[],2));
    foo=conv2(r1,filt,'same');
    foo=foo./repmat(max(foo,[],1),size(foo,1),1);
    
    gg=foo';
    [~,II]=max(gg,[],2);
%     mask=make_mask(r1);  %%% manual MASK
    mask=ones(size(r1));
    
    % NB=II<30|II>size(r1,1)-30;
    % mask(:,NB)=NaN;
    % mask=mask.*fliplr(mask);
    
    %masking macroscopic errors
    

    RR=floor(size(r1,1));
    reg=(1:RR)';
    
    
    % ST=size(r1,1)-II';
    ST=RR-II';
    
    IND=repmat(ST,length(reg),1)+repmat(reg,1,size(r1,2));
    INDX=repmat(1:size(IND,2),size(IND,1),1);
    
    rrsh=nan(max(IND(:)),size(IND,2));
    fooind=sub2ind(size(rrsh),IND,INDX);
    rrsh(fooind)=r1.*mask;
    % MN=repmat(min(rrsh,[],1),size(rrsh,1),1);
    % MX=repmat(max(rrsh,[],1),size(rrsh,1),1);
    % rrsh=(rrsh-MN)./(MX-MN);
    
%     show_img(1,rrsh)
    
    
    template=nanmean(rrsh,2);
    template=template(~isnan(template));
    [~,mmx]=max(abs(gradient(template)));
%     template=template(max([mmx-8,1]):min([mmx+8,length(template)]));

    L2=round(length(template)/2);
    template=template(max([L2-8,1]):min([L2+8,length(template)]));
end
end


function maskd=make_mask(m)

figure(1)
clf
imagesc(m);
axis image
title({'MASKING the region that you want to exclude from the template calcularion';'Press Shift+s to save the mask'; 'or Shift+e to exit'})
set(gca,'Units','normalized')
[sx sy]=size(m);
maskd=ones(size(m));
while 1
    k=waitforbuttonpress;
    if strcmp(get(gcf,'CurrentCharacter'),'S')==1
        disp('END')
        break;
    end
    if strcmp(get(gcf,'CurrentCharacter'),'E')==1
        disp('EXIT')
        break;
    end
    p1=get(gca,'CurrentPoint');
    finrct=rbbox;
    p2=get(gca,'CurrentPoint');
    %m(p1(1,2):p2(2,2),p1(1,1):p2(1,1))=0;
%     x1=fix(min([p1(1,2) p2(2,2)]));
%     x2=fix(max([p1(1,2) p2(2,2)]));
    x1=1;
    x2=sx;
    y1=fix(min([p1(1,1) p2(1,1)]));
    y2=fix(max([p1(1,1) p2(1,1)]));
     
    if y1<1
        y1=1;
    end
    
    if y2>sy
        y2=sy;
    end
    
    maskd(x1:x2,y1:y2)=NaN;
%     maskd=maskd.*fliplr(maskd);
    
    imagesc(m.*maskd);
    axis image;
    title({'Press Shift+s to save the mask'; 'or Shift+e to exit'})
end

close(1)
end



function [contour,ringcorr]=find_contour(rr,radius_limits,template)

rr(isnan(rr))=nanmean(nanmean(rr));

%correlation matrix (template matching)
c=conv2(rr,template,'same');
h=conv2(rr,ones(size(template)),'same');
cx=c./h;
cx=cx./max(max(cx));


L=length(template);
L2=floor(L/2);
ringcorr=cx(L2:end-L2,:);

% finding the maximum correlation using a parabolic fitting
% procedure.
[~,II]=nanmin(ringcorr);
II=clean_cont(II,5,size(ringcorr,1));
if numel(II)==0
    II=round(size(ringcorr,1)/2)*ones(1,size(ringcorr,2));
end
p=zeros(1,size(ringcorr,2));


for k=1:size(ringcorr,2);
    xx=II(k)+(-5:1:5)';
    xx=xx(xx>0&xx<size(ringcorr,1));
    yy=ringcorr(xx,k);
    p(k)=lsq_parab_vertex(xx,yy);
end
isnan_parab=isnan(mean(p,1));
p(:,isnan_parab)=repmat(nanmean(p,2),1,sum(isnan_parab));

contour=radius_limits(1)-p-L2/2;


end


%% FINE CONTOUR GLOBAL FITTING ###########################
function [c,c0,ee]=fit_contour(ringcorr)

[~,II]=nanmin(ringcorr);
II=clean_cont(II,5,size(ringcorr,1)-5);

L=length(II);

x=(-5:1:5)';

indx=bsxfun(@plus,II,x);
indx(indx<1)=1;
indx(indx>size(ringcorr,1))=size(ringcorr,1);
indy=repmat(1:L,length(x),1);
indy(indy>size(ringcorr,2))=size(ringcorr,2);

IND=sub2ind(size(ringcorr),indx,indy);

rc=ringcorr(IND);

rcnorm=bsxfun(@rdivide,rc,min(rc,[],1))-1;

x0=zeros(1,L);
parab=zeros(3,L);
for k=1:L;
    [x0(k),~,parab(:,k)]=lsq_parab_vertex(x,rcnorm(:,k));
end
x0(isnan(x0))=nanmean(x0);
isnan_parab=isnan(mean(parab,1));
parab(:,isnan_parab)=repmat(nanmean(parab,2),1,length(find(isnan_parab)));



nn=smooth(abs(gradient(x0+II)),3)'>1;
pnn=find(nn);
for k=pnn
    x0(k)=mean(x0(~nn));
    II(k)=round(mean(II(~nn)));
    parab(:,k)=mean(parab(:,~nn),2);
end


xbest=x0;
N=36;
options=optimset('fminsearch');
options=optimset(options,'MaxFunEvals',320,'MaxIter',80,'Display','off');
% options=optimset('fminunc');
% options=optimset(options,'Algorithm','interior-point');
nn=0;
for k=[1:N,1:N];
    if k==1
        nn=nn+1;
    end
    range=mod((k-1)*L/N+(1:L/N)+N/2*(nn-1),L);
    range(range==0)=L;
    f=@(x) energy(x,II,parab,range,x0);
%         [xb,ee]=fminunc(f,xbest(range),options);
    [xb,ee]=fminsearch(f,xbest(range),options);
%     [xb,ee]=patternsearch(f,xbest(range));
    xbest(range)=xb;
end

c=II+xbest;
c0=II+x0;
% c=clean_cont(c,5,size(ringcorr,1));
end

function ee=energy(xb,II,parab,range,x0)
x=x0;
x(range)=xb;
curv=curvature(II+x);
k=4e-1;

ee=sum(parab(1,:).*x.^2+parab(2,:).*x+parab(3,:)+k*curv);
end

function d2=curvature(cn)
c1=[cn(3:end),cn(1:2)];
c2=[cn(2:end),cn(1)];
c3=[cn(end),cn(1:end-1)];
c4=[cn(end-1:end),cn(1:end-2)];

d2=abs(-c1+16*c2-30*cn+16*c3-c4)/12;
end


%% DATA ANALYSIS ##########################################

function chisq = globalfit(par,data)
k=par(1);
kd=par(2);
s=par(3);
gamma=par(4);
eta=par(5);

modelstatic=PSDq(s,k,gamma,data.T,data.R,data.static.q(data.static.range));
modeldynamic=PSDw(s,kd,eta,data.T,data.R,data.dynamic.omega(data.dynamic.range));
modelcorr=tauq(eta,s,kd,gamma,data.R,data.corr.q(data.corr.range));

chiq=(log10(data.static.ps(data.static.range))-log10(modelstatic))./length(modelstatic);
chiw=(log10(data.dynamic.ps(data.dynamic.range))-log10(modeldynamic))./length(modeldynamic);
chit=(log10(data.corr.tau(data.corr.range))-log10(modelcorr))./length(modelcorr);

chiq=chiq./(data.static.err(data.static.range)./data.static.ps(data.static.range));
chiw=chiw./(data.dynamic.err(data.dynamic.range)./data.dynamic.ps(data.dynamic.range));
chit=chit./(data.corr.err(data.corr.range)./data.corr.tau(data.corr.range));

% chiq=(data.static.ps(data.static.range)-modelstatic)./data.static.err(data.static.range)./sqrt(length(modelstatic));
% chiw=(data.dynamic.ps(data.dynamic.range)-modeldynamic)./data.dynamic.err(data.dynamic.range)./sqrt(length(modeldynamic));
% chit=(data.corr.tau(data.corr.range)-modelcorr)./data.corr.err(data.corr.range)./sqrt(length(modelcorr));

chisq=[chiq;chiw;chit];

end

function model=tauq(eta,s,k,g,R,qm)
eta=10.^eta;
s=10.^s;
k=10.^k;
etaM=1e-9;
etaext=1e-3;
model=2*(etaM/R^2+qm*(eta+etaext))./(2*g+s*qm.^2+k*qm.^4);

end

function spectrum=PSDq(s,k,g,T,L,x)
s=10.^s;
k=10.^k;
kb=1.38e-23; % SI

spectrum=kb*T/(L) *sqrt(k/(2*(s^2-4*k*g))) *...
    ( 1./sqrt(2*k*x.^2+s - sqrt(s^2-4*k*g)) - 1./sqrt(2*k*x.^2+s +  sqrt(s^2-4*k*g)) );
end

function spectrum=PSDw(s,k,eta,T,R,omega)
s=10.^s;
k=10.^k;
eta=10.^eta;

kb=1.38e-23; % SI

spectrum=zeros(size(omega));
for l=2:130
    vl=-l:l;
    Zlm=Zl(vl);
    sZl=sum(Zlm(~isinf(Zlm)));
    foo=kb.*T./R./eta./sZl.*(2*l+1)/2/pi./(wl(s,k,eta,R,l).^2+omega.^2);
    spectrum=spectrum+foo;
end
end

function wl=wl(s,k,eta,R,l)
wl=(k*(l+2).*(l-1).*l.*(l+1)+ s*R^2*(l+2).*(l-1))./(eta*R^3*Zl(l));
end

function Zl=Zl(l)
Zl=(2*l+1).*(2*l.^2+2.*l-1)./(l.*(l+1));
end


%% AUTOCOR ALGORITHM (multi-tau, vectorized algorithm)
function g2=autocor(I,dwell_time,show)
%	Function autocorrelate the data using multiple tau method.

global buf G num cts cur nobuf nolev

I=I./mean(I(:));

nobuf=8;  %must be even!
nolev=8;

timedelay=delays(dwell_time);
%initialize all the arrays
Ng2=size(I,1);
buf=zeros(nolev,nobuf,Ng2); %matrix of buffers
cts=zeros(nolev,1);
cur=nobuf*ones(nolev,1);
G=zeros((nolev+1)*nobuf/2,Ng2);
num=zeros(nolev,1);


if show
    figure(1)
    clf
    h=loglog(timedelay,G(:,1),'b.-');
    xlim([timedelay(2),timedelay(end)])
    
end
tic;
L=size(I,2);
for n=1:L
    insertimg(I(:,n))
    %timer,el1
    %et=el1-el
    if show
        fprintf('.:: processed frame: %d(%d)  - ',n,L)
        fprintf(' elapsed time: %.4f \r',toc);
    end
    if show && (mod(n,30)==0)
        set(h,'Ydata',G(:,1),'Xdata',timedelay);
        pause(0.01)
    end
    
end
if show
    disp('Done!');
    disp(' ');
end
norm=repmat(mean(I,2).^2,1,length(timedelay))';
g2=[timedelay,G./norm];
end

function dly=delays(time)
%    return array of delays.
%    KEYWORD:  time: scale delays by time ( should be time between frames)

global nolev nobuf

dly=zeros((nolev+1)*nobuf/2,1);
for i=1:nolev
    if i==1
        imin=2;
    else
        imin=nobuf/2+1;
    end
    ptr=(i-1)*nobuf/2+(imin:nobuf);
    dly(ptr)= ((imin-1):(nobuf-1)).*2^(i-1);
end
dly=dly*time;
end

function insertimg(In)
%   read and insert image, n, into accumulator buffers,
%   then increments appropriate correlators for new image.
global buf cts cur nobuf nolev


cur(1)=1+mod(cur(1),nobuf);  %increment buffer

buf(1,cur(1),:)=In;

process(1,cur(1));

processing=1;
lev=2;
while processing
    % either add previous two images then process or set flag next one
    if(cts(lev))
        prev=1+mod((cur(lev-1)-1-1+nobuf),nobuf);
        cur(lev)=1+mod(cur(lev),nobuf);
        buf(lev,cur(lev),:) = (buf(lev-1,prev,:)+buf(lev-1,cur(lev-1),:))./2;
        cts(lev)=0;
        process(lev,cur(lev));
        lev=lev+1;
        %Since this level finished, test if there is a next level for processing
        processing = (lev<=nolev);
    else
        cts(lev)=1; % set flag to process next time
        processing=0; % can stop until more images are accumulated
    end
end
end

function process(lev,bufno)
%	The buffer bufno at level lev has just been added so update
%	the correlation and intensity averages it affects
global buf G num nobuf


num(lev)=num(lev)+1; % one more in average
if (lev==1)
    imin=1;
else
    imin=1+nobuf/2;% no need for lower half of level > 1
end

for i=imin:min([num(lev),nobuf])%loop over delays
    ptr=(lev-1)*nobuf/2+i;
    delayno=1+mod((bufno-(i-1)-1+nobuf),nobuf); %cyclic buffers
    IP=squeeze(buf(lev,delayno,:));
    IF=squeeze(buf(lev,bufno,:));
    foo=(IF.*IP)';
    G(ptr,:) = G(ptr,:)+(foo-G(ptr,:))/(num(lev)-i+1);
end

end


