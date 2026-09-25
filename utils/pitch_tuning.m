function pitch_tuning(ele,clu)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short description
% long description
%
% USAGE:
%       [out] = template(in,in2)
%
% INPUT:
%       in - input 1
%       in2 - input 2
%
% OUTPUT:
%       p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 21/02/23 Initial release

%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK  
    config.outname = 'klustest';
    skipfigs = 0; % 1 = skip making a figure if it exists already    
    save_figs = 1; % 1 = make and save figures (VERY time consuming, takes about 90% of klustest's time)
    fast_figures = 1; % 1 = 4-5x faster figure saving, but at a slightly lower quality
    warning('off','MATLAB:table:RowsAddedExistingVars');

    % settings
    mapset.drive_height_mm = 20;  
    min_dwell = 0.05;
    thresh = 0.5; % threshold for field detection (prop.)
    colormap_dwell = 'viridis';
    colormap_rate = 'turbo';

    % plot 3D tuning curves with depth, i.e. 3D tuning curves will have a 'shape' 
    % rather that just plotted as spherical
    % plot_3D_tuning_curves = 1; 
    plot_3D_tuning_curves = 0; 

    % degrees to rotate 3D tuning curves, if plot_3D_tuning_curves is set to 1, 
    % this should probably be set to 45, otherwise azimuthal cells face directly
    % on to the camera and they are not easy to see
    % offset_3d_tuning_curves = 45;
    offset_3d_tuning_curves = 0;

    % normalise colormaps across all sessions from 0Hz to the max frate of all sessions
    normalise_all_colormaps = 1;

    % 0 to disable, otherwise center sessions to average peak firing rate bin
    % The value of center_all_2D_maps gives the sessions to average    
    % FOr example, [1 2 3] means get the peak angle for sessions 1,2,3 (arena 1,
    % hills, arena 2), average these and shift ALL maps so this average angle is
    % in the center of all plots
    center_all_2D_maps = [1 2 3]; 
    center_all_3D_maps = [1 2 3];    
    centering_method = 4; % 1 = mean peaks, 2 = median peaks, 3 = mean max angle, 4 = median max angle

    % plotting
    % elevation of 3D plots, i.e. do you look 'down' on the 3D tuning curve or
    % not?
    plot_el = 20; 
    plot_field_peak = 1; % set to 1 to plot a marker on the receptive field peak

    if ~exist('ele','var') || isempty(ele)
        ele = [];
    end
    if ~exist('clu','var') || isempty(clu)
        clu = [];
    end

    %% Data formats
    % % if cluster cut with kwiktint
    % formats.pos         = 'Neuralynx';
    % formats.clu         = 'Neuralynx';
    % formats.tet         = 'Tint';
    % formats.set         = 'Neuralynx';    
    % formats.spk         = 'Neuralynx';
    % formats.lfp         = 'Neuralynx';
    % formats.iso         = 'klustakwik';
    
    % if cluster cut with kwikcut
    formats.pos         = 'Neuralynx';
    formats.clu         = 'Kwikcut';
    formats.tet         = 'Kwikcut';
    formats.set         = 'Kwikcut';    
    formats.spk         = 'Kwikcut';
    formats.lfp         = 'Neuralynx';
    formats.iso         = 'Kwikcut';   
    config.cname        = 'kwikcut';
    
    formats.pos         = 'reconstruction';
    formats.front_led_color = 1; % 1 = red, 2 = green, 3 = blue
    formats.back_led_color = 2; % 1 = red, 2 = green, 3 = blue
    formats.led_angle_offset = 0; % CCW offset of LEDs on the head

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PREPARE DATA
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Tetrodes and sessions
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    stk = dbstack;
    tnow = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',stk.name,tnow))

%% >>>>>>>>>> Retrieve the data
    sname = [pwd '\' config.outname '\sdata.mat'];
    disp(sprintf('Loading sdata: %s',sname))    
    load(sname,'sdata'); % load saved session data
    pdata = sdata.Properties.CustomProperties.pdata;
    part_config = pdata.part_config;
    nparts = size(part_config,1);
    disp(sprintf('\t...%d sessions',size(pdata.sessions,1)));   
    disp(sprintf('\t...recording date: %s',pdata.date));            
    disp(sprintf('\t...done'))   

%% >>>>>>>>>> find which tetrodes are available
    % When sessions are cluster cut together (which you should do with multiple sessions recording the same cells)
    % kilocut saves the filenames of the individual sessions in the kilo.mat file as well as the tetrodes analysed
    disp(sprintf('Assessing data...'))
    [~,snames,data_dirs] = get_tets_for_klustest(formats.tet,config);  

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Position data
    disp(sprintf('\t...positions'));
    pos_srate = 50; % desired sampling rate of position data (Hz)

    [pos,~,tstart] = get_pos_for_klustest(formats,data_dirs,snames,pos_srate,mapset); % directories
    pdata.pos = pos;
    pdata.pos_srate = pos_srate;
    pdata.tstart = tstart;

%% >>>>>>>>>> Analyse trajectory and add 3D HD
    disp(sprintf('Loading 3D trajectory...'))   
    if 0
        for pp = 1:nparts % for every part     
            part_now = part_config.part_names{pp};
            if pp==1
                disp(sprintf('\t%s',part_now))                   
            else
                disp(sprintf('\b | %s',part_now))       
            end
    
            % cut the position data to include only this part
            part_times = part_config.part_times{pp};                
            pos = pdata.pos;
            pot = pos.pot;
            pindax = logical(sum(pot' >= part_times(:,1) & pot' <= part_times(:,2),1));
    
            pdata.(part_now).pot = pos.pot(pindax,1); % pos time for this part            
            pdata.(part_now).pox = pos.pox(pindax,1); % pos x for this part
            pdata.(part_now).poy = pos.poy(pindax,1); % pos y for this part
            pdata.(part_now).poz = pos.poz(pindax,1); % pos y for this part            
            pdata.(part_now).yaw = pos.poh(pindax,1); % yaw HD for this part
            pdata.(part_now).pit = pos.pitch(pindax,1); % pitch HD for this part
            pdata.(part_now).rol = pos.roll(pindax,1); % roll HD for this part
        end
        assignin('base','pdata',pdata)
    end
    disp(sprintf('\t...done'))    
    pdata = evalin('base', 'pdata');

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Session data and figure
    % run through every part
    all_pos = [];
    if 0
        % open figure
        fig_dwell = figure('Units','pixels','Position',[100 50 1200 1800],'visible','on');
        set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
        set(gcf,'color','w'); % makes the background colour white 
    
        ann_str = sprintf('Rat: %s, Date: %s, Analysed: %s',sdata.rat{1},sdata.date{1},datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',10,'LineStyle','none','interpreter','none'); 
    
        map_now = 'turbo';

        xsiz = 250;
        xbuff = 130;
        ysiz = xsiz/2;
        ybuff = 100;
        yvec = [780 560 30];
        xnow = 150;

        xt = 1:(nparts+1);
        xvec = [xnow xnow+(xsiz*xt)+(xbuff*xt)];
        set(fig_dwell,'Position',[100 100 max(xvec)+50 950])

        dwell_maps_3d = cell(1,nparts+1);
        for pp = 1:(nparts+1) % for every part               
            if pp==nparts+1
                part_now = 'combined';
            else
                part_now = part_config.part_names{pp};
            end
            if pp==1
                disp(sprintf('\t\t%s',part_now))      
            else
                disp(sprintf('\b, %s',part_now))      
            end
               
            if pp==nparts+1
                % get position data
                ppox = all_pos(:,1); % pos x for this part, in mm
                ppoy = all_pos(:,2); % pos y for this part, in mm
                ppoz = all_pos(:,3); % pos z for this part, in mm  
                ppot = all_pos(:,4); % pos t 
                ppit = all_pos(:,5); % pitch for this part, in rads   
                pyaw = all_pos(:,6); % pitch for this part, in rads
            else
                % get position data
                ppox = double( pdata.(part_now).pox ).*10; % pos x for this part, in mm
                ppoy = double( pdata.(part_now).poy ).*10; % pos y for this part, in mm
                ppoz = double( pdata.(part_now).poz ).*10; % pos z for this part, in mm  
                ppot = double( pdata.(part_now).pot ); % pos t 
                ppit = double( pdata.(part_now).pit ); % pitch for this part, in rads   
                pyaw = double( pdata.(part_now).yaw ); % pitch for this part, in rads
                all_pos = [all_pos; ppox ppoy ppoz ppot ppit pyaw];
            end
    
            % plot trajectory
            ax = axes('Units','pixels','Position',[xvec(pp) yvec(1) xsiz ysiz]);
                dindax = abs([0; diff(pyaw)])>deg2rad(90);
                yaw_plot = pyaw;
                yaw_plot(dindax,:) = NaN;
    
                plot(yaw_plot,ppit,'k'); hold on;
                xlabel(sprintf('Azimuth (%c)',176))
                ylabel(sprintf('Pitch (%c)',176))
                axis xy tight
                daspect([1 1 1])
                view(0,90);
                set(gca,'FontSize',12)
    
                ax.XLim = [-180 180];
                ax.YLim = [-90 90];
                ax.XTick = -180:90:180;
                ax.YTick = -90:45:90;
    
                duration = length(ppot)*(1/50);
                text(ax,1,1,sprintf('%s | %.fs',part_now,duration),'units','normalized','FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom');

            % plot dwell map
            ax = axes('Units','pixels','Position',[xvec(pp) yvec(2) xsiz ysiz]);    
                if isempty(dwell_maps_3d{1,pp})
                    [az_map, tilt_map, F1] = projectTORUS(pyaw, ppit);
                    dwell_maps_3d{1,pp} = F1;
                end   
                F1 = dwell_maps_3d{1,pp};
                F1(F1<min_dwell) = NaN;

                imagesc(az_map(1,:), tilt_map(:,1), F1,'AlphaData',~isnan(F1));
                xlabel(sprintf('Azimuth (%c)',176))
                ylabel(sprintf('Pitch (%c)',176))
                axis xy tight
                daspect([1 1 1])
                set(gca,'FontSize',12)
                ax.CLim = [0,max(F1(:))];
                colormap(gca,colormap_dwell)
                ax.XTick = -180:90:180;
                ax.YTick = -90:45:90;
    
                mval = max(F1,[],"all",'omitmissing');
                text(ax,1,1,sprintf('%.2fs',mval),'units','normalized','FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom');

                if pp==nparts+1
                    axc = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+20 ax.Position(2)+10 12 80]);
                        mat = (linspace(0,max(F1(:)),100))';
                        imagesc([1 1],[min(mat(:)) max(mat(:))],mat);
                        colormap(axc,colormap_dwell);
                        axis xy
        
                        axc.YTick = [];
                        axc.XTick = [];
                        axc.YAxisLocation = 'right';
                        text(0.5,1.35,sprintf('Dwell\ntime (s)'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,1.1,sprintf('Max'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,-0.1,sprintf('0 Hz'),'FontSize',8,'HorizontalAl','center','Units','normalized') 
                end
        end
        assignin('base','dwell_maps_3d',dwell_maps_3d)     

        % Save the figure  
        [~,~,~] = mkdir([pwd '\' pdata.outname '\part_figures']); % create a folder to hold outputs 
        fname = [pwd '\' pdata.outname '\part_figures\pitch_dwell.png'];
        if fast_figures
            frame = getframe(fig_dwell); % fig is the figure handle to save
            [raster, raster_map] = frame2im(frame); % raster is the rasterized image, raster_map is the colormap
            if isempty(raster_map)
                imwrite(raster, fname);
            else
                imwrite(raster, raster_map, fname); % fig_file is the path to the image
            end            
        else
            exportgraphics(fig_dwell,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',350);  
        end
    end
    dwell_maps_3d = evalin('base', 'dwell_maps_3d');

    % close(fig_clust);   
    % keyboard

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Run through clusters
    disp(sprintf('Analysing clusters...'))
    ucis = unique(sdata.uci); % list of unique cells in sdata   
    if ~isempty(ele)
        if ~isempty(clu)
            ucis = unique(sdata.uci(ismember(sdata.tetrode,ele) & ismember(sdata.cluster,clu))); % list of unique cells in sdata   
        else
            ucis = unique(sdata.uci(ismember(sdata.tetrode,ele))); % list of unique cells in sdata   
        end
    end

    %% run through every cell
    loopout = looper(length(ucis));  
    sdata.hd_map_3d = cell(size(sdata,1),1);
    sdata.hd_3d_curves = cell(size(sdata,1),1);
    sdata.hd_3d_info = NaN(size(sdata,1),4);
    
    dwell_maps_2d = cell(2,nparts+1);
    all_dat = [];
    for uu = 1:length(ucis)     
        dat = table;

        % open figure
        fig_clust = figure('Units','pixels','Position',[100 100 1400 1800],'visible','off');
        set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
        set(gcf,'color','w'); % makes the background colour white

        % add an annotation to the figure with some important info
        uci = ucis{uu};
        disp(sprintf('\t%s',uci))  

        idx = find( ismember(sdata.uci,uci) & sdata.partn==1 );
        ann_str = sprintf('Cell: %s, Rat: %s, Date: %s, Tetrode: %d, Cluster: %d, Analysed: %s',sdata.uci{idx},sdata.rat{idx},sdata.date{idx},sdata.tetrode(idx),sdata.cluster(idx),datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',10,'LineStyle','none','interpreter','none');      

        xsiz = 250;
        xbuff = 160;
        ysiz = xsiz/2;
        ybuff = 100;
        yvec = [780 560 30];
        xnow = 150;

        xt = (1:nparts+1);
        xvec = [xnow xnow+(xsiz*xt)+(xbuff*xt)];
        set(fig_clust,'Position',[100 100 max(xvec) 950])

        % run through every part
        all_pos = [];
        all_sindx = [];
        all_axes = cell(5,nparts+1);
        center_angles = NaN(1,nparts+1);
        for pp = 1:(nparts+1) % for every part   
            if pp==nparts+1
                part_now = 'combined';
            else
                part_now = part_config.part_names{pp};
            end
            if pp==1
                disp(sprintf('\t\t%s',part_now))      
            else
                disp(sprintf('\b, %s',part_now))      
            end
               
            dat.rat(pp,1) = sdata.rat(idx);
            dat.date(pp,1) = sdata.date(idx);
            dat.uci(pp,1) = sdata.uci(idx);
            dat.tetrode(pp,1) = sdata.tetrode(idx);
            dat.cluster(pp,1) = sdata.cluster(idx);

            if pp==nparts+1
                % get position data
                ppox = all_pos(:,1); % pos x for this part, in mm
                ppoy = all_pos(:,2); % pos y for this part, in mm
                ppoz = all_pos(:,3); % pos z for this part, in mm  
                ppot = all_pos(:,4); % pos t 
                ppit = all_pos(:,5); % pitch for this part, in rads   
                pyaw = all_pos(:,6); % pitch for this part, in rads
                sindx = all_sindx;
            else
                % get position data
                ppox = double( pdata.(part_now).pox ).*10; % pos x for this part, in mm
                ppoy = double( pdata.(part_now).poy ).*10; % pos y for this part, in mm
                ppoz = double( pdata.(part_now).poz ).*10; % pos z for this part, in mm  
                ppot = double( pdata.(part_now).pot ); % pos t 
                ppit = double( pdata.(part_now).pit ); % pitch for this part, in rads   
                pyaw = double( pdata.(part_now).yaw ); % pitch for this part, in rads
                n_offset = size(all_pos,1);
                all_pos = [all_pos; ppox ppoy ppoz ppot ppit pyaw];
    
                % get spike data
                idx = find( ismember(sdata.uci,uci) & sdata.partn==pp );
                dat.partn(pp,1) = sdata.partn(idx);          
                dat.nspikes(pp,1) = sdata.nspikes(idx);
                dat.frate(pp,1) = sdata.frate(idx);

                sindx = sdata.spt_pot_index{idx};       
                all_sindx = [all_sindx; sindx+n_offset];
            end
            dindax = abs([0; diff(pyaw)])>deg2rad(300);
            yaw_plot = pyaw;
            yaw_plot(dindax,:) = NaN;

            pspx = ppox(sindx);
            pspy = ppoy(sindx);
            pspz = ppoz(sindx);
            pspt = ppot(sindx);
            psyaw = pyaw(sindx);
            pspit = ppit(sindx);

            % plot spikes and positions
            ax = axes('Units','pixels','Position',[xvec(pp) yvec(1) xsiz ysiz]);
                all_axes{1,pp} = gca;
                if isempty(pspt)
                    continue
                end
    
                % triplicate data so we can center plot later
                % positions
                yp = [yaw_plot(:)-360; NaN; yaw_plot(:); NaN; yaw_plot(:)+360];
                tp = [ppit(:); NaN; ppit(:); NaN; ppit(:)];
                plot(yp,tp,'k'); hold on
                % spikes
                yp = [psyaw(:)-360; NaN; psyaw(:); NaN; psyaw(:)+360];
                tp = [pspit(:); NaN; pspit(:); NaN; pspit(:)];
                plot(yp,tp,'r.','MarkerSize',6); hold on;                     
                
                xlabel(sprintf('Azimuth (%c)',176))
                ylabel(sprintf('Pitch (%c)',176))
                axis xy tight
                daspect([1 1 1])
                view(0,90);
                set(gca,'FontSize',10)
    
                ax.XLim = [-180 180];
                ax.YLim = [-90 90];
                ax.XTick = -540:45:540;
                ax.XTickLabel = string(wrapTo180(ax.XTick));
                ax.YTick = -90:45:90;

                duration = length(ppot)*(1/50);
                text(ax,1,1,sprintf('%s | %.fs',part_now,duration),'units','normalized','FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom');

            % plot 3D tuning curve (flat projection)
            ax = axes('Units','pixels','Position',[xvec(pp) yvec(2) xsiz ysiz]);    
                all_axes{2,pp} = gca;
                F1 = dwell_maps_3d{1,pp};                
                [az_map, tilt_map, F2] = projectTORUS(pyaw(sindx), ppit(sindx));
                F3 = F2 ./ F1;
                F3(F1<min_dwell) = NaN;
                torus_stats = quantifyTORUS(az_map, tilt_map, F3, 'threshold', thresh);

                % triplicate data so we can center plot later
                x_plot = [az_map(1,:)'-360; az_map(1,:)'; az_map(1,:)'+360];
                y_plot = [tilt_map(:,1); tilt_map(:,1); tilt_map(:,1)];
                c_plot = [F3 F3 F3];
                imagesc('XData',x_plot,'YData',y_plot,'CData',c_plot,'AlphaData',~isnan(c_plot)); hold on;

                xlabel(sprintf('Azimuth (%c)',176))
                ylabel(sprintf('Pitch (%c)',176))
                axis xy tight
                daspect([1 1 1])
                view(0,90);
                set(gca,'FontSize',10)
                pfrate = max(F3,[],'all','omitmissing');
                ax.CLim = [0,max(0.001,pfrate,'omitmissing')];
                ax.XLim = [-180 180];
                ax.YLim = [-90 90];
                colormap(gca,colormap_rate)
                ax.XTick = [];
                ax.YTick = [];
                box on

                dat.m3D_HD_dwellmap(pp,1) = { F1 };
                dat.m3D_HD_ratemap(pp,1) = { F3 };

                if centering_method==1 || centering_method==2
                    center_angles(1,pp) = torus_stats.peak_azimuth;
                elseif centering_method==3 || centering_method==4
                    center_angles(1,pp) = torus_stats.peak_total_azimuth;
                end

                text(ax,1,1,sprintf('%.2f Hz',pfrate),'units','normalized','FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom');

                if plot_field_peak
                    hold on;
                    plot(ax,[torus_stats.peak_azimuth-360,torus_stats.peak_azimuth,torus_stats.peak_azimuth+360],repmat(torus_stats.peak_tilt,1,3),'+k','MarkerSize',15,'LineStyle','none');
                    xline(ax,[torus_stats.peak_total_azimuth torus_stats.peak_total_azimuth-360 torus_stats.peak_total_azimuth+360],'k');
                end

                if pp==nparts+1
                    axc = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+20 ax.Position(2)+50 12 80]);
                        mat = (linspace(0,max(F3(:)),100))';
                        imagesc([1 1],[min(mat(:)) max(mat(:))],mat);
                        colormap(axc,colormap_rate);
                        axis xy
        
                        axc.YTick = [];
                        axc.XTick = [];
                        axc.YAxisLocation = 'right';
                        text(0.5,1.35,sprintf('Frate\n(Hz)'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,1.1,sprintf('Max'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,-0.1,sprintf('0 Hz'),'FontSize',8,'HorizontalAl','center','Units','normalized') 
                end

            % azimuthal tuning curve
            ax_az = axes('Units','pixels','Position',[ax.Position(1) ax.Position(2)-60 ax.Position(3) 55]);
                all_axes{3,pp} = gca;
                edg = linspace(-180,180,60);
                xi = movmean(edg,2,'EndPoints','discard');
                if isempty(dwell_maps_2d{1,pp})
                    d = histcounts(pyaw,edg);
                    dwell_maps_2d{1,pp} = d;
                else
                    d = dwell_maps_2d{1,pp};
                end
                s = histcounts(pyaw(sindx),edg);
                ratemap = s ./ (d .* (1/pdata.pos_srate));
                xi = movmean(edg,2,'Endpoints','discard');

                % triplicate data so we can center plot later
                e = [xi-360 xi xi+360];
                r = [ratemap ratemap ratemap];
                bar(e,r,1,'k')

                xlabel(sprintf('Azimuth (%c)',176))
                ylabel('Firing Rate (Hz)')  
                set(gca,'FontSize',10)
                box off
                ax_az.YAxisLocation = 'right';
                ax_az.XTick = [-540:45:540];
                ax_az.XTickLabel = string(wrapTo180(ax_az.XTick));
                ax_az.FontSize = 10;
                ax_az.XLim = [-180 180];

                % head direction analyses
                hd3n = ratemap ./ max(ratemap); % normalise cell hd
                hd3n = hd3n(:);            
                mx1 = xi(hd3n == max(hd3n)); % preferred angle (location of max frate)
                if length(mx1)>1
                    mx1 = mx1(1);
                elseif isempty(mx1)
                    mx1 = NaN;
                end
                hold on
                line([mx1-360 mx1-360],ax_az.YLim,'Color','r')
                line([mx1 mx1],ax_az.YLim,'Color','r')
                line([mx1+360 mx1+360],ax_az.YLim,'Color','r')

                % tuning width
                original_length = length(ratemap);
                [peaks,locs,widths,~] = findpeaks(r, 'WidthReference', 'halfprom');

                % filter to keep ONLY the peaks in the middle section
                middle_idx = (locs > original_length) & (locs <= 2 * original_length);
                middle_peaks = peaks(middle_idx);
                middle_widths = widths(middle_idx);

                % find the absolute maximum out of the middle peaks
                [~, max_idx] = max(middle_peaks);
                bin_size = median(diff(edg)); % bin size in degrees
                final_peak_width = middle_widths(max_idx) * bin_size;

                % directional info shuffle
                % observed values
                dmap = (d .* (1/pdata.pos_srate));
                smetric = get_spatial_info(dmap,ratemap,'metrics',{'spatial_info','kld'});
                vals = [smetric.kldivergence smetric.skaggs_si_bits_per_spike];
            
                % shuffle values
                rng(999); % for reproducibility
                iti = 1000;
                vals_shuff = NaN(iti,1);
                for ii = 1:iti
                    [~,sindx2] = simply_shuffle_spike_train(ppot,pspt,'spindx',sindx);
                    s_shuff = histcounts(pyaw(sindx2),edg);
                    ratemap_shuff = s_shuff ./ dmap;
                    smetric = get_spatial_info(dmap,ratemap_shuff,'metrics',{'spatial_info'});
                    vals_shuff(ii) = smetric.skaggs_si_bits_per_spike;
                end
                vals_z = (vals - mean(vals_shuff,1,'omitnan')) ./ std(vals_shuff,1,'omitnan');
                vals_p = (sum(vals_shuff >= vals, 1) + 1) / (iti + 1);

                % accumulate
                dat.m2D_azimuth(pp,1) = { xi };
                dat.m2D_azimuth_dwellmap(pp,1) = { d };
                dat.m2D_azimuth_ratemap(pp,1) = { ratemap };
                dat.m2D_azimuth_pfd(pp,:) = mx1;  
                dat.m2D_azimuth_pfd_frate(pp,:) = max(ratemap(:));  
                dat.m2D_azimuth_avg_frate(pp,:) = mean(ratemap(:),'omitnan');                  
                dat.m2D_azimuth_info(pp,:) = vals;
                dat.m2D_azimuth_info_zscored(pp,:) = vals_z;
                dat.m2D_azimuth_info_p(pp,:) = vals_p;
                dat.m2D_azimuth_tuning_width(pp,:) = final_peak_width;
               
                text(ax,0,-1.1,sprintf('Azimuth (N shuff = %d)\nspatial info %.2f b/spike (z = %.2f, p = %.3f)\nTuning width: %.2f%c',iti,vals(1),vals_z(1),vals_p(1),final_peak_width,176),'units','normalized','FontSize',10);

            % pitch tuning curve
            ax_pi = axes('Units','pixels','Position',[ax.Position(1)-60 ax.Position(2) 55 ax.Position(4)]);
                all_axes{4,pp} = gca;
                edg = linspace(-90,90,60);
                xi = movmean(edg,2,'EndPoints','discard');
                if isempty(dwell_maps_2d{2,pp})
                    d = histcounts(ppit,edg);
                    dwell_maps_2d{2,pp} = d;
                else
                    d = dwell_maps_2d{2,pp};
                end
                s = histcounts(ppit(sindx),edg);
                ratemap = s ./ (d .* (1/pdata.pos_srate));
                xi = movmean(edg,2,'Endpoints','discard');
                                
                barh(xi,ratemap,1,'k')
                ylabel(sprintf('Pitch (%c)',176))
                xlabel('Firing Rate (Hz)')  
                set(gca,'FontSize',10)
                box off
                ax_pi.XAxisLocation = 'top';
                ax_pi.YTick = [-90:45:90];
                ax_pi.FontSize = 10;
                ax.YLim = [-90 90];
                ax.XLim = [-180 180];

                % head direction analyses
                hd3n = ratemap ./ max(ratemap); % normalise cell hd
                hd3n = hd3n(:);            
                mx2 = xi(hd3n == max(hd3n)); % preferred angle (location of max frate)
                if length(mx2)>1
                    mx2 = mx2(1);
                elseif isempty(mx2)
                    mx2 = NaN;
                end
                hold on
                line(ax_pi.XLim,[mx2 mx2],'Color','r')

                % tuning width
                original_length = length(ratemap);
                [peaks,locs,widths,~] = findpeaks(r, 'WidthReference', 'halfprom');

                % filter to keep ONLY the peaks in the middle section
                middle_idx = (locs > original_length) & (locs <= 2 * original_length);
                middle_peaks = peaks(middle_idx);
                middle_widths = widths(middle_idx);

                % find the absolute maximum out of the middle peaks
                [~, max_idx] = max(middle_peaks);
                bin_size = median(diff(edg)); % bin size in degrees
                final_peak_width = middle_widths(max_idx) * bin_size;

                % directional info shuffle
                % observed values
                dmap = (d .* (1/pdata.pos_srate));
                smetric = get_spatial_info(dmap,ratemap,'metrics',{'spatial_info'});
                vals = [smetric.skaggs_si_bits_per_spike];
                        
                % shuffle values
                rng(999); % for reproducibility
                iti = 1000;
                vals_shuff = NaN(iti,1);
                if ~isempty(pspt) && ~all(isnan(pspt))
                    for ii = 1:iti
                        [~,sindx2] = simply_shuffle_spike_train(ppot,pspt,'spindx',sindx);
                        s_shuff = histcounts(ppit(sindx2),edg);
                        ratemap_shuff = s_shuff ./ dmap;
                        smetric = get_spatial_info(dmap,ratemap_shuff,'metrics',{'spatial_info'});
                        vals_shuff(ii) = smetric.skaggs_si_bits_per_spike;
                    end
                end
                vals_z = (vals - mean(vals_shuff,1,'omitnan')) ./ std(vals_shuff,1,'omitnan');
                vals_p = (sum(vals_shuff >= vals, 1) + 1) / (iti + 1);

                % accumulate
                dat.m2D_pitch(pp,1) = { xi };
                dat.m2D_pitch_dwellmap(pp,1) = { d };
                dat.m2D_pitch_ratemap(pp,1) = { ratemap };
                dat.m2D_pitch_pfd(pp,:) = mx2;    
                dat.m2D_pitch_pfd_frate(pp,:) = max(ratemap(:));  
                dat.m2D_pitch_avg_frate(pp,:) = mean(ratemap(:),'omitnan');             
                dat.m2D_pitch_info(pp,:) = vals;
                dat.m2D_pitch_info_zscored(pp,:) = vals_z;
                dat.m2D_pitch_info_p(pp,:) = vals_p;
                dat.m2D_pitch_tuning_width(pp,:) = final_peak_width;

                text(ax,0,-1.5,sprintf('Pitch (N shuff = %d)\nspatial info %.2f b/spike (z = %.2f, p = %.3f)\nTuning width: %.2f%c',iti,vals(1),vals_z(1),vals_p(1),final_peak_width,176),'units','normalized','FontSize',10);

            % 3D tuning curve
            % ax_3 = nexttile(pp*3);
            ax_3 = axes('Units','pixels','Position',[xvec(pp) yvec(3) 350 350]);  
                all_axes{5,pp} = gca;
                az_bins = az_map(1,:);
                tilt_bins = tilt_map(:,1);
                [AZ, TILT] = meshgrid(az_bins(:)', tilt_bins(:)');

                X_sph = cosd(TILT) .* cosd(AZ);
                Y_sph = cosd(TILT) .* sind(AZ);
                Z_sph = sind(TILT);
                
                if plot_3D_tuning_curves
                    [azimuth,elevation,~] = cart2sph(X_sph,Y_sph,Z_sph);
                    r = F3 ./ max(F3,[],"all",'omitmissing');
                    [X_sph,Y_sph,Z_sph] = sph2cart(azimuth,elevation,r);
                end

                surf(X_sph, Y_sph, Z_sph, F3, 'EdgeColor', 'none');
                axis equal off
                ax_3.CLim = [0,max([0.001 max(F3(:))],[],'omitmissing')];
                colormap(gca,colormap_rate)
                view(0,plot_el); % Nice isometric starting view

                % Add X, Y, and Z axis indicators
                hold on;
                plot3([-1.5 1.5], [0 0], [0 0], 'k-', 'LineWidth', 1.5); % X-axis
                text(1.65, 0, 0, 'X', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                
                plot3([0 0], [-1.5 1.5], [0 0], 'k-', 'LineWidth', 1.5); % Y-axis
                text(0, 1.65, 0, 'Y', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                
                plot3([0 0], [0 0], [-1.5 1.5], 'k-', 'LineWidth', 1.5); % Z-axis
                text(0, 0, 1.65, 'Z', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

                % accumulate
                dat.m3D_stats_azimuth_width(pp,1) = torus_stats.az_width_deg;
                dat.m3D_stats_azimuth_sigma(pp,1) = torus_stats.az_sigma_deg;
                dat.m3D_stats_tilt_width(pp,1) = torus_stats.tilt_width_deg;
                dat.m3D_stats_tilt_sigma(pp,1) = torus_stats.tilt_sigma_deg;

                text(ax_3,-0.1,0.85,sprintf('Azimuth:\nwidth %.2f%c\nsigma %.2f%c\nPitch:\nwidth %.2f%c\nsigma %.2f%c\n',torus_stats.az_width_deg,176,torus_stats.az_sigma_deg,176,torus_stats.tilt_width_deg,176,torus_stats.tilt_sigma_deg,176),'units','normalized','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','top');
                text(ax_3,-0.1,0.14,sprintf('Peak azimuth: %.f%c',center_angles(1,pp),176),'units','normalized','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom');
        end

        % make all the colormaps the same (zero to max of all maps)
        if normalise_all_colormaps
            cmax = 0;
            for ii = 1:size(all_axes,2)
                cmax = max([all_axes{2,ii}.CLim cmax],[],'all','omitmissing');
            end            
            for ii = 1:size(all_axes,2)
                if isempty(all_axes{2,ii})
                    continue
                end
                all_axes{2,ii}.CLim = [0 cmax];
            end

            cmax = 0;
            for ii = 1:size(all_axes,2)
                cmax = max([all_axes{5,ii}.CLim cmax],[],'all','omitmissing');
            end             
            for ii = 1:size(all_axes,2)
                if isempty(all_axes{5,ii})
                    continue
                end
                all_axes{5,ii}.CLim = [0 cmax];
            end
        end

        if any(center_all_2D_maps,"all")
            peak_ang = round( rad2deg( circ_mean( deg2rad(center_angles(center_all_2D_maps))' ) ) );
            for ii = 1:size(all_axes,2)
                all_axes{1,ii}.XLim = [peak_ang-180 peak_ang+180];                
                all_axes{2,ii}.XLim = [peak_ang-180 peak_ang+180];
                all_axes{3,ii}.XLim = [peak_ang-180 peak_ang+180];
            end   
        end

        if any(center_all_3D_maps,"all")
            if centering_method==1 || centering_method==3
                peak_ang = round( rad2deg( circ_mean( deg2rad(center_angles(center_all_3D_maps))' ) ) );
            elseif centering_method==2 || centering_method==4
                peak_ang = round( rad2deg( circ_median( deg2rad(center_angles(center_all_3D_maps))' ) ) );
            end

            % MATLAB's camera azimuth is offset by 90 degrees from standard cartesian
            % so we add another 90 degrees no matter what
            % lastly, the user has the option to offset the view angle by a
            % consistent offset_3d_tuning_curves value
            cam_az = peak_ang + 90 + offset_3d_tuning_curves; 
            for ii = 1:size(all_axes,2)
                all_axes{5,ii}.View = [cam_az,plot_el];    
                text(all_axes{5,ii},-0.1,0.08,sprintf('View azimuth: %.f%c',cam_az-90,176),'units','normalized','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','bottom');
            end  
        end

% set(gcf,'visible','on')  
% keyboard
% return
%% >>>>>>>>>> Test stability between the arenas
        % data for arena 1
        partn = 1;
        part_now = part_config.part_names{partn};        
        ppot1 = double( pdata.(part_now).pot ); % pos t 
        ppit1 = double( pdata.(part_now).pit ); % pitch for this part, in rads   
        pyaw1 = double( pdata.(part_now).yaw ); % pitch for this part, in rads            
        idx1 = find( ismember(sdata.uci,uci) & sdata.partn==partn );
        sindx1 = sdata.spt_pot_index{idx1}; 
        pspt1 = ppot1(sindx1);
        d1y = dwell_maps_2d{1,partn}; % yaw dwellmap
        d1p = dwell_maps_2d{2,partn}; % pitch dwellmap

        % data for hills
        partn = 2;
        part_now = part_config.part_names{partn};        
        ppot2 = double( pdata.(part_now).pot ); % pos t 
        ppit2 = double( pdata.(part_now).pit ); % pitch for this part, in rads   
        pyaw2 = double( pdata.(part_now).yaw ); % pitch for this part, in rads            
        idx2 = find( ismember(sdata.uci,uci) & sdata.partn==partn );
        sindx2 = sdata.spt_pot_index{idx2}; 
        pspt2 = ppot2(sindx2);
        d2y = dwell_maps_2d{1,partn}; % yaw dwellmap
        d2p = dwell_maps_2d{2,partn}; % pitch dwellmap

        % data for arena 2
        partn = 3;
        part_now = part_config.part_names{partn};        
        ppot3 = double( pdata.(part_now).pot ); % pos t 
        ppit3 = double( pdata.(part_now).pit ); % pitch for this part, in rads   
        pyaw3 = double( pdata.(part_now).yaw ); % pitch for this part, in rads            
        idx3 = find( ismember(sdata.uci,uci) & sdata.partn==partn );
        sindx3 = sdata.spt_pot_index{idx3}; 
        pspt3 = ppot3(sindx3);
        d3y = dwell_maps_2d{1,partn}; % yaw dwellmap
        d3p = dwell_maps_2d{2,partn}; % pitch dwellmap

        % observed arena 1 vs arena 2 value
        dat.arena1_vs_arena2_azimuth_corr = NaN(nparts+1,1);
        dat.arena1_vs_arena2_pitch_corr = NaN(nparts+1,1);
        dat.arena1_vs_arena2_map_corr = NaN(nparts+1,1);
        dat.arena1_vs_hills_azimuth_corr = NaN(nparts+1,1);
        dat.arena1_vs_hills_pitch_corr = NaN(nparts+1,1);  
        dat.arena1_vs_hills_map_corr = NaN(nparts+1,1);  

        m1 = reshape(dat.m2D_azimuth_ratemap{1,1},[],1);
        m2 = reshape(dat.m2D_azimuth_ratemap{3,1},[],1);
        if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
            dat.arena1_vs_arena2_azimuth_corr(1,1) = corr(m1,m2,'type','Pearson','rows','pairwise');
        end

        m1 = reshape(dat.m2D_pitch_ratemap{1,1},[],1);
        m2 = reshape(dat.m2D_pitch_ratemap{3,1},[],1);
        if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
            dat.arena1_vs_arena2_pitch_corr(1,1) = corr(m1,m2,'type','Pearson','rows','pairwise');
        end   

        m1 = reshape(dat.m3D_HD_ratemap{1,1},[],1);
        m2 = reshape(dat.m3D_HD_ratemap{3,1},[],1);
        if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
            dat.arena1_vs_arena2_map_corr(1,1) = corr(m1,m2,'type','Pearson','rows','pairwise');
        end  

        m1 = reshape(dat.m2D_azimuth_ratemap{1,1},[],1);
        m2 = reshape(dat.m2D_azimuth_ratemap{2,1},[],1);
        if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
            dat.arena1_vs_hills_azimuth_corr(1,1) = corr(m1,m2,'type','Pearson','rows','pairwise');        
        end   

        m1 = reshape(dat.m2D_pitch_ratemap{1,1},[],1);
        m2 = reshape(dat.m2D_pitch_ratemap{2,1},[],1);
        if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
            dat.arena1_vs_hills_pitch_corr(1,1) = corr(m1,m2,'type','Pearson','rows','pairwise');
        end

        m1 = reshape(dat.m3D_HD_ratemap{1,1},[],1);
        m2 = reshape(dat.m3D_HD_ratemap{2,1},[],1);
        if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
            dat.arena1_vs_hills_map_corr(1,1) = corr(m1,m2,'type','Pearson','rows','pairwise');
        end  

        % shuffle values
        rng(999); % for reproducibility
        iti = 1000;
        r_shuff = NaN(iti,4);
        for ii = 1:iti
            % shuffled arena 1 map yaw
            ratemap_shuff1 = NaN;
            if ~isempty(pspt1)
                [~,shindx1] = simply_shuffle_spike_train(ppot1,pspt1,'spindx',sindx1);
                s_shuff1 = histcounts(pyaw1(shindx1),edg);
                ratemap_shuff1 = s_shuff1 ./ (d1y .* (1/pdata.pos_srate));
            end

            % shuffled hills map yaw
            ratemap_shuff2 = NaN;
            if ~isempty(pspt2)
                [~,shindx2] = simply_shuffle_spike_train(ppot2,pspt2,'spindx',sindx2);
                s_shuff2 = histcounts(pyaw2(shindx2),edg);
                ratemap_shuff2 = s_shuff2 ./ (d2y .* (1/pdata.pos_srate));
            end

            % shuffled arena 2 map yaw
            ratemap_shuff3 = NaN;
            if ~isempty(pspt3)
                [~,shindx3] = simply_shuffle_spike_train(ppot3,pspt3,'spindx',sindx3);
                s_shuff3 = histcounts(pyaw3(shindx3),edg);
                ratemap_shuff3 = s_shuff3 ./ (d3y .* (1/pdata.pos_srate));
            end

            m1 = ratemap_shuff1(:);
            m2 = ratemap_shuff3(:);
            if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
                r_shuff(ii,1) = corr(m1,m2,'type','Pearson','rows','pairwise');% arena1 vs arena2 yaw
            end
            m2 = ratemap_shuff2(:);
            if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
                r_shuff(ii,2) = corr(m1,m2,'type','Pearson','rows','pairwise');% arena1 vs hills yaw
            end

            % shuffled arena 1 map pitch
            ratemap_shuff1 = NaN;
            if ~isempty(pspt1)     
                [~,shindx1] = simply_shuffle_spike_train(ppot1,pspt1,'spindx',sindx1);
                s_shuff1 = histcounts(ppit1(shindx1),edg);
                ratemap_shuff1 = s_shuff1 ./ (d1p .* (1/pdata.pos_srate));
            end

            % shuffled hills map pitch
            ratemap_shuff2 = NaN;
            if ~isempty(pspt2)     
                [~,shindx2] = simply_shuffle_spike_train(ppot2,pspt2,'spindx',sindx2);
                s_shuff2 = histcounts(ppit2(shindx2),edg);
                ratemap_shuff2 = s_shuff2 ./ (d2p .* (1/pdata.pos_srate));
            end

            % shuffled arena 2 map pitch
            ratemap_shuff3 = NaN;
            if ~isempty(pspt3)         
                [~,shindx3] = simply_shuffle_spike_train(ppot3,pspt3,'spindx',sindx3);
                s_shuff3 = histcounts(ppit3(shindx3),edg);
                ratemap_shuff3 = s_shuff3 ./ (d3p .* (1/pdata.pos_srate));
            end

            m1 = ratemap_shuff1(:);
            m2 = ratemap_shuff3(:);
            if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
                r_shuff(ii,3) = corr(m1,m2,'type','Pearson','rows','pairwise');% arena1 vs arena2 pitch
            end
            m2 = ratemap_shuff2(:);
            if ~isempty(m1) && ~isempty(m2) && ~all(isnan(m1)) && ~all(isnan(m2))
                r_shuff(ii,4) = corr(m1,m2,'type','Pearson','rows','pairwise'); % arena1 vs hills pitch
            end
        end

        dat.arena1_vs_arena2_azimuth_corr_zscored = NaN(nparts+1,1);
        dat.arena1_vs_arena2_pitch_corr_zscored = NaN(nparts+1,1);
        dat.arena1_vs_hills_azimuth_corr_zscored = NaN(nparts+1,1);
        dat.arena1_vs_hills_pitch_corr_zscored = NaN(nparts+1,1);
        dat.arena1_vs_arena2_azimuth_corr_zscored(1,1) = (dat.arena1_vs_arena2_azimuth_corr(1,1) - mean(r_shuff(:,1),1,'omitnan')) ./ std(r_shuff(:,1),1,'omitnan');
        dat.arena1_vs_arena2_pitch_corr_zscored(1,1) = (dat.arena1_vs_arena2_pitch_corr(1,1) - mean(r_shuff(:,3),1,'omitnan')) ./ std(r_shuff(:,3),1,'omitnan');
        dat.arena1_vs_hills_azimuth_corr_zscored(1,1) = (dat.arena1_vs_hills_azimuth_corr(1,1) - mean(r_shuff(:,2),1,'omitnan')) ./ std(r_shuff(:,2),1,'omitnan');
        dat.arena1_vs_hills_pitch_corr_zscored(1,1) = (dat.arena1_vs_hills_pitch_corr(1,1) - mean(r_shuff(:,4),1,'omitnan')) ./ std(r_shuff(:,4),1,'omitnan');

        dat.arena1_vs_arena2_azimuth_corr_pvalue = NaN(nparts+1,1);
        dat.arena1_vs_arena2_pitch_corr_pvalue = NaN(nparts+1,1);
        dat.arena1_vs_hills_azimuth_corr_pvalue = NaN(nparts+1,1);
        dat.arena1_vs_hills_pitch_corr_pvalue = NaN(nparts+1,1);
        pfun = @(val,shuffs,iti) (sum(shuffs >= val, 1) + 1) / (iti + 1);
        dat.arena1_vs_arena2_azimuth_corr_pvalue(1,1) = pfun(dat.arena1_vs_arena2_azimuth_corr(1,1),r_shuff(:,1),iti);
        dat.arena1_vs_arena2_pitch_corr_pvalue(1,1) = pfun(dat.arena1_vs_arena2_pitch_corr(1,1),r_shuff(:,3),iti);
        dat.arena1_vs_hills_azimuth_corr_pvalue(1,1) = pfun(dat.arena1_vs_hills_azimuth_corr(1,1),r_shuff(:,2),iti);
        dat.arena1_vs_hills_pitch_corr_pvalue(1,1) = pfun(dat.arena1_vs_hills_pitch_corr(1,1),r_shuff(:,4),iti);

        ann_str = sprintf('N shuff = %d | arena 1 vs arena 2 azimuth r = %.2f (z = %.2f, p = %.5f) | arena 1 vs hills azimuth r = %.2f (z = %.2f, p = %.5f) | arena 1 vs arena 2 pitch r = %.2f (z = %.2f, p = %.5f) | arena 1 vs hills pitch r = %.2f (z = %.2f, p = %.5f)',iti...
            ,dat.arena1_vs_arena2_azimuth_corr(1,1),dat.arena1_vs_arena2_azimuth_corr_zscored(1,1),dat.arena1_vs_arena2_azimuth_corr_pvalue(1,1)...
            ,dat.arena1_vs_hills_azimuth_corr(1,1),dat.arena1_vs_hills_azimuth_corr_zscored(1,1),dat.arena1_vs_hills_azimuth_corr_pvalue(1,1)...
            ,dat.arena1_vs_arena2_pitch_corr(1,1),dat.arena1_vs_arena2_pitch_corr_zscored(1,1),dat.arena1_vs_arena2_pitch_corr_pvalue(1,1)...
            ,dat.arena1_vs_hills_pitch_corr(1,1),dat.arena1_vs_hills_pitch_corr_zscored(1,1),dat.arena1_vs_hills_pitch_corr_pvalue(1,1));
        annotation('textbox',[0.05, 0.018, 1, 0],'string',ann_str,'FontSize',10,'LineStyle','none','interpreter','none','VerticalAlignment','bottom','HorizontalAlignment','left');      

        ann_str = sprintf('arena 1 vs arena 2 map r = %.2f | arena 1 vs hills map r = %.2f '...
            ,dat.arena1_vs_arena2_map_corr(1,1),dat.arena1_vs_hills_map_corr(1,1));
        annotation('textbox',[0.05, 0, 1, 0],'string',ann_str,'FontSize',10,'LineStyle','none','interpreter','none','VerticalAlignment','bottom','HorizontalAlignment','left');      


% set(gcf,'visible','on')  
% keyboard
% return        
%% >>>>>>>>>> Save the overall figure
        % Save the figure  
        [~,~,~] = mkdir([pwd '\' pdata.outname '\part_figures']); % create a folder to hold outputs 
        fname = [pwd '\' pdata.outname '\part_figures\' uci '_pitch.png'];
        if fast_figures
            frame = getframe(fig_clust); % fig is the figure handle to save
            [raster, raster_map] = frame2im(frame); % raster is the rasterized image, raster_map is the colormap
            if isempty(raster_map)
                imwrite(raster, fname);
            else
                imwrite(raster, raster_map, fname); % fig_file is the path to the image
            end        
        else
            exportgraphics(fig_clust,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',300);  
        end
        close(fig_clust);  

        % save the temporary dat table
        all_dat = [all_dat; dat];

    end
    save([pwd '\' pdata.outname '\part_figures\pitch_tuning_analysis.mat'],'all_dat')

end











