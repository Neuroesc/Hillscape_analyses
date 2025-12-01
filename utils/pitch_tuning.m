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
    disp(sprintf('\t...done'))    
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Session data and figure
    % open figure
    fig_clust = figure('Units','pixels','Position',[100 50 1200 800],'visible','off');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white 

    ann_str = sprintf('Rat: %s, Date: %s, Analysed: %s',sdata.rat{1},sdata.date{1},datestr(now,'yyyy-mm-dd-HH-MM-SS'));
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',10,'LineStyle','none','interpreter','none'); 

    map_now = 'turbo';

    % run through every part
    if 1
        dwell_maps_3d = cell(1,nparts);
        for pp = 1:nparts % for every part               
            part_now = part_config.part_names{pp};
            if pp==1
                disp(sprintf('\t\t%s',part_now))      
            else
                disp(sprintf('\b, %s',part_now))      
            end
               
            % get position data
            ppox = double( pdata.(part_now).pox ).*10; % pos x for this part, in mm
            ppoy = double( pdata.(part_now).poy ).*10; % pos y for this part, in mm
            ppoz = double( pdata.(part_now).poz ).*10; % pos z for this part, in mm          
            ppit = double( pdata.(part_now).pit ); % pitch for this part, in rads   
            pyaw = double( pdata.(part_now).yaw ); % pitch for this part, in rads
            % pos = [ppox ppoy ppoz ppit pyaw];
    
            % plot trajectory
            ix = sub2ind([nparts 3],pp,1); % index for this subplot
            ax = subplot(3,nparts,ix);
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
    
            % plot dwell map
            ix = sub2ind([nparts 3],pp,2); % index for this subplot
            ax = subplot(3,nparts,ix,'Units','pixels');
                vadjust = 0;    
                if isempty(dwell_maps_3d{1,pp})
                    npoints = numel(pyaw); % number of points
                    [X,Y,Z] = sph2cart(deg2rad(pyaw),deg2rad(ppit+vadjust),ones(npoints,1)); % convert to XYZ
                    [Xd,Yd,Zd,F1] = projectEIGS([X,Y,Z],64,10); % project these onto sphere     
                    F1 = (F1 ./ sum(F1(:))) .* npoints .* (1/50);
                    dwell_maps_3d{1,pp} = F1;
                end   
                [TH,PHI,R] = cart2sph(Xd,Yd,Zd);
    
                % surf('XData',TH,'YData',PHI,'ZData',ones(size(PHI)),'CData',F1,'EdgeColor','none'); hold on;
                F1 = dwell_maps_3d{1,pp};
                F1(F1<min_dwell) = NaN;
                imagesc([-180 180],[-90 90],F1,'alphadata',~isnan(F1))
                xlabel(sprintf('Azimuth (%c)',176))
                ylabel(sprintf('Pitch (%c)',176))
                axis xy tight
                daspect([1 1 1])
                view(0,90);
                set(gca,'FontSize',12)
                ax.CLim = [0,max(F1(:))];
                colormap(gca,map_now)
                ax.XTick = -180:90:180;
                ax.YTick = -90:45:90;
    
                if pp==nparts
                    axc = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+20 ax.Position(2)+10 12 ax.Position(4)-10]);
                        mat = (linspace(0,max(F1(:)),100))';
                        imagesc([1 1],[min(mat(:)) max(mat(:))],mat);
                        colormap(axc,map_now);
                        axis xy
        
                        axc.YTick = [];
                        axc.XTick = [];
                        axc.YAxisLocation = 'right';
                        text(0.5,1.25,sprintf('Dwell\ntime (s)'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,1.1,sprintf('Max'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,-0.1,sprintf('0 Hz'),'FontSize',8,'HorizontalAl','center','Units','normalized') 
                end
    
    % keyboard
        end
    end

%% >>>>>>>>>> Save the overall figure
    % Save the figure  
    [~,~,~] = mkdir([pwd '\' pdata.outname '\part_figures']); % create a folder to hold outputs 
    fname = [pwd '\' pdata.outname '\part_figures\pitch_dwell.png'];
    if fast_figures
        frame = getframe(fig_clust); % fig is the figure handle to save
        [raster, raster_map] = frame2im(frame); % raster is the rasterized image, raster_map is the colormap
        if isempty(raster_map)
            imwrite(raster, fname);
        else
            imwrite(raster, raster_map, fname); % fig_file is the path to the image
        end            
    else
        exportgraphics(fig_clust,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',150);  
    end
    close(fig_clust);    


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
    
    dwell_maps_2d = cell(2,nparts);
    for uu = 1:length(ucis)     
        dat = table;

        % open figure
        fig_clust = figure('Units','pixels','Position',[100 100 1400 800],'visible','off');
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
        yvec = [615 370];
        xnow = 150;

        xt = 1:nparts-1;
        xvec = [xnow xnow+(xsiz*xt)+(xbuff*xt)];
        set(fig_clust,'Position',[100 100 max(xvec)+xsiz+50 800])

        % run through every part
        for pp = 1:nparts % for every part       
            part_now = part_config.part_names{pp};
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

            % get position data
            ppox = double( pdata.(part_now).pox ).*10; % pos x for this part, in mm
            ppoy = double( pdata.(part_now).poy ).*10; % pos y for this part, in mm
            ppoz = double( pdata.(part_now).poz ).*10; % pos z for this part, in mm  
            ppot = double( pdata.(part_now).pot ); % pos t 
            ppit = double( pdata.(part_now).pit ); % pitch for this part, in rads   
            pyaw = double( pdata.(part_now).yaw ); % pitch for this part, in rads
            % pos = [ppox ppoy ppoz ppit pyaw];

            % get spike data
            idx = find( ismember(sdata.uci,uci) & sdata.partn==pp );
            dat.partn(pp,1) = sdata.partn(idx);          
            dat.nspikes(pp,1) = sdata.nspikes(idx);
            dat.frate(pp,1) = sdata.frate(idx);

            sindx = sdata.spt_pot_index{idx};            
            pspx = ppox(sindx);
            pspy = ppoy(sindx);
            pspz = ppoz(sindx);
            pspt = ppot(sindx);
            psyaw = pyaw(sindx);
            pspit = ppit(sindx);

            % % plot trajectory
            % ix = sub2ind([nparts 4],pp,1); % index for this subplot
            % subplot(4,nparts,ix)
            %     plot3(ppox,ppoy,ppoz,'k'); hold on;
            %     plot3(pspx,pspy,pspz,'r.','MarkerSize',10);
            %     daspect([1 1 1])
            %     axis xy off tight
            %     view(0,90);
            %     set(gca,'FontSize',12)

            % plot spikes and positions
            ax = axes('Units','pixels','Position',[xvec(pp) yvec(1) xsiz ysiz]);
                if pp==1
                    ax_1 = ax;
                end

                if isempty(pspt)
                    continue
                end
                
                dindax = abs([0; diff(pyaw)])>deg2rad(90);
                yaw_plot = pyaw;
                yaw_plot(dindax,:) = NaN;
    
                plot(yaw_plot,ppit,'k'); hold on;
                plot(psyaw,pspit,'r.','MarkerSize',10); hold on;                
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

            % plot 3D tuning curve (flat projection)
            ax = axes('Units','pixels','Position',[xvec(pp) yvec(2) xsiz ysiz]);
                F1 = dwell_maps_3d{1,pp};

                npoints = numel(sindx); % number of points                
                [X,Y,Z] = sph2cart(deg2rad(pyaw(sindx)),deg2rad(ppit(sindx)+vadjust),ones(npoints,1)); % convert to XYZ
                [~,~,~,F2] = projectEIGS([X,Y,Z],64,10); % project these onto sphere 
                F2 = (F2 ./ sum(F2(:))) .* npoints;
                F3 = F2 ./ F1;
                F3(F1<min_dwell) = NaN;                
            
                imagesc([-180 180],[-90 90],F3,'alphadata',~isnan(F3))
                xlabel(sprintf('Azimuth (%c)',176))
                ylabel(sprintf('Pitch (%c)',176))
                axis xy tight
                daspect([1 1 1])
                view(0,90);
                set(gca,'FontSize',12)
                ax.CLim = [0,max([0.001 max(F3(:))],[],'omitmissing')];
                colormap(gca,map_now)
                ax.XTick = [];
                ax.YTick = [];

                dat.m3D_HD_dwellmap(pp,1) = { F1 };
                dat.m3D_HD_ratemap(pp,1) = { F3 };

            % azimuthal tuning curve
            ax_az = axes('Units','pixels','Position',[ax.Position(1) ax.Position(2)-60 ax.Position(3) 55]);
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
                ratemap_az = interp1(xi,ratemap,edg,'linear',NaN);
                
                bar(edg,ratemap_az,1,'k')
                xlabel(sprintf('Azimuth (%c)',176))
                ylabel('Firing Rate (Hz)')  
                set(gca,'FontSize',12)
                box off
                ax_az.YAxisLocation = 'right';
                ax_az.XTick = [-180:90:180];
                ax_az.FontSize = 10;

                % head direction analyses
                hd3n = ratemap_az ./ max(ratemap_az); % normalise cell hd
                hd3n = hd3n(:);            
                mx1 = edg(hd3n == max(hd3n)); % preferred angle (location of max frate)
                if length(mx1)>1
                    mx1 = mx1(1);
                elseif isempty(mx1)
                    mx1 = NaN;
                end
                % text(0,1.05,sprintf('r: %.2f, PFD%c: %.2f',r1,176,mx1),'Units','normalized','FontSize',10,'HorizontalAlignment','left')
                hold on
                line([mx1 mx1],ax_az.YLim,'Color','r')

                % directional info shuffle
                % observed values
                smetric = spatialMETRICS(ratemap,d);
                vals = [smetric.kl_divergence smetric.spatial_information smetric.signal_to_noise];
            
                % shuffle values
                rng(999); % for reproducibility
                iti = 100;
                vals_shuff = NaN(iti,3);
                for ii = 1:iti
                    [~,sindx2] = shift_spike_train(ppot,pspt,'spindx',sindx);
                    s_shuff = histcounts(pyaw(sindx2),edg);
                    ratemap_shuff = s_shuff ./ (d .* (1/pdata.pos_srate));
                    smetric = spatialMETRICS(ratemap_shuff,d);
                    vals_shuff(ii,:) = [smetric.kl_divergence smetric.spatial_information smetric.signal_to_noise];
                end
                % vals_z = (vals - mean(vals_shuff,1,'omitnan')) ./ std(vals_shuff,1,'omitnan');
                [vals_z,vals_p] = stats_z_probability(vals,vals_shuff);
keyboard
                % accumulate
                dat.m2D_azimuth_dwellmap(pp,1) = { interp1(xi,d,edg,'linear',NaN) };
                dat.m2D_azimuth_ratemap(pp,1) = { ratemap_az };
                dat.m2D_azimuth_pfd(pp,:) = mx1;  
                dat.m2D_azimuth_pfd_frate(pp,:) = max(ratemap_az(:));  
                dat.m2D_azimuth_avg_frate(pp,:) = mean(ratemap_az(:),'omitnan');                  
                dat.m2D_azimuth_info(pp,:) = vals;
                dat.m2D_azimuth_info_zscored(pp,:) = vals_z;
                dat.m2D_azimuth_info_p(pp,:) = vals_p;

                text(ax,0,-1.3,sprintf('Azimuth\nkl divergence %.2f (z = %.2f, p = %.3f)\nspatial info %.2f (z = %.2f, p = %.3f)\nSNR %.2f (z = %.2f, p = %.3f)\n',vals(1),vals_z(1),vals_p(1),vals(2),vals_z(2),vals_p(2),vals(3),vals_z(3),vals_p(3)),'units','normalized','FontSize',10);

            % pitch tuning curve
            ax_pi = axes('Units','pixels','Position',[ax.Position(1)-60 ax.Position(2) 55 ax.Position(4)]);
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
                ratemap_pit = interp1(xi,ratemap,edg,'linear',NaN);
                
                barh(edg,ratemap_pit,1,'k')
                ylabel(sprintf('Pitch (%c)',176))
                xlabel('Firing Rate (Hz)')  
                set(gca,'FontSize',12)
                box off
                ax_pi.XAxisLocation = 'top';
                ax_pi.YTick = [-90:45:90];
                ax_pi.FontSize = 10;

                % head direction analyses
                hd3n = ratemap_pit ./ max(ratemap_pit); % normalise cell hd
                hd3n = hd3n(:);            
                mx2 = edg(hd3n == max(hd3n)); % preferred angle (location of max frate)
                if length(mx2)>1
                    mx2 = mx2(1);
                elseif isempty(mx2)
                    mx2 = NaN;
                end
                % text(0,1.05,sprintf('r: %.2f, PFD%c: %.2f',r2,176,mx2),'Units','normalized','FontSize',10,'HorizontalAlignment','left')
                hold on
                line(ax_pi.XLim,[mx2 mx2],'Color','r')

                % directional info shuffle
                % observed values
                smetric = spatialMETRICS(ratemap,d);
                vals = [smetric.kl_divergence smetric.spatial_information smetric.signal_to_noise];
            
                % shuffle values
                rng(999); % for reproducibility
                iti = 100;
                vals_shuff = NaN(iti,3);
                if ~isempty(pspt) && ~all(isnan(pspt))
                    for ii = 1:iti
                        [~,sindx2] = shift_spike_train(ppot,pspt,'spindx',sindx);
                        s_shuff = histcounts(ppit(sindx2),edg);
                        ratemap_shuff = s_shuff ./ (d .* (1/pdata.pos_srate));
                        smetric = spatialMETRICS(ratemap_shuff,d);
                        vals_shuff(ii,:) = [smetric.kl_divergence smetric.spatial_information smetric.signal_to_noise];
                    end
                end
                % vals_z = (vals - mean(vals_shuff,1,'omitnan')) ./ std(vals_shuff,1,'omitnan');
                [vals_z,vals_p] = stats_z_probability(vals,vals_shuff);

                % accumulate
                dat.m2D_pitch_dwellmap(pp,1) = { interp1(xi,d,edg,'linear',NaN) };
                dat.m2D_pitch_ratemap(pp,1) = { ratemap_pit };
                dat.m2D_pitch_pfd(pp,:) = mx2;    
                dat.m2D_pitch_pfd_frate(pp,:) = max(ratemap_pit(:));  
                dat.m2D_pitch_avg_frate(pp,:) = mean(ratemap_pit(:),'omitnan');             
                dat.m2D_pitch_info(pp,:) = vals;
                dat.m2D_pitch_info_zscored(pp,:) = vals_z;
                dat.m2D_pitch_info_p(pp,:) = vals_p;

                text(ax,0,-1.9,sprintf('Pitch\nkl divergence %.2f (z = %.2f, p = %.3f)\nspatial info %.2f (z = %.2f, p = %.3f)\nSNR %.2f (z = %.2f, p = %.3f)\n',vals(1),vals_z(1),vals_p(1),vals(2),vals_z(2),vals_p(2),vals(3),vals_z(3),vals_p(3)),'units','normalized','FontSize',10);
                if pp==1
                    ax_1 = ax;
                end

                % %% add to sdata
                % sdata.hd_map_3d(idx) = { single(F3) }; % azimuth x pitch tuning curve
                % sdata.hd_3d_curves(idx) = { single([ratemap_az(:)'; ratemap_pit(:)']) }; % azimuth tuning curve; pitch tuning curve
                % sdata.hd_3d_info(idx,:) = single([r1 r2 mx1 mx2]); % azimuth rayleigh, pitch rayleigh, azimuth pfd, pitch pfd
        end

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
        dat.arena1_vs_arena2_azimuth_corr = NaN(nparts,1);
        dat.arena1_vs_arena2_pitch_corr = NaN(nparts,1);
        dat.arena1_vs_hills_azimuth_corr = NaN(nparts,1);
        dat.arena1_vs_hills_pitch_corr = NaN(nparts,1);  

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

        % shuffle values
        rng(999); % for reproducibility
        iti = 100;
        r_shuff = NaN(iti,4);
        for ii = 1:iti
            % shuffled arena 1 map yaw
            ratemap_shuff1 = NaN;
            if ~isempty(pspt1)
                [~,shindx1] = shift_spike_train(ppot1,pspt1,'spindx',sindx1);
                s_shuff1 = histcounts(pyaw1(shindx1),edg);
                ratemap_shuff1 = s_shuff1 ./ (d1y .* (1/pdata.pos_srate));
            end

            % shuffled hills map yaw
            ratemap_shuff2 = NaN;
            if ~isempty(pspt2)
                [~,shindx2] = shift_spike_train(ppot2,pspt2,'spindx',sindx2);
                s_shuff2 = histcounts(pyaw2(shindx2),edg);
                ratemap_shuff2 = s_shuff2 ./ (d2y .* (1/pdata.pos_srate));
            end

            % shuffled arena 2 map yaw
            ratemap_shuff3 = NaN;
            if ~isempty(pspt3)
                [~,shindx3] = shift_spike_train(ppot3,pspt3,'spindx',sindx3);
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
                [~,shindx1] = shift_spike_train(ppot1,pspt1,'spindx',sindx1);
                s_shuff1 = histcounts(ppit1(shindx1),edg);
                ratemap_shuff1 = s_shuff1 ./ (d1p .* (1/pdata.pos_srate));
            end

            % shuffled hills map pitch
            ratemap_shuff2 = NaN;
            if ~isempty(pspt2)                        
                [~,shindx2] = shift_spike_train(ppot2,pspt2,'spindx',sindx2);
                s_shuff2 = histcounts(ppit2(shindx2),edg);
                ratemap_shuff2 = s_shuff2 ./ (d2p .* (1/pdata.pos_srate));
            end

            % shuffled arena 2 map pitch
            ratemap_shuff3 = NaN;
            if ~isempty(pspt3)             
                [~,shindx3] = shift_spike_train(ppot3,pspt3,'spindx',sindx3);
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

        dat.arena1_vs_arena2_azimuth_corr_zscored = NaN(nparts,1);
        dat.arena1_vs_arena2_pitch_corr_zscored = NaN(nparts,1);
        dat.arena1_vs_hills_azimuth_corr_zscored = NaN(nparts,1);
        dat.arena1_vs_hills_pitch_corr_zscored = NaN(nparts,1);
        dat.arena1_vs_arena2_azimuth_corr_zscored(1,1) = (dat.arena1_vs_arena2_azimuth_corr(1,1) - mean(r_shuff(:,1),1,'omitnan')) ./ std(r_shuff(:,1),1,'omitnan');
        dat.arena1_vs_arena2_pitch_corr_zscored(1,1) = (dat.arena1_vs_arena2_pitch_corr(1,1) - mean(r_shuff(:,3),1,'omitnan')) ./ std(r_shuff(:,3),1,'omitnan');
        dat.arena1_vs_hills_azimuth_corr_zscored(1,1) = (dat.arena1_vs_hills_azimuth_corr(1,1) - mean(r_shuff(:,2),1,'omitnan')) ./ std(r_shuff(:,2),1,'omitnan');
        dat.arena1_vs_hills_pitch_corr_zscored(1,1) = (dat.arena1_vs_hills_pitch_corr(1,1) - mean(r_shuff(:,4),1,'omitnan')) ./ std(r_shuff(:,4),1,'omitnan');

        text(ax_1,0,-2.5,sprintf('arena 1 vs arena 2 azimuth r %.2f (z = %.2f)\narena 1 vs hills azimuth r %.2f (z = %.2f)\n arena 1 vs arena 2 pitch r %.2f (z = %.2f)\narena 1 vs hills pitch r %.2f (z = %.2f)'...
            ,dat.arena1_vs_arena2_azimuth_corr(1,1),dat.arena1_vs_arena2_azimuth_corr_zscored(1,1),dat.arena1_vs_hills_azimuth_corr(1,1),dat.arena1_vs_hills_azimuth_corr_zscored(1,1)...
            ,dat.arena1_vs_arena2_pitch_corr(1,1),dat.arena1_vs_arena2_pitch_corr_zscored(1,1),dat.arena1_vs_hills_pitch_corr(1,1),dat.arena1_vs_hills_pitch_corr_zscored(1,1)),'units','normalized','FontSize',10);

%% >>>>>>>>>> Save the overall figure
% set(gcf,'visible','on')
% keyboard
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
            exportgraphics(fig_clust,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',150);  
        end
        close(fig_clust);  

        % save the temporary dat table
        save([pwd '\' pdata.outname '\part_figures\' uci '_pitch.mat'],'dat')

    end

% %% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
%     save([pwd '\' pdata.outname '\sdata.mat'],'sdata'); % save session data
%     analysis_log({'pitch_modulation_analysis'},1,'version',{'v1.0.0'});
% 
%     % finish up
%     toc1 = toc/60;
%     disp(sprintf('%s has finished. It took %0.3f seconds or %0.3f minutes',stk.name,toc,toc1)) % Stop counting time and display results
%     disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
%     disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',[pwd '\' pdata.outname],' &'');">','klustest folder','</a>'])
%     disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
% 











































