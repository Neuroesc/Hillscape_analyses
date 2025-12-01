function PIT_pitch_tuning(data_dir,rnow,dnow)
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
    skipfigs = 1; % 1 = skip making a figure if it exists already    
    save_figs = 1; % 1 = make and save figures (VERY time consuming, takes about 90% of klustest's time)
    fast_figures = 1; % 1 = 4-5x faster figure saving, but at a slightly lower quality
    
    % settings 
    min_dwell = 0.05;
    mapset.drive_height_mm = 20;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PREPARE DATA
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Tetrodes and sessions
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    stk = dbstack;
    tnow = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',stk.name,tnow))
    
%% >>>>>>>>>> Prepare the data
    sname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    disp(sprintf('Loading cluma: %s',sname))    
    load(sname,'cluma'); % load saved session data
    pdata = cluma.Properties.CustomProperties.pdata;
    part_names = {'arena1','hills','arena2'};
    nparts = size(pdata.session_times,1);
    disp(sprintf('\t...%d sessions',nparts));   
    disp(sprintf('\t...recording date: %s',dnow));            
    disp(sprintf('\t...done')) 

%% >>>>>>>>>> Analyse trajectory and add 3D HD
    pos = pdata.pos;
    rpoxn = pos.rx;
    rpoyn = -pos.ry;
    rpozn = -pos.rz;
    bpoxn = pos.bx;
    bpoyn = -pos.by;
    bpozn = -pos.bz;
    gpoxn = pos.gx;
    gpoyn = -pos.gy;
    gpozn = -pos.gz;

    % calculate roll pitch yaw
    disp(sprintf('Loading 3D trajectory...'))                                                             
    P0 = [rpoxn(:) rpoyn(:) rpozn(:)]; % red position
    P1 = [bpoxn(:) bpoyn(:) bpozn(:)]; % blue position
    P2 = [gpoxn(:) gpoyn(:) gpozn(:)]; % green position  
    mpos = [mean([rpoxn(:) gpoxn(:) bpoxn(:)],2,'omitnan') mean([rpoyn(:) gpoyn(:) bpoyn(:)],2,'omitnan') mean([rpozn(:) gpozn(:) bpozn(:)],2,'omitnan')]; % mean position
    % the red, green and blue LEDs form a plane, this next line will find the surface normal of this
    % i.e. the vector which points directly 'up' perpendicular from the LED array
    % we can use this to correct the position data and get the animal's head position
    % note that if the surface normal z coordinate is below the medium point the animal must be inverted
    led_normal = normalize( cross(P0-P1, P0-P2),2,'norm' ); 
    
    % correct the mean position point so it is at the animal's head
    hpos = mpos - (led_normal*mapset.drive_height_mm); % hpos = head position
    
    % the red, mean position and surface normal form a plane, this next line will find the surface normal of this
    % i.e. the vector which points directly 'left-to-right' along the LED array      
    % we can use this to calculate the animal's head roll
    roll_normal = normalize( cross(P0-mpos, P0-(mpos+led_normal)),2,'norm' );  
    % negative roll means left ear is lower than the right
    
    % the animal's azimuth is defined by the line from the mean position to the red LED
    % i.e. the vector which points directly 'front-to-back' along the LED array
    % we can use this to calculate the animal's yaw and pitch
    nose_normal = normalize( P0-mpos,2,'norm' );  
    % yaw: [-y 0] = -90, [+y 0] = 90, [+x 0] = 0, [-x 0] = 180
    
    % yaw, pitch and roll from these
    [yaw,pitch,~] = cart2sph(nose_normal(:,1),nose_normal(:,2),nose_normal(:,3));
    [~,roll,~] = cart2sph(roll_normal(:,1),roll_normal(:,2),roll_normal(:,3));
    roll(led_normal(:,3)<0) = roll(led_normal(:,3)<0)+pi; % if the animal is inverted add 180 degrees to roll
    pos.yaw = rad2deg(yaw(:));
    pos.pitch = rad2deg(pitch(:));
    pos.roll = rad2deg(roll(:));
    disp(sprintf('\t...done'))    
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Session data and figure
    % open figure
    fig_clust = figure('Units','pixels','Position',[100 50 1200 800],'visible','off');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white 
    ann_str = sprintf('Rat: %s, Date: %s, Analysed: %s',rnow,dnow,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',10,'LineStyle','none','interpreter','none'); 

    map_now = 'turbo';

    % run through every part
    dwell_maps_3d = cell(1,nparts);
    for pp = 1:nparts % for every part               
        session_times = pdata.session_times;
        session_now = pp;
        idx = pos.pot > session_times(session_now,1) & pos.pot < session_times(session_now,2); % index for position data in this part
           
        % get position data
        ppox = pos.pox(idx);
        ppoy = pos.poy(idx);
        ppoz = pos.poz(idx);       
        ppit = pos.pitch(idx);
        pyaw = pos.yaw(idx);   

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
            ax.CLim = [0,max([0.1 max(F1(:),[],'omitnan')],[],'omitnan')];
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
    end
    pdata.dwell_maps_3d = dwell_maps_3d;

%% >>>>>>>>>> Save the overall figure
    fname = [data_dir '\' rnow '\' dnow '\3Danalysis\HD_dwell.png'];       
    print(gcf,'-dpng','-r250',fname)
    close(gcf); 
    disp(sprintf('\t...done'))  

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Run through clusters
    disp(sprintf('Analysing clusters...'))
    cluma.hd_3d_ratemap = cell(size(cluma,1),1);
    cluma.hd_3d_curves = cell(size(cluma,1),1);
    cluma.hd_3d_info = NaN(size(cluma,1),10);
    dwell_maps_2d = cell(2,nparts);

    ucis = unique(cluma.uci); % list of unique cells in sdata    
    for uu = 1:length(ucis)    
        uci = ucis{uu};
        disp(sprintf('\tCell %d of %d (%.f%%): %s',uu,length(ucis),uu/length(ucis)*100,uci))

        % open figure
        fig_clust = figure('Units','pixels','Position',[100 100 1400 800],'visible','off');
        set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
        set(gcf,'color','w'); % makes the background colour white

        % add an annotation to the figure with some important info
        ann_str = sprintf('Cell: %s, Rat: %s, Date: %s, Tetrode: %d, Cluster: %d, Analysed: %s',uci,rnow,dnow,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',10,'LineStyle','none','interpreter','none');      

        xsiz = 250;
        xbuff = 160;
        ysiz = xsiz/2;
        ybuff = 100;
        yvec = [615 370];
        xnow = 150;
        xvec = [xnow xnow+(xsiz*1)+(xbuff*1) xnow+(xsiz*2)+(xbuff*2)];

        for pp = 1:nparts % for every part               
            part_now = part_names{pp};
            disp(sprintf('\t\t%s',part_now))                   
            
            idx = find( ismember(cluma.uci,uci) & cluma.partn==pp );
            session_times = pdata.session_times;
            pidx = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part
           
            ppot = pos.pot(pidx);
            ppit = pos.pitch(pidx); % in mm            
            pyaw = pos.yaw(pidx); % in mm            

            pspt = cluma.spike_times_s{idx};   
            if isempty(pspt)
                continue
            end
            sindx = knnsearch(ppot,pspt);          
            psyaw = pyaw(sindx);
            pspit = ppit(sindx);

            % plot spikes and positions
            ax = axes('Units','pixels','Position',[xvec(pp) yvec(1) xsiz ysiz]);
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
                ax.CLim = [0,max([0.1 max(F3(:),[],'omitnan')],[],'omitnan')];
                colormap(gca,map_now)
                ax.XTick = [];
                ax.YTick = [];

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
                ratemap = s ./ (d .* (1/50));
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
                hold on
                line([mx1 mx1],ax_az.YLim,'Color','r')

                % directional info shuffle
                % observed values
                smetric = spatialMETRICS(ratemap,d);
                az_vals = [smetric.spatial_information smetric.signal_to_noise];
            
                % shuffle values
                rng(999); % for reproducibility
                iti = 100;
                vals_shuff = NaN(iti,2);
                for ii = 1:iti
                    [~,sindx2] = shift_spike_train(ppot,pspt,'spindx',sindx);
                    s_shuff = histcounts(pyaw(sindx2),edg);
                    ratemap_shuff = s_shuff ./ (d .* (1/50));
                    smetric = spatialMETRICS(ratemap_shuff,d);
                    vals_shuff(ii,:) = [smetric.spatial_information smetric.signal_to_noise];
                end
                az_vals_z = (az_vals - mean(vals_shuff,1,'omitnan')) ./ std(vals_shuff,1,'omitnan');

                text(ax,0,-1.3,sprintf('Azimuth\nspatial info %.2f (z = %.2f)\nSNR %.2f (z = %.2f)\n',az_vals(1),az_vals_z(1),az_vals(2),az_vals_z(2)),'units','normalized','FontSize',10);

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
                ratemap = s ./ (d .* (1/50));
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
                hold on
                line(ax_pi.XLim,[mx2 mx2],'Color','r')

                % directional info shuffle
                % observed values
                smetric = spatialMETRICS(ratemap,d);
                pi_vals = [smetric.spatial_information smetric.signal_to_noise];
            
                % shuffle values
                rng(999); % for reproducibility
                iti = 100;
                vals_shuff = NaN(iti,2);
                for ii = 1:iti
                    [~,sindx2] = shift_spike_train(ppot,pspt,'spindx',sindx);
                    s_shuff = histcounts(ppit(sindx2),edg);
                    ratemap_shuff = s_shuff ./ (d .* (1/50));
                    smetric = spatialMETRICS(ratemap_shuff,d);
                    vals_shuff(ii,:) = [smetric.spatial_information smetric.signal_to_noise];
                end
                pi_vals_z = (pi_vals - mean(vals_shuff,1,'omitnan')) ./ std(vals_shuff,1,'omitnan');
                text(ax,0,-1.9,sprintf('Pitch\nspatial info %.2f (z = %.2f)\nSNR %.2f (z = %.2f)\n',pi_vals(1),pi_vals_z(1),pi_vals(2),pi_vals_z(2)),'units','normalized','FontSize',10);

                % accumulate
                cluma.hd_3d_info(idx,:) = [mx1 mx2 az_vals az_vals_z pi_vals pi_vals_z];
                cluma.hd_3d_curves(idx,1) = { single([ratemap_az; ratemap_pit]) };
                cluma.hd_3d_ratemap(idx,1) = { single(F3) };                
        end

%% >>>>>>>>>> Save the overall figure
        % Save the figure  
        fname = [data_dir '\' rnow '\' dnow '\3Danalysis\part_figures\' uci '_HD_3D.png'];
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
    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    cluma = rmprop(cluma,'pdata');
    cluma = addprop(cluma,{'pdata'},{'table'});
    pdata.pos = pos;
    cluma.Properties.CustomProperties.pdata = pdata;
    cname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    save(cname,'cluma'); % save session data
    analysis_log({'pitch_tuning'},1,'version',{'v1.0.0'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('hill_analysis has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');









































