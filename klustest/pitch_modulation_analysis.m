function pitch_modulation_analysis
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
    close all;

    outname = 'klustest';
    skipfigs = 0; % 1 = skip making a figure if it exists already    
    save_figs = 1; % 1 = make and save figures (VERY time consuming, takes about 90% of klustest's time)
    fast_figures = 1; % 1 = 4-5x faster figure saving, but at a slightly lower quality
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    stk = dbstack;
    tnow = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',stk.name,tnow))
    
%% >>>>>>>>>> Prepare the data
    sname = [pwd '\' outname '\sdata.mat'];
    disp(sprintf('Loading sdata: %s',sname))    
    load(sname,'sdata'); % load saved session data
    pdata = sdata.Properties.CustomProperties.pdata;
    part_config = pdata.part_config;
    nparts = size(part_config,1);
    disp(sprintf('\t...%d sessions',size(pdata.sessions,1)));   
    disp(sprintf('\t...recording date: %s',pdata.date));            
    disp(sprintf('\t...done'))    
    
    %% Map settings    
    mapset = pdata.mapset;

%% >>>>>>>>>> Analyse trajectory and add 3D HD
    disp(sprintf('Loading 3D trajectory...'))       
    for pp = 1:nparts % for every part     
        part_now = part_config.part_names{pp};
        if pp==1
            disp(sprintf('\t\t%s',part_now))                   
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
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Run through clusters
    disp(sprintf('Analysing clusters...'))
    ucis = unique(sdata.uci); % list of unique cells in sdata    

    %% run through every cell
    loopout = looper(length(ucis));  
    sdata.hd_map_3d = cell(size(sdata,1),1);
    sdata.hd_3d_curves = cell(size(sdata,1),1);
    sdata.hd_3d_info = NaN(size(sdata,1),4);
    
    for uu = 1:length(ucis) 
        % open figure
        fig_clust = figure('Units','pixels','Position',[100 100 1200 900],'visible','off');
        set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
        set(gcf,'color','w'); % makes the background colour white

        % add an annotation to the figure with some important info
        uci = ucis{uu};
        disp(sprintf('\t%s',uci))  

        idx = find( ismember(sdata.uci,uci) & sdata.partn==1 );
        ann_str = sprintf('Cell: %s, Rat: %s, Date: %s, Tetrode: %d, Cluster: %d, Analysed: %s',sdata.uci{idx},sdata.rat{idx},sdata.date{idx},sdata.tetrode(idx),sdata.cluster(idx),datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',10,'LineStyle','none','interpreter','none');      

        % run through every part
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
            ppit = rad2deg(double( pdata.(part_now).pit )); % pitch for this part, in rads   
            pyaw = rad2deg(double( pdata.(part_now).yaw )); % pitch for this part, in rads
            pos = [ppox ppoy ppoz ppit pyaw];

            % get spike data
            idx = find( ismember(sdata.uci,uci) & sdata.partn==pp );
            sindx = sdata.spt_pot_index{idx};            
            pspx = ppox(sindx);
            pspy = ppoy(sindx);
            pspz = ppoz(sindx);
            psph = pyaw(sindx);

            % plot trajectory
            ix = sub2ind([nparts 4],pp,1); % index for this subplot
            subplot(4,nparts,ix)
                plot3(ppox,ppoy,ppoz,'k'); hold on;
                plot3(pspx,pspy,pspz,'r.','MarkerSize',15);
                daspect([1 1 1])
                axis xy off tight
                view(0,90);
                set(gca,'FontSize',12)

            % plot 3D tuning curve (flat projection)
            ix = sub2ind([nparts 4],pp,2); % index for this subplot
            subplot(4,nparts,ix)
                npoints = numel(pyaw); % number of random points
                vadjust = 0;    
                [X,Y,Z] = sph2cart(deg2rad(pyaw),deg2rad(ppit+vadjust),ones(npoints,1)); % convert to XYZ
                [Xd,Yd,Zd,F1] = projectEIGS([X,Y,Z],64,10); % project these onto sphere       

                npoints = numel(sindx); % number of random points                
                [X,Y,Z] = sph2cart(deg2rad(pyaw(sindx)),deg2rad(ppit(sindx)+vadjust),ones(npoints,1)); % convert to XYZ
                [Xs,Ys,Zs,F2] = projectEIGS([X,Y,Z],64,10); % project these onto sphere  
                F3 = F2 ./ F1 .* (1/pdata.pos_srate);
                F3(F1<(1/pdata.pos_srate)) = NaN;                
                [TH,PHI,R] = cart2sph(Xs,Ys,Zs);
            
                surf('XData',TH,'YData',PHI,'ZData',ones(size(PHI)),'CData',F3,'EdgeColor','none'); hold on;
                xlabel('Azimuth (radians)')
                ylabel('Pitch (radians)')
                axis xy tight
                daspect([1 1 1])
                view(0,90);
                set(gca,'FontSize',12)

            % plot azimuthal tuning curve
            ix = sub2ind([nparts 4],pp,3); % index for this subplot
            subplot(4,nparts,ix)
                edg = linspace(-180,180,60);
                xi = movmean(edg,2,'EndPoints','discard');
                d = histcounts(pyaw,edg);
                s = histcounts(pyaw(sindx),edg);
                ratemap = s ./ (d .* (1/pdata.pos_srate));
                ratemap_az = interp1(xi,ratemap,edg,'linear',NaN);
                
                bar(edg,ratemap_az,1,'k')
                xlabel('Azimuth (deg)')
                ylabel('Firing Rate (Hz)')  
                set(gca,'FontSize',12)

                % head direction analyses
                hd3n = ratemap_az ./ max(ratemap_az); % normalise cell hd
                hd3n = hd3n(:);            
                r1 = circ_r(deg2rad(edg(:)),hd3n); % rayleigh vector length
                mx1 = edg(hd3n == max(hd3n)); % preferred angle (location of max frate)
                if length(mx1)>1
                    mx1 = mx1(1);
                elseif isempty(mx1)
                    mx1 = NaN;
                end
                text(0,1.05,sprintf('r: %.2f, PFD%c: %.2f',r1,176,mx1),'Units','normalized','FontSize',10,'HorizontalAlignment','left')
                hold on
                ax = gca;
                line([mx1 mx1],ax.YLim,'Color','r')
                ax.XTick = [-180 -90 0 90 180];

            % plot pitch tuning curve
            ix = sub2ind([nparts 4],pp,4); % index for this subplot
            subplot(4,nparts,ix)
                edg = linspace(-90,90,60);
                xi = movmean(edg,2,'EndPoints','discard');
                d = histcounts(ppit,edg);
                s = histcounts(ppit(sindx),edg);
                ratemap = s ./ (d .* (1/pdata.pos_srate));
                ratemap_pit = interp1(xi,ratemap,edg,'linear',NaN);
                
                bar(edg,ratemap_pit,1,'k')
                xlabel('Pitch (deg)')
                ylabel('Firing Rate (Hz)')  
                set(gca,'FontSize',12)

                % head direction analyses
                hd3n = ratemap_pit ./ max(ratemap_pit); % normalise cell hd
                hd3n = hd3n(:);            
                r2 = circ_r(deg2rad(edg(:).*2),hd3n); % angle doubled rayleigh vector length
                mx2 = edg(hd3n == max(hd3n)); % preferred angle (location of max frate)
                if length(mx2)>1
                    mx2 = mx2(1);
                elseif isempty(mx2)
                    mx2 = NaN;
                end
                text(0,1.05,sprintf('r: %.2f, PFD%c: %.2f',r2,176,mx2),'Units','normalized','FontSize',10,'HorizontalAlignment','left')
                hold on
                ax = gca;
                line([mx2 mx2],ax.YLim,'Color','r')
                ax.XTick = [-90 -45 0 45 90];     

                %% add to sdata
                sdata.hd_map_3d(idx) = { single(F3) }; % azimuth x pitch tuning curve
                sdata.hd_3d_curves(idx) = { single([ratemap_az(:)'; ratemap_pit(:)']) }; % azimuth tuning curve; pitch tuning curve
                sdata.hd_3d_info(idx,:) = single([r1 r2 mx1 mx2]); % azimuth rayleigh, pitch rayleigh, azimuth pfd, pitch pfd
        end

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
            exportgraphics(fig_clust,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',150);  
        end
        close(fig_clust);    

    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    save([pwd '\' pdata.outname '\sdata.mat'],'sdata'); % save session data
    analysis_log({'pitch_modulation_analysis'},1,'version',{'v1.0.0'});

    % finish up
    toc1 = toc/60;
    disp(sprintf('%s has finished. It took %0.3f seconds or %0.3f minutes',stk.name,toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',[pwd '\' pdata.outname],' &'');">','klustest folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');







% 
%     for pp = 1:nparts % for every part   
% pp = 1      
%         part_now = part_config.part_names{pp};
%         disp(sprintf('\t%s',part_now))      
% 
%         % get position data
%         ppox = double( pdata.(part_now).pox ).*10; % pos x for this part, in mm
%         ppoy = double( pdata.(part_now).poy ).*10; % pos y for this part, in mm
%         ppoz = double( pdata.(part_now).poz ).*10; % pos z for this part, in mm          
%         ppit = double( pdata.(part_now).pit ); % pitch for this part, in rads   
%         pyaw = double( pdata.(part_now).yaw ); % pitch for this part, in rads
%         pos = [ppox ppoy ppoz ppit pyaw];
% 
%         % filter by pitch angle
%         med_val = median(ppit(:),'omitnan');   
%         pos_up = pos;
%         pos_up(ppit<med_val,:) = NaN; % keep only data where nose is pitched above med_val
%         pos_dn = pos;
%         pos_dn(ppit>med_val,:) = NaN; % keep only data where nose is pitched below med_val
% 
%         % filter by yaw
%         yaw_val = deg2rad(90);
%         pos_left = pos;
%         pos_left(abs(pyaw)<yaw_val,:) = NaN; % keep only data where nose is pitched above med_val
%         pos_right = pos;
%         pos_right(abs(pyaw)>yaw_val,:) = NaN; % keep only data where nose is pitched below med_val
% 
%         loopout = looper(length(ucis));
%         for uu = 1:length(ucis)     
% uu = 5            
%             uci = ucis{uu};
%             disp(sprintf('\t%s',uci))      
% 
%             idx = find( ismember(sdata.uci,uci) & sdata.partn==pp );
%             sindx = sdata.spt_pot_index{idx};
%             pspx = ppox(sindx);
%             pspy = ppoy(sindx);
%             pspz = ppoz(sindx);
%             psph = pyaw(sindx);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% figure('Units','pixels','Position',[100 100 1200 700])
% set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
% set(gcf,'color','w'); % makes the background colour white
% 
% subplot(1,2,1)
% c = cline(ppox,ppoy,ppoz,pyaw); hold on;
% daspect([1 1 1])
% axis xy
% 
% subplot(1,2,2)
% scatter3(pspx,pspy,pspz,20,psph,'filled')
% daspect([1 1 1])
% axis xy
% 
%       npoints = numel(pyaw); % number of random points
%       vadjust = 0;
%       [X,Y,Z] = sph2cart(deg2rad(pyaw),deg2rad(ppit+vadjust),ones(npoints,1)); % convert to XYZ
%       [Xd,Yd,Zd,F1] = projectEIGS([X,Y,Z],64,10); % project these onto sphere
% 
%       npoints = numel(psph); % number of random points
%       [X,Y,Z] = sph2cart(deg2rad(pyaw(sindx)),deg2rad(ppit(sindx)+vadjust),ones(npoints,1)); % convert to XYZ
%       [Xs,Ys,Zs,F2] = projectEIGS([X,Y,Z],64,10); % project these onto sphere
%       F3 = F2 ./ F1 .* (1/pdata.pos_srate);
%       F3(F1<(1/pdata.pos_srate)) = NaN;
% 
% figure('Units','pixels','Position',[100 100 1200 700])
% set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
% set(gcf,'color','w'); % makes the background colour white
% colormap(gcf,'turbo')
%       subplot(2,3,1)
%       plot3(X,Y,Z,'ko')
%       daspect([1 1 1])
%       axis xy tight
% 
%       subplot(2,3,2)
%       surf(Xs,Ys,Zs,F1,'EdgeColor','none');
%       daspect([1 1 1])
%       axis xy tight
%       title('Dwell time')
% 
%       subplot(2,3,3)
%       surf(Xs,Ys,Zs,F2,'EdgeColor','none');
%       daspect([1 1 1])
%       axis xy tight
%       title('Spikes')
% 
%       subplot(2,3,4)
%       surf(Xs,Ys,Zs,F3,'EdgeColor','none');
%       daspect([1 1 1])
%       axis xy tight
%       title('Rate')
% 
%         subplot(2,3,5)
%         [x1,y1,z1] = sphere(32);
%         surf(x1,y1,z1,ones(size(x1)),'EdgeColor',[0 0 0],'EdgeAlpha',0.2,'FaceColor',[0 0 0],'FaceAlpha',0.1); hold on;
%         daspect([1 1 1])
%         axis off xy tight
%         line([-1 1],[0 0],[0 0],'Color','k')
%         line([0 0],[-1 1],[0 0],'Color','k')
%         line([0 0],[0 0],[-1 1],'Color','k')
%         [TH,PHI,R] = cart2sph(Xs,Ys,Zs);
%         [X,Y,Z] = sph2cart(TH,PHI,F3./max(F3(:)));
%         surf(X,Y,Z,F3,'EdgeColor','none','FaceColor','flat','AlphaData',F3); hold on;
%         ax = gca;
%         ax.CLim = [0 max(F3(:))];
%         view(-120,20);
% 
%         subplot(2,3,6)
%         surf('XData',TH,'YData',PHI,'ZData',ones(size(PHI)),'CData',F3,'EdgeColor','none'); hold on;
%         xlabel('Azimuth (radians)')
%         ylabel('Pitch (radians)')
%         axis xy tight
%         daspect([1 1 1])
%         view(0,90);
% 
% exportgraphics(gcf, "C:\Users\F004KS7\Downloads\cell1_fig1.png",'Resolution',300)
% 
% 
% 
% %       keyboard
% exportgraphics(gcf, "C:\Users\F004KS7\Downloads\cell1_fig2.png",'Resolution',300)
% 
% 
% 
% 
% 
% keyboard




























% 
% 
% %  keyboard           
%             idx = find( ismember(sdata.uci,uci) & sdata.partn==pp );
% 
% %% >>>>>>>>>> Firing rate maps            
%             %% pitch up ratemap     
%             rmset = mapset;
%             lx = pdata.maze_frame(:,1) .* 10; % in mm
%             ly = pdata.maze_frame(:,2) .* 10; % in mm
%             lz = pdata.maze_frame(:,3) .* 10; % in mm
% 
%             spk_up = pos_up(sdata.spt_pot_index{idx},1:2);
%             rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
%             rmset.srate = pdata.pos_srate;
%             if isfield(pdata.(part_now),'pitch_up_speedlift') % if a previously computed dwellmap exists, use this to save time
%                 speedlift = pdata.(part_now).pitch_up_speedlift{1};
%                 [ratemap_up,dwellmap_up, ~, ~, ~] = rate_mapper(pos_up(:,1:2),spk_up,rmset,speedlift);             
%             else
%                 [ratemap_up,dwellmap_up, ~, ~,speedlift] = rate_mapper(pos_up(:,1:2),spk_up,rmset);              
%                 pdata.(part_now).pitch_up_speedlift = { single(speedlift) };
%             end
%             %sdata.ratemap_pitch_up(idx) = { single(ratemap_up) };
%             %pdata.(part_now).dwellmap_pitch_up = single(dwellmap_up);
% 
%             %% pitch down ratemap     
%             spk_dn = pos_dn(sdata.spt_pot_index{idx},1:2);
%             rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
%             rmset.srate = pdata.pos_srate;
%             if isfield(pdata.(part_now),'pitch_down_speedlift') % if a previously computed dwellmap exists, use this to save time
%                 speedlift = pdata.(part_now).pitch_down_speedlift{1};
%                 [ratemap_dn,dwellmap_dn, ~, ~, ~] = rate_mapper(pos_dn(:,1:2),spk_dn,rmset,speedlift);             
%             else
%                 [ratemap_dn,dwellmap_dn, ~, ~,speedlift] = rate_mapper(pos_dn(:,1:2),spk_dn,rmset);              
%                 pdata.(part_now).pitch_down_speedlift = { single(speedlift) };
%             end
%             %sdata.ratemap_pitch_down(idx) = { single(ratemap_dn) };
%             %pdata.(part_now).dwellmap_pitch_down = single(dwellmap_dn);        
% 
%             %% yaw ratemap 1
%             spk_left = pos_left(sdata.spt_pot_index{idx},1:2);
%             rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
%             rmset.srate = pdata.pos_srate;
%             if isfield(pdata.(part_now),'pitch_up_speedlift') % if a previously computed dwellmap exists, use this to save time
%                 speedlift = pdata.(part_now).pitch_up_speedlift{1};
%                 [ratemap_left,dwellmap_left, ~, ~, ~] = rate_mapper(pos_left(:,1:2),spk_left,rmset,speedlift);             
%             else
%                 [ratemap_left,dwellmap_left, ~, ~,speedlift] = rate_mapper(pos_left(:,1:2),spk_left,rmset);              
%                 pdata.(part_now).pitch_up_speedlift = { single(speedlift) };
%             end
%             %sdata.ratemap_left(idx) = { single(ratemap_left) };
%             %pdata.(part_now).dwellmap_left = single(dwellmap_left);
% 
%             %% yaw ratemap 2   
%             spk_right = pos_right(sdata.spt_pot_index{idx},1:2);
%             rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
%             rmset.srate = pdata.pos_srate;
%             if isfield(pdata.(part_now),'pitch_down_speedlift') % if a previously computed dwellmap exists, use this to save time
%                 speedlift = pdata.(part_now).pitch_down_speedlift{1};
%                 [ratemap_right,dwellmap_right, ~, ~, ~] = rate_mapper(pos_right(:,1:2),spk_right,rmset,speedlift);             
%             else
%                 [ratemap_right,dwellmap_right, ~, ~,speedlift] = rate_mapper(pos_right(:,1:2),spk_right,rmset);              
%                 pdata.(part_now).pitch_down_speedlift = { single(speedlift) };
%             end
%             %sdata.ratemap_right(idx) = { single(ratemap_right) };
%             %pdata.(part_now).dwellmap_right = single(dwellmap_right);             
% 
%             %% ratemap shuffles
%             iti = 100;
%             rs = NaN(iti,1);
%             for ii = 1:iti
%                 rindx = randperm(size(pos,1),size(pos,1)); % random index into position data
%                 rindx = rindx(:);
%                 rindx1 = rindx(1:2:end); % one half of random data
%                 rindx2 = rindx(2:2:end); % second half of random data
%                 spike_indx = sdata.spt_pot_index{idx};
%                 sindx = ismember(spike_indx,rindx1);
% 
%                 pos1 = pos(rindx1,1:2);
%                 spk1 = pos(spike_indx(sindx),1:2);
%                 pos2 = pos(rindx2,1:2);
%                 spk2 = pos(spike_indx(~sindx),1:2);
% 
%                 [ratemap1,~, ~, ~, ~] = rate_mapper(pos1,spk1,rmset);             
%                 [ratemap2,~, ~, ~, ~] = rate_mapper(pos2,spk2,rmset);             
%                 rs(ii) = corr(ratemap1(:),ratemap2(:),'rows','pairwise','type','Pearson');
%             end
% 
%             %% add correlations to sdata
%             map1 = ratemap_up;
%             map2 = ratemap_dn;
%             r = corr(map1(:),map2(:),'rows','pairwise','type','Pearson');  
%             z = (r-mean(rs(:),'omitnan')) / std(rs(:),'omitnan');
%             if ~any(ismember(sdata.Properties.VariableNames,'pitch_mod_r')) % if the column(s) do not exist yet
%                 sdata.pitch_mod_r = NaN(size(sdata,1),2); % preallocate
%             end               
%             sdata.pitch_mod_r(idx,:) = [r z];        
% 
%             map1 = ratemap_left;
%             map2 = ratemap_right;
%             r = corr(map1(:),map2(:),'rows','pairwise','type','Pearson');  
%             z = (r-mean(rs(:),'omitnan')) / std(rs(:),'omitnan');
%             if ~any(ismember(sdata.Properties.VariableNames,'yaw_mod_r')) % if the column(s) do not exist yet
%                 sdata.yaw_mod_r = NaN(size(sdata,1),2); % preallocate
%             end               
%             sdata.yaw_mod_r(idx,:) = [r z];  
% keyboard
%             loopout = looper(loopout);            
%         end
% keyboard        
% %         disp(sprintf('\t...saving figure'))              
% %         fname = [pwd '\' pdata.outname '\' part_now '_repetition.png'];
% %         if (exist(fname,'file') && skipfigs) || ~save_figs  % if the figure exists and we don't want to overwrite it
% %             % do nothing
% %         else
% %             klustfig_repetition(pdata,sdata,pp,fast_figures,'off');
% %         end    
% %         disp(sprintf('\t...done'))                      
%     end
% 
% %% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
%     sdata = rmprop(sdata,'pdata');
%     sdata = addprop(sdata,{'pdata'},{'table'});
%     sdata.Properties.CustomProperties.pdata = pdata;
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






















