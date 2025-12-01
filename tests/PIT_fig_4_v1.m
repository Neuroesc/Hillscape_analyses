%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% PIT_fig_1_v1  figure script written for:
% Grieves, Duvelle and Taube (202X) 
%
% See also: GIT_audit

% HISTORY:
% version 1.0.0, Release 06/11/23 Code conception
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2023 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Field anisotropy example cell
    % collect data
    ucis = unique(clumaa.uci(pidx));
    datn = cell(1,3);
    for pp=1:3    
        dat = [];
        for uu = 1:length(ucis) % for each cell
            idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==pp;
            if pp==2
                fdata = clumaa.surficial_field_data{idx}; % field data for this cell in this session in this part
                rmap = clumaa.ratemap_surficial{idx};       
            else
                fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
                rmap = clumaa.ratemap_planar{idx};                       
            end
            if isempty(fdata)
                continue
            end

            if pp==2
                scale_x = 116 ./ size(rmap,2); % we need to scale field centroids to match a common map size
            else
                scale_x = 96 ./ size(rmap,2); % we need to scale field centroids to match a common map size
            end

            b = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0]+0.5; 
            a = [90 0 90 0];
            gpoly = [];
            res = 100;
            for jj = 1:size(b,1)-1
                gpoly = [gpoly; linspace(b(jj,1),b(jj+1,1),res)' linspace(b(jj,2),b(jj+1,2),res)' repmat(a(jj),res,1)];
            end

            ps = linspace(0,size(rmap,2),7);
            tops = ps(2:2:end);
            b2 = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0;tops(1) 0;tops(1) size(rmap,1);tops(2) size(rmap,1);tops(2) 0;tops(3) 0;tops(3) size(rmap,1)]+0.5; 
            a2 = [90 0 90 0 0 90 0 90 0 90];
            gpoly2 = [];
            res = 100;
            for jj = 1:size(b2,1)-1
                gpoly2 = [gpoly2; linspace(b2(jj,1),b2(jj+1,1),res)' linspace(b2(jj,2),b2(jj+1,2),res)' repmat(a2(jj),res,1)];
            end

            ds = NaN(size(fdata,1),4);
            for jj = 1:size(fdata,1)
                [i,d] = knnsearch(gpoly(:,[1:2]),[fdata.Centroid(jj,1),fdata.Centroid(jj,2)],'K',1);
                ds(jj,1:2) = [d gpoly(i,3)];

                [i,d] = knnsearch(gpoly2(:,[1:2]),[fdata.Centroid(jj,1),fdata.Centroid(jj,2)],'K',1);
                ds(jj,3:4) = [d gpoly2(i,3)];                
            end
            dat = [dat; fdata.MajorAxisLength(:,1) fdata.MinorAxisLength(:,1) fdata.Orientation(:,1) fdata.WeightedCentroid(:,1).*scale_x fdata.WeightedCentroid(:,2) fdata.Area(:,1) ds(:,1) repmat(uu,size(fdata,1),1) ds(:,2) ds(:,3:4)];
        end
        datn(pp) = { dat };
    end
    
    xnow = 20;
    ynow = 600;
    xbuff = 10;
    xbuff2 = 310;
    xsiz = 125;
    ysiz = 70;

    uci = 'RG26_230119_t3_c5';
    m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
    m2 = clumaa.ratemap_surficial{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
    m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('d',sprintf('Example place field anisotropy'),'yoffset',0,'xoffset',20,'width',400);  
    
        imagesc(m1,'alphadata',ones(size(m1)).*0.4); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,turbo)   
        ax1.CLim = [0 max([max(m1(:),[],'omitnan') 1])];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        % text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        % text(0,1,sprintf('Example fields'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        
        % ellipses
        cm = cmocean('balance');
        % cm = cmocean('thermal');                
        cm = cm(32:224,:);
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        cidx = linspace(-180,180,size(cmap,1));
        f = clumaa.planar_field_data{(ismember(clumaa.uci,uci) & clumaa.partn==1)};
        for ff = 1:size(f,1)
            e = ellipse(f.MajorAxisLength(ff,1)./2,f.MinorAxisLength(ff,1)./2,-deg2rad(f.Orientation(ff,1)),f.Centroid(ff,1),f.Centroid(ff,2),'Clipping','off');
            e.LineWidth = 2.5;
            c_idx = knnsearch(cidx(:),f.Orientation(ff,1));
            e.Color = cmap(c_idx,:);
        end

    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
        imagesc(m2,'alphadata',ones(size(m2)).*0.4); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,turbo)   
        ax2.CLim = [0 max([max(m2(:),[],'omitnan') 1])];
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(m2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

        % text
        % text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

        % ellipses
        f = clumaa.surficial_field_data{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
        for ff = 1:size(f,1)
            e = ellipse(f.MajorAxisLength(ff,1)./2,f.MinorAxisLength(ff,1)./2,-deg2rad(f.Orientation(ff,1)),f.Centroid(ff,1),f.Centroid(ff,2),'Clipping','off');
            e.LineWidth = 2.5;
            c_idx = knnsearch(cidx(:),f.Orientation(ff,1));
            e.Color = cmap(c_idx,:);
        end

%% >>>>>>>>>> Field anisotropy map
    xnow = xnow+xbuff2;
    
    % get raw data
    var = 'field_anisotropy';
    v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
    d1 = cat(1,v1{:});

    v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data
    d2 = cat(1,v2{:});

    v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 2 data
    d3 = cat(1,v3{:});

    % map data
    % arena 1
    dist_cutoff = 160;
    mapset.binsize = 32; % (mm) firing rate map bin size
    
    dist_cutoff_bins = dist_cutoff ./ mapset.binsize;    
    anisotropy_map_1 = NaN(48,96);
    mapXY = d1(:,3:4);
    for bb = 1:numel(anisotropy_map_1)
        [r,c] = ind2sub(size(anisotropy_map_1),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_1(bb) = mean(d1(box_idx,1),1,'omitnan');
    end
    anisotropy_map_1 = imgaussfilt(anisotropy_map_1,1);

    % hills  
    anisotropy_map_2 = NaN(48,116);
    mapXY = d2(:,3:4);
    for bb = 1:numel(anisotropy_map_2)
        [r,c] = ind2sub(size(anisotropy_map_2),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_2(bb) = mean(d2(box_idx,1),1,'omitnan');
    end
    anisotropy_map_2 = imgaussfilt(anisotropy_map_2,1);

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
        ah = add_panel_title('f',sprintf('Avg. place field anisotropy'),'yoffset',0,'xoffset',20,'width',400);  
    
        imagesc(anisotropy_map_1,'alphadata',~isnan(anisotropy_map_1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,cmocean('balance'))
        % colormap(ax1,cmocean('thermal'))                
        ax1.CLim = [-0.5 0.5];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,1,sprintf('N = %d fields',size(d1,1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
        imagesc(anisotropy_map_2,'alphadata',~isnan(anisotropy_map_2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)
        ax2.CLim = [-0.5 0.5];
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(anisotropy_map_2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

        % text
        text(0,1,sprintf('N = %d fields',size(d2,1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+200 ynow+80 70 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Anisotropy score'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')
























        return
%% >>>>>>>>>> Field anisotropy analysis (if needed)
    if ~any(ismember(clumaa.Properties.VariableNames,'field_anisotropy')) % if the column(s) do not exist yet
        clumaa.field_anisotropy = cell(size(clumaa,1),1);

        for ss = 1:size(posdata,1) % for each recording session
            disp(sprintf('\tSession %d of %d (%.f%%)',ss,size(posdata,1),ss/size(posdata,1)*100))
            sidx = ismember(clumaa.pos_idx(:),posdata.pos_idx(ss));
            session_times = posdata.session_times{ss};
    
            ucis = unique(clumaa.uci(sidx));
    
            for pp = 1:size(session_times,1) % for each part
                for uu = 1:length(ucis) % for each cell recorded in this session
                    idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==pp;
                    if pp==2
                        fdata = clumaa.surficial_field_data{idx}; % field data for this cell in this session in this part
                        rmap = clumaa.ratemap_surficial{idx};
                        bmap = posdata.anisotropy_map{ss,pp};
                    else
                        fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
                        rmap = clumaa.ratemap_planar{idx};       
                        bmap = posdata.anisotropy_map{ss,pp};
                    end
                    if isempty(fdata)
                        continue
                    end
    
                    % get field anisotropy
                    fanisotropy = (fdata.BoundingBox(:,4) - fdata.BoundingBox(:,3)) ./  (fdata.BoundingBox(:,4) + fdata.BoundingBox(:,3)); % height - width / height + width
    
                    % get behaviour anisotropy
                    % fcents = fdata(:,2:3); % centroids
                    fcents = fdata.WeightedCentroid(:,1:2); % weighted centroids
                    % banisotropy = interp2(bmap,fcents(:,1),fcents(:,2),'nearest');
    
                    banisotropy = NaN(size(fcents,1),1);
                    for ff = 1:size(fdata,1) % for every field
                        pxls = fdata.PixelIdxList{ff,1};
                        banisotropy(ff) = median(reshape(bmap(pxls),[],1),'omitnan');
                    end
                    
                    % accumulate
                    if pp==2
                        scale_x = 116 ./ size(bmap,2); % we need to scale field centroids to match a common map size
                    else
                        scale_x = 1;
                    end
                    clumaa.field_anisotropy(idx) = { single([fanisotropy banisotropy fdata.WeightedCentroid(:,1).*scale_x fdata.WeightedCentroid(:,2)]) };
    
                    if 0
                        figure
                        subplot(1,2,1)
                        imagesc(rmap); hold on;
                        plot(fcents(:,1),fcents(:,2),'ko')
                        daspect([1 1 1])
                        axis xy
        
                        subplot(1,2,2)
                        imagesc(bmap); hold on;
                        plot(fcents(:,1),fcents(:,2),'ko')
                        daspect([1 1 1])
                        axis xy
                        keyboard
                    end
                end
            end
        end
    end

%% >>>>>>>>>> Behaviour anisotropy example trajectory
    xnow = -5;
    ynow = 680;
    xbuff = 210;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+50 220 90]);
        ah = add_panel_title('a',sprintf('Example trajectory'),'yoffset',-10,'xoffset',45,'width',400);  

        % get data
        eg_session = 4;
        pos = posdata.pos{eg_session};
        stimes = posdata.session_times{eg_session};
        mf = posdata.maze_frame{eg_session};

        pp = 1;
        idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
        pox = pos.pox_planar(idx);
        poy = pos.poy_planar(idx);
        poz = pos.poz_planar(idx);
        poh = pos.yaw(idx);

        % plot data
        plot3(mf(:,1),mf(:,2),mf(:,3),'k'); hold on;
        c = cline(pox,poy,poz,poh);
        lw = 1.5;
        set(c,'LineWidth',lw);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        cm = cmocean('balance');
        % cm = cmocean('thermal');        
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        colormap(ax1,cmap)
        view(0,90)
        ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax1.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax1.CLim = [-180 180];

        % text
        % text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        % text(0,1,sprintf('Example session'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
        % get data
        pp = 2;
        idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
        pox = pos.pox_surficial(idx);
        poy = pos.poy_surficial(idx);
        poz = pos.poz_surficial(idx);
        poh = pos.yaw(idx);

        % plot data
        plot3(mf(:,4),mf(:,2),mf(:,3),'k'); hold on;        
        c = cline(pox,poy,poz,poh);
        set(c,'LineWidth',lw);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)
        view(0,90)
        ax2.XLim = [min(mf(:,4)) max(mf(:,4))];
        ax2.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax2.CLim = [-180 180];

        % additional plots
        ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

    % plot arena 2
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+xbuff ax2.Position(2) ax1.Position(3) ax1.Position(4)]);
        % get data
        pp = 3;
        idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
        pox = pos.pox_planar(idx);
        poy = pos.poy_planar(idx);
        poz = pos.poz_planar(idx);
        poh = pos.yaw(idx);

        % plot data
        plot3(mf(:,1),mf(:,2),mf(:,3),'k'); hold on;
        c = cline(pox,poy,poz,poh);
        set(c,'LineWidth',lw);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        cm = cmocean('balance');
        % cm = cmocean('thermal');        
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        colormap(ax3,cmap)
        view(0,90)
        ax3.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax3.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax3.CLim = [-180 180];

    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax2.Position(1)+20 ax1.Position(2)+100 80 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Direction (%c)',176),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')







return

% %% check 3D heading maps exist
%     [clumaa,posdata] = PIT_3D_heading(config,clumaa,posdata,0);
% 
% %% >>>>>>>>>> Directional sampling in each maze, example sessions
%     xnow = 40;
%     ynow = 770;
%     xbuff = -10;
%     xsiz = 150;
%     ysiz = 80;
% 
%     % plot example arena session
%     ax1a = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
%         ah = add_panel_title('a',sprintf('Example trajectory'),'yoffset',-20,'xoffset',0,'width',400);  
% 
%         % get data
%         eg_session = 3;
%         pos = posdata.pos{eg_session};   
%         stimes = posdata.session_times{eg_session};
%         mf = posdata.maze_frame{eg_session};
% 
%         pp = 1;
%         idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
%         pox = pos.pox_surficial(idx);
%         poy = pos.poy_surficial(idx);
%         poz = pos.poz_surficial(idx);
%         poh = pos.yaw(idx);
%         mf = posdata.maze_frame{eg_session};
% 
%         % plot data
%         plot3(mf(:,1),mf(:,2),mf(:,3),'k'); hold on;        
%         c = cline(pox,poy,poz,poh);
%         set(c,'LineWidth',1);
% 
%         % axis settings
%         axis xy off tight
%         daspect([1 1 1])
%         cm = cmocean('balance');
%         cmap = [cm; flipud(cm);cm; flipud(cm)];
%         colormap(ax1a,cmap)
%         view(0,90)
%         ax1a.XLim = [min(mf(:,4)) max(mf(:,4))];
%         ax1a.YLim = [min(mf(:,2)) max(mf(:,2))];
%         ax1a.FontSize = 8;
% 
%         % additional plots
%         ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
%         tops = ps(2:2:end);
%         plot([tops; tops],ax1a.YLim,'k:')
% 
%         % text
%         text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',ax1a.FontSize,'Color','k','VerticalAlignment','bottom','rotation',90)
% 
%     % plot example hills session
%     xnow = xnow+240;
%     ax2a = axes('Units','pixels','Position',[ax1a.Position(1)+ax1a.Position(3)+xbuff ax1a.Position(2) ax1a.Position(3) ax1a.Position(4)]);
%         pp = 2;
%         idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
%         pox = pos.pox_planar(idx);
%         poy = pos.poy_planar(idx);
%         poz = pos.poz_planar(idx);
%         poh = pos.yaw(idx);
% 
%         % plot data
%         plot3(mf(:,1),mf(:,2),mf(:,3),'k'); hold on;        
%         c = cline(pox,poy,poz,poh);
%         set(c,'LineWidth',1);
% 
%         % axis settings
%         axis xy off tight
%         daspect([1 1 1])
%         cm = cmocean('balance');
%         cmap = [cm; flipud(cm);cm; flipud(cm)];
%         colormap(ax2a,cmap)
%         view(0,90)
%         ax2a.XLim = [min(mf(:,4)) max(mf(:,4))];
%         ax2a.YLim = [min(mf(:,2)) max(mf(:,2))];
%         ax2a.YTick = [];
%         ax2a.FontSize = 8;
% 
%         % additional plots
%         ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
%         tops = ps(2:2:end);
%         plot([tops; tops],ax2a.YLim,'k:')
% 
%         % text
%         text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',ax2a.FontSize,'Color','k','VerticalAlignment','bottom','rotation',90)

    % % circular colorbar
    % axc = axes('Units','pixels','Position',[ax2a.Position(1)+ax2a.Position(3)-60 ynow+80 30 30]);
    %     m = ones(101,101);
    %     [rr,cc] = ndgrid(1:101,1:101);
    %     c = [51 51];
    %     n = numel(rr(:));
    %     v1 = [cc(:) rr(:)];
    %     v2 = repmat(c,size(v1,1),1);
    %     theta = angle2Points(v1,v2);
    %     a = reshape(theta,size(rr));
    %     a = rad2deg(a);
    % 
    %     dmat = zeros(size(m));
    %     dmat(c(1),c(2)) = 1;
    %     d = bwdist(dmat);
    %     a(d>50) = NaN;
    % 
    %     % plot image
    %     imagesc(a,'alphadata',~isnan(a))
    %     colormap(axc,ax1a.Colormap);
    % 
    %     % axis settings
    %     axc.FontSize = 7;
    %     axis off tight
    %     daspect([1 1 1])
    % 
    %     % text
    %     text(0.5,1,sprintf('Direction'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

% %% >>>>>>>>>> Directional sampling in each maze, 3D HD heatmaps
%     xbuff = -5;
%     ybuff = 90;
%     xsiz = 125;
%     ysiz = 70;
% 
%     % plot arena 3D HD map
%     ax1b = axes('Units','pixels','Position',[ax1a.Position(1) ax1a.Position(2)-ybuff xsiz ysiz]);
%         % get data
%         % dwellmap = posdata.hd_3d_dwell{eg_session,1}; % HD
%         dwellmap = posdata.mov_3d_dwell{eg_session,1}; % movement 
% 
%         % plot data
%         imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))
% 
%         % axis settings
%         xlabel(sprintf('Azimuth (%c)',176))
%         ylabel(sprintf('Pitch (%c)',176))
%         axis xy tight
%         view(0,90);
%         ax1b.FontSize = 8;
%         ax1b.CLim = [0,max(dwellmap(:))];
%         colormap(gca,turbo)
%         ax1b.XTick = -180:90:180;
%         ax1b.YTick = -90:45:90;
% 
%     % horizontal colormap
%     axc = axes('Units','pixels','Position',[ax1b.Position(1)+20 ax1b.Position(2)+75 80 8]);
%         x = linspace(ax1a.CLim(1),ax1a.CLim(2),100);
%         imagesc(x,'XData',ax1a.CLim,'YData',[0 1]);
%         colormap(axc,ax1a.Colormap);
%         axis on
%         axc.XTick = [];
%         axc.YTick = [];
%         axc.FontSize = 7;
%         text(0,0.5,sprintf('%.f ',ax1a.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
%         text(1,0.5,sprintf(' %.f',ax1a.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
%         text(0.5,1,sprintf('Direction (%c)',176),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')
% 
%     % plot hills 3D HD map
%     ax2b = axes('Units','pixels','Position',[ax2a.Position(1) ax2a.Position(2)-ybuff xsiz ysiz]);
%         % get data
%         % dwellmap = posdata.hd_3d_dwell{eg_session,2}; % HD
%         dwellmap = posdata.mov_3d_dwell{eg_session,2}; % movement
% 
%         % plot data
%         imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))
% 
%         % axis settings
%         xlabel(sprintf('Azimuth (%c)',176))
%         axis xy tight
%         view(0,90);
%         ax2b.FontSize = 8;
%         ax2b.CLim = [0,max(dwellmap(:))];
%         colormap(gca,turbo)
%         ax2b.XTick = ax1b.XTick;
%         ax2b.YTick = ax1b.YTick;        
%         ax2b.YTickLabel = {};
% 
%     % horizontal colormap
%     axc = axes('Units','pixels','Position',[ax2b.Position(1)+20 ax2b.Position(2)+75 70 8]);
%         x = linspace(ax1b.CLim(1),ax1b.CLim(2),100);
%         imagesc(x,'XData',ax1b.CLim,'YData',[0 1]);
%         colormap(axc,ax1b.Colormap);
%         axis on
%         axc.XTick = [];
%         axc.YTick = [];
%         axc.FontSize = 7;
%         text(0,0.5,sprintf('%.1f ',ax1b.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
%         text(1,0.5,sprintf(' %.1f',ax1b.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
%         text(0.5,1,sprintf('Dwell time (s)'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')
% 
% %% >>>>>>>>>> Directional sampling in each maze, averaged
%     xnow = xnow+80;
%     xbuff = 10;
%     ybuff = 80;
%     xsiz = 125;
%     ysiz = 70;
%     ysiz2 = 60;
%     ybuff2 = 10;
% 
%     % plot arena 3D HD map (averaged)
%     ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
%         ah = add_panel_title('b',sprintf('Avg. directional sampling'),'yoffset',-10,'xoffset',0,'width',400);  
% 
%         % get data
%         % dwellmap = posdata.hd_3d_dwell(:,1); % HD
%         dwellmap = posdata.mov_3d_dwell(:,1); % movement
%         dwellmap = cat(3,dwellmap{:});
%         for ii = 1:size(dwellmap,3)
%             dwellmap(:,:,ii) = (dwellmap(:,:,ii) - mean(dwellmap(:,:,ii),'all','omitmissing')) ./ std(dwellmap(:,:,ii),[],'all','omitmissing');
%         end
%         dwellmap_a = mean(dwellmap,3,'omitmissing');
% 
%         % plot data
%         imagesc([-180 180],[-90 90],dwellmap_a,'alphadata',~isnan(dwellmap_a))
% 
%         % axis settings
%         ylabel(sprintf('Pitch (%c)',176))
%         axis xy tight
%         view(0,90);
%         ax1.FontSize = 8;
%         ax1.CLim = [-1,6];
%         colormap(gca,turbo)
%         ax1.XTick = [];
%         ax1.YTick = -90:45:90;
% 
%         % text
%         text(0,1,sprintf('N = %d sessions',size(dwellmap,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','w')
% 
%         % plot yaw only
%         ax1c = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-ysiz2-ybuff2 ax1.Position(3) ysiz2]); 
%             dwellmap = posdata.hd_3d_dwell(:,1); 
%             v1 = cat(3,dwellmap{:});
%             vm = squeeze(sum(v1,1,'omitmissing'))';
%             vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
% 
%             x = linspace(-180,180,size(vm,2));
%             m = mean(vm,1,'omitmissing');
%             s = std(vm,[],1,'omitmissing');
% 
%             b = boundedline(x,m,s,'k','cmap',rgb('Black'));
% 
%             % axis settings
%             ax1c.XTick = -180:90:180;
%             ax1c.XLim = [-180 180];
%             xlabel(sprintf('Azimuth (%c)',176))
%             ax1c.YLim = [-2 3];
%             ylabel(sprintf('Time (z)'))
% 
%             % additional plots
%             line(ax1c.XLim,[0 0],'Color',[.5 .5 .5]);
%             line([0 0],ax1c.YLim,'Color','k');
%             line([90 90],ax1c.YLim,'Color','k','LineStyle',':');
%             line([-90 -90],ax1c.YLim,'Color','k','LineStyle',':');   
% 
%     % plot hills 3D HD map (averaged)
%     ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz ysiz]);
%         % get data
%         % dwellmap = posdata.hd_3d_dwell(:,2); % HD
%         dwellmap = posdata.mov_3d_dwell(:,2); % movement        
%         dwellmap = cat(3,dwellmap{:});
%         for ii = 1:size(dwellmap,3)
%             dwellmap(:,:,ii) = (dwellmap(:,:,ii) - mean(dwellmap(:,:,ii),'all','omitmissing')) ./ std(dwellmap(:,:,ii),[],'all','omitmissing');
%         end
%         dwellmap_a = mean(dwellmap,3,'omitmissing'); 
% 
%         % plot data
%         imagesc([-180 180],[-90 90],dwellmap_a,'alphadata',~isnan(dwellmap_a))
% 
%         % axis settings
%         % xlabel(sprintf('Azimuth (%c)',176))
%         axis xy tight
%         view(0,90);
%         ax2.FontSize = 8;
%         ax2.CLim = ax1.CLim;
%         colormap(gca,turbo)
%         ax2.XTick = [];
%         ax2.YTick = ax1.YTick;        
%         ax2.YTickLabel = {};
% 
%         % plot yaw only
%         ax2c = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)-ysiz2-ybuff2 ax2.Position(3) ysiz2]); 
%             dwellmap = posdata.hd_3d_dwell(:,2); 
%             v1 = cat(3,dwellmap{:});
%             vm = squeeze(sum(v1,1,'omitmissing'))';
%             vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
% 
%             x = linspace(-180,180,size(vm,2));
%             m = mean(vm,1,'omitmissing');
%             s = std(vm,[],1,'omitmissing');
% 
%             b = boundedline(x,m,s,'k','cmap',rgb('Black'));
% 
%             % axis settings
%             ax2c.XTick = -180:90:180;
%             ax2c.XLim = [-180 180];
%             xlabel(sprintf('Azimuth (%c)',176))
%             ax2c.YLim = ax1c.YLim;
%             ax2c.YTick = [];
%             % ylabel(sprintf('Time (z)'))
% 
%             % additional plots
%             line(ax2c.XLim,[0 0],'Color',[.5 .5 .5]);
%             line([0 0],ax2c.YLim,'Color','k');
%             line([90 90],ax2c.YLim,'Color','k','LineStyle',':');
%             line([-90 -90],ax2c.YLim,'Color','k','LineStyle',':'); 



% %% >>>>>>>>>> Behaviour anisotropy map
%     xnow = 330;
%     ynow = ynow+yplot_buff;
% 
%     var = 'anisotropy_map';
%     v1 = posdata.(var)(:,1); % arena 1 data
%     all_amaps_a1 = cat(3,v1{:});
%     a1 = mean(all_amaps_a1,3,'omitnan');
% 
%     v2 = posdata.(var)(:,2); % hills data (surficial)
%     v2 = cellfun(@(x) imresize(x,[48 116],'bilinear'),v2,'UniformOutput',false);
%     all_amaps_a2 = cat(3,v2{:});
%     a2 = mean(all_amaps_a2,3,'omitnan');
% 
%     v3 = posdata.(var)(:,3); % arena 2 data
%     all_amaps_a3 = cat(3,v3{:});
%     a3 = mean(all_amaps_a3,3,'omitnan');
% 
%     % plot arena
%     ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
%         ah = add_panel_title('e',sprintf('Avg. behavioural anisotropy'),'yoffset',-20,'xoffset',20,'width',400);  
% 
%         imagesc(a1,'alphadata',~isnan(a1)); hold on;
%         axis xy on
%         daspect([1 1 1])
%         colormap(ax1,cmocean('balance'))
%         % colormap(ax1,cmocean('thermal'))        
%         ax1.CLim = [-0.8 0.8];
%         ax1.XTick = [];
%         ax1.YTick = [];
% 
%         % text
%         text(0,1,sprintf('N = %d sessions',size(all_amaps_a1,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')
% 
%     % plot hills
%     ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
%         imagesc(a2,'alphadata',~isnan(a2)); hold on;
%         axis xy on
%         daspect([1 1 1])
%         colormap(ax2,ax1.Colormap)
%         ax2.CLim = ax1.CLim;
%         ax2.XTick = [];
%         ax2.YTick = [];
% 
%         % additional plots
%         ps = linspace(0,size(a2,2),7);
%         tops = ps(2:2:end);
%         plot([tops; tops],ax2.YLim,'k:')
% 
%         % text
%         text(0,1,sprintf('N = %d sessions',size(all_amaps_a2,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')
% 
%     % horizontal colormap
%     axc = axes('Units','pixels','Position',[xnow+210 ynow+80 70 8]);
%         x = linspace(ax1.CLim(1),ax1.CLim(2),100);
%         imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
%         colormap(axc,ax1.Colormap);
%         axis on
%         axc.XTick = [];
%         axc.YTick = [];
%         axc.FontSize = 7;
%         text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
%         text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
%         text(0.5,1,sprintf('Anisotropy score'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')


%% >>>>>>>>>> Anisotropy scatter density
    xnow = 50;
    ynow = ynow-165;

    % get data
    v1 = [anisotropy_map_1(:) a1(:)];
    v2 = [anisotropy_map_2(:) a2(:)];      
    bin_res = 256;

    xbuff = 15;
    xsiz = 125;
    ysiz = 125;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('g',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400);  
    
        % get distribution
        [F,ye,xe] = histcounts2(v1(:,2),v1(:,1),-1:0.1:1,-1:0.1:1,'Normalization','probability');
        F = ( imgaussfilt(F,1) );
    
        % plot data
        contour(movmean(xe,2,'EndPoints','discard'),movmean(ye,2,'EndPoints','discard'),F,'Color',[.5 .5 .5]); hold on;        
        msiz = 8; 
        s1 = scatter(v1(:,1),v1(:,2),msiz,plot_set{1,1},'filled','Marker','o','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none'); hold on;   

        s1f = polyfit(v1(:,1),v1(:,2),1);
        rf1 = refline(s1f(1),s1f(2));
        set(rf1,'Color',plot_set{1,1},'LineWidth',1);

        % axis settings
        axis xy
        box off
        ax1.XLim = [-0.5 0.5];
        ax1.YLim = [-0.8 0.8];        
        xlabel('Field anisotropy')
        ylabel('Behavioural anisotropy')       
        colormap(ax1,'turbo')
        ytickformat('%.1f')

        % text
        text(0.5,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        [r1,p1] = corr(v1(:,2),v1(:,1),'type','Pearson','rows','pairwise');     
        text(ax1,1,0,sprintf('{\\itr} = %.2f\n{\\itp} < .001',r1),'Units','normalized','FontSize',7,'Color','k','VerticalAlignment','bottom','HorizontalAlignment','right')

    % plot hills
    ax2 = axes('Units','pixels','Position',[xnow++ax1.Position(3)+xbuff ynow xsiz ysiz]);
        % get distribution
        [F,ye,xe] = histcounts2(v2(:,2),v2(:,1),-1:0.1:1,-1:0.1:1,'Normalization','probability');
        F = ( imgaussfilt(F,1) );
    
        % plot data
        contour(movmean(xe,2,'EndPoints','discard'),movmean(ye,2,'EndPoints','discard'),F,'Color',[.5 .5 .5]); hold on;        
        s1 = scatter(v2(:,1),v2(:,2),msiz,plot_set{1,2},'filled','Marker','o','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none'); hold on;      

        s1f = polyfit(v2(:,1),v2(:,2),1);
        rf1 = refline(s1f(1),s1f(2));
        set(rf1,'Color',plot_set{1,2},'LineWidth',1);

        % axis settings
        axis xy
        box off
        ax2.XLim = [-0.5 0.5];
        ax2.YLim = [-0.8 0.8];   
        ax2.YTick = [];
        xlabel('Field anisotropy')
        colormap(ax2,'turbo')
        ytickformat('%.1f')

        % text
        text(0.5,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        [r2,p2] = corr(v2(:,2),v2(:,1),'type','Pearson','rows','pairwise');     
        text(ax2,1,0,sprintf('{\\itr} = %.2f\n{\\itp} < .001',r2),'Units','normalized','FontSize',7,'Color','k','VerticalAlignment','bottom','HorizontalAlignment','right')

% return
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig 5.png']; 
        if fast_figs
            frame = getframe(gcf); % fig is the figure handle to save
            [raster, raster_map] = frame2im(frame); % raster is the rasterized image, raster_map is the colormap
            if isempty(raster_map)
                imwrite(raster, fname);
            else
                imwrite(raster, raster_map, fname); % fig_file is the path to the image
            end            
        else
            exportgraphics(gcf,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',350);  
        end
        close(gcf);   
    end







