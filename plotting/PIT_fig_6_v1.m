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
%% >>>>>>>>>> Behaviour anisotropy example trajectory
    xnow = 30;
    ynow = 600;
    ybuff = 80;
    xsiz = 130;
    ysiz = 80;
    xbuff = 140;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+ybuff xsiz ysiz]);
        ah = add_panel_title('b',sprintf('Behavioural anisotropy'),'yoffset',0,'xoffset',10,'width',400);  

        % get data
        eg_session = 2;
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
        set(c,'LineWidth',1);

        % axis settings
        ax1 = gca;
        axis xy off tight
        daspect([1 1 1])
        cm = cmocean('balance');
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        colormap(ax1,cmap)
        view(0,90)
        ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax1.YLim = [min(mf(:,2)) max(mf(:,2))];

        % text
        text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        text(0,1,sprintf('Example session'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        % get data
        pp = 2;
        idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
        pox = pos.pox_surficial(idx);
        poy = pos.poy_surficial(idx);
        poz = pos.poz_surficial(idx);
        poh = pos.yaw(idx);
        mf = posdata.maze_frame{eg_session};

        % plot data
        plot3(mf(:,4),mf(:,2),mf(:,3),'k'); hold on;        
        c = cline(pox,poy,poz,poh);
        set(c,'LineWidth',1);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        colormap(ax2,cmap)
        view(0,90)
        ax2.XLim = [min(mf(:,4)) max(mf(:,4))];
        ax2.YLim = [min(mf(:,2)) max(mf(:,2))];

        % additional plots
        ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

        % text
        text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        
    % colorbar
    axc = axes('Units','pixels','Position',[ax2.Position(1)+(ax2.Position(3).*0.125) ax2.Position(2)-10 ax2.Position(3).*0.75 10]);
        mat = (linspace(0,100,100));
        imagesc(ones(size(mat)),mat,mat);
        colormap(axc,ax2.Colormap);
        axis xy

        axc.YTick = [];
        axc.XTick = [];
        text(0.5,1.6,sprintf('Direction (%c)',176),'FontSize',7,'HorizontalAl','center','Units','normalized')
        text(0,0,sprintf('-180'),'FontSize',7,'HorizontalAl','left','Units','normalized','VerticalAl','top')
        text(0.5,0,sprintf('0'),'FontSize',7,'HorizontalAl','center','Units','normalized','VerticalAl','top')        
        text(1,0,sprintf('180'),'FontSize',7,'HorizontalAl','right','Units','normalized','VerticalAl','top')

%% >>>>>>>>>> Behaviour anisotropy map
    xnow = xnow+xbuff;
    
    var = 'anisotropy_map';
    v1 = posdata.(var)(:,1); % arena 1 data
    all_amaps_a1 = cat(3,v1{:});
    a1 = mean(all_amaps_a1,3,'omitnan');

    v2 = posdata.(var)(:,2); % hills data (surficial)
    v2 = cellfun(@(x) imresize(x,[48 116],'bilinear'),v2,'UniformOutput',false);
    all_amaps_a2 = cat(3,v2{:});
    a2 = mean(all_amaps_a2,3,'omitnan');

    v3 = posdata.(var)(:,3); % arena 2 data
    all_amaps_a3 = cat(3,v3{:});
    a3 = mean(all_amaps_a3,3,'omitnan');

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+ybuff xsiz ysiz]);    
        imagesc(a1,'alphadata',~isnan(a1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,cmocean('balance'))
        ax1.CLim = [-0.8 0.8];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,1,sprintf('All sessions'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

    % plot hills
    ax2 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        imagesc(a2,'alphadata',~isnan(a2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)
        ax2.CLim = ax1.CLim;
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(a2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

    % colorbar
    axc = axes('Units','pixels','Position',[ax2.Position(1)+(ax2.Position(3).*0.125) ax2.Position(2)-10 ax2.Position(3).*0.75 10]);
        mat = (linspace(0,100,100));
        imagesc(ones(size(mat)),mat,mat);
        colormap(axc,ax2.Colormap);
        axis xy

        axc.YTick = [];
        axc.XTick = [];
        text(0.5,1.6,sprintf('Behavioural anisotropy'),'FontSize',7,'HorizontalAl','center','Units','normalized')
        text(0,0,sprintf('%.1f',ax2.CLim(1)),'FontSize',7,'HorizontalAl','left','Units','normalized','VerticalAl','top')
        text(0.5,0,sprintf('0'),'FontSize',7,'HorizontalAl','center','Units','normalized','VerticalAl','top')        
        text(1,0,sprintf('%.1f',ax2.CLim(2)),'FontSize',7,'HorizontalAl','right','Units','normalized','VerticalAl','top')

%% >>>>>>>>>> Directional dwell time analysis (if necessary)
    dir_bins = 360;
    k = 0.01;
    ri = linspace(-pi,pi,dir_bins)';
    xi = movmean(ri,2,'Endpoints','discard');

    dmaps = NaN(size(posdata,1),numel(xi),3);
    for ss = 1:size(posdata,1) % for every session
        for pp = 1:3 % for every maze
            pos = posdata.pos{ss,1};
            mf = posdata.maze_frame{ss,1};
            
            d = wrapToPi(deg2rad(pos.yaw(pos.session==pp)));
            dmaps(ss,:,pp) = circ_ksdensity(d,xi,[],k);

        end
    end

    xnow = xnow+260;
    xsiz = 140;
    ysiz = 140;

    % plot arena
    ax1 = polaraxes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
        alph = 0.80;
        ri = ri(:)';
        f1 = squeeze(mean(dmaps(:,:,1),1,'omitnan'));
        f2 = squeeze(mean(dmaps(:,:,2),1,'omitnan'));

        p1 = polarhistogram('BinEdges',ri,'BinCounts',f1,'DisplayStyle','stairs','EdgeColor',plot_set{1,1},'FaceColor','none'); hold on;
        p2 = polarhistogram('BinEdges',ri,'BinCounts',f2,'DisplayStyle','stairs','EdgeColor',plot_set{1,2},'FaceColor','none'); hold on;

        % axis settings
        ax1.RTick = [];
        ax1.RAxisLocation = 180;
        ax1.ThetaZeroLocation = 'right';
        ax1.LineWidth = 1.5;

        % legend
        axt = axes('Units','pixels','Position',[xnow ynow 200 150],'Color','none');
            axt.XLim = [0 1];
            axt.YLim = [0 1];
            axis off
            p1 = patch(repmat(-10,4,1),repmat(-10,4,1),'k','FaceColor','none','EdgeColor',p1.EdgeColor); hold on;
            p2 = patch(repmat(-10,4,1),repmat(-10,4,1),'k','FaceColor','none','EdgeColor',p2.EdgeColor);            
            [~,leg] = legendflex([p1,p2],{'Arena 1','Hills'},'anchor',{'n','n'},'ncol',2,'box','off','buffer',[-35,35],'xscale',.5,'fontsize',8); 
            leg(3).FaceAlpha = p1.FaceAlpha;
            leg(4).FaceAlpha = p2.FaceAlpha;

% return


%% >>>>>>>>>> Rearing analysis (if necessary)
    rear_cutoff = 3;
    rear_duration_min = 0.2;

    all_rears = cell(size(posdata,1),3);
    dmaps = cell(size(posdata,1),3);
    for ss = 1:size(posdata,1) % for every session
        for pp = 1:3 % for every maze
            pos = posdata.pos{ss,1};
            mf = posdata.maze_frame{ss,1};
            
            idx = pos.session==pp;
            pox = pos.pox(idx) - min(mf(:,1));
            poy = pos.poy(idx) - min(mf(:,2));
            poz = -pos.poz(idx);
        
            mf(:,1) = mf(:,1)-min(mf(:,1));
            mf(:,2) = mf(:,2)-min(mf(:,2));

            if pp==2
                new_z = poz+pos.poz_curve(idx); % z-position minus the curvature of the maze
            else
                new_z = poz;
            end
            zscore_z = (new_z - mean(new_z,"all",'omitnan')) ./ std(new_z,[],"all",'omitnan');
            rears = zscore_z>rear_cutoff;
    
            b = bwlabel(rears);
            r = regionprops('table',b,zscore_z,'Area','Centroid','WeightedCentroid','MaxIntensity');
            rindx = r.Area > rear_duration_min/(1/50);
            r = r(rindx,:);
            if isempty(r)
                continue
            end
    
            r_points = round(r.WeightedCentroid(:,2));
            r_pos = [pox(r_points) poy(r_points) new_z(r_points) poz(r_points)];
            all_rears(ss,pp) = { r_pos };

            % density
            bs = 10;
            xi = min(mf(:,1)) : bs : max(mf(:,1));
            f = histcounts(r_pos(:,1),xi,'Normalization','count');
            dmaps(ss,pp) = { f };

% keyboard
% figure
% subplot(1,2,1)
%     cline(pox,poy,new_z,zscore_z);
%     daspect([1 1 1])
%     colorbar
% subplot(1,2,2)
%     r_points = round(r.WeightedCentroid(:,2));
%     scatter3(pox(r_points),poy(r_points),new_z(r_points),50,zscore_z(r_points),'filled')
%     daspect([1 1 1])
% keyboard
        end
    end

%% >>>>>>>>>> Behaviour anisotropy example trajectory
    xnow = [30 30+xbuff];
    ynow = 300;
    ybuff = 80;
    xsiz = 130;
    ysiz = 80;

    % plot arena
    for ii = 1:2
        ax1 = axes('Units','pixels','Position',[xnow(ii) ynow xsiz ysiz]);
            % ah = add_panel_title('b',sprintf('Behavioural anisotropy'),'yoffset',0,'xoffset',10,'width',400);  
            dat = all_rears(:,ii); % hills data
            idx = cell2mat(cellfun(@isempty,dat,'UniformOutput',false));
            dat = dat(~idx);
            dat = cat(1,dat{:});
    
            mf = posdata.maze_frame{1,1};
            mf(:,1) = mf(:,1)-min(mf(:,1));
            mf(:,2) = mf(:,2)-min(mf(:,2));
    
            dat = dat./1e3;
            mf = mf./1e3;
    
            plot(dat(:,1),dat(:,2),'ko','MarkerFaceColor','k','MarkerSize',2); hold on;
            plot(mf(:,1),mf(:,2),'k');
    
            daspect([1 1 1])
            axis xy off tight
            ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
            ax1.YLim = [min(mf(:,2)) max(mf(:,2))]; 
            ax1.XTickLabel = {};
    
            % additional plots
            ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax1.YLim,'r')        
    
        % plot arena
        ax2 = axes('Units','pixels','Position',[xnow(ii) ynow-30 xsiz 30]);
            bs = 50./1e3;
            xi = min(mf(:,1)) : bs : max(mf(:,1));
            f = histcounts(dat(:,1),xi,'Normalization','probability') .* 100;
    
            x = movmean(xi,2,'Endpoints','discard');
            bar(x,f,1,'k'); hold on;
    
            ax2.YScale = 'log';
            ax2.XLim = ax1.XLim;
            xlabel('Position (m)')
            ax2.YLim = [0.1 20];
            ax2.YTick
            if ii==2
                ax2.YTick = [];
            else
                ylabel('Rearing (%)')
            end

            % additional plots
            plot([tops; tops],ax2.YLim,'r')  

    end

return











%% >>>>>>>>>> Save the overall figure
    if 0
        fname = [config.fig_dir '\Fig 6.png']; 
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







