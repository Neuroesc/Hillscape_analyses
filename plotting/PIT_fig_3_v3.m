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
    warning('off','MATLAB:legend:IgnoringExtraEntries')

    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white
    fs = [15 10];

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> elongation (distribution)
    xnow = 50;
    ynow = 720;

    % collect data
    ucis = unique(clumaa.uci(pidx));
    dat = cell(1,3);
    prt = [1 2 2];
    for pp=1:3    
        datn = [];
        for uu = 1:length(ucis) % for each cell
            idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prt(pp);
            if pp==3
                fdata = clumaa.surficial_field_data{idx}; % field data for this cell in this session in this part
                rmap = clumaa.ratemap_surficial{idx};       
            else
                fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
                rmap = clumaa.ratemap_planar{idx};                       
            end
            if isempty(fdata)
                continue
            end
            maj_axis = fdata.MajorAxisLength(:,1);
            min_axis = fdata.MinorAxisLength(:,1);
            e = sqrt(1 - (min_axis ./ maj_axis));
            % datn = [datn; e fdata.MajorAxisLength(:,1) fdata.MinorAxisLength(:,1)];
            datn = [datn; mean(e,'all','omitmissing') mean(fdata.MajorAxisLength(:,1),'all','omitmissing') mean(fdata.MinorAxisLength(:,1),'all','omitmissing')];            
        end
        dat{pp} = datn;
    end

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 270 125]);
        ah = add_panel_title('A',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);  

        v1 = dat{1}(:,1); % arena 1
        v2 = dat{2}(:,1); % hills
        v3 = dat{3}(:,1); % hills surficial

        xi = 0:0.01:1;
        bw = 0.03;
        f1 = ksdensity(v1(:),xi(:),"Bandwidth",bw);
        f2 = ksdensity(v2(:),xi(:),"Bandwidth",bw);
        f3 = ksdensity(v3(:),xi(:),"Bandwidth",bw);

        % main plot
        alph = 0.7;  
        a1 = area(xi,f1,'FaceColor',plot_set{1,1},'EdgeColor','none','FaceAlpha',alph); hold on;
        a2 = area(xi,f2,'FaceColor',plot_set{1,2},'EdgeColor','none','FaceAlpha',alph);
        % a3 = area(xi,f3,'FaceColor','none','EdgeColor',plot_set{1,2});

        % axis settings
        ax.XLim = [0 1]; 
        box off
        xlabel(sprintf('Elongation'))    
        ylabel(sprintf('PDF'))    
        xtickformat('%.1f');
        ytickformat('%.1f');
        ax.FontSize = 8;

        % text
        text(0.05,0.5,sprintf('N = %d cells',numel(v1)),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',ax.FontSize,'Units','normalized','Color',plot_set{1,1})
        text(0.05,0.4,sprintf('N = %d cells',numel(v2)),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',ax.FontSize,'Units','normalized','Color',plot_set{1,2})

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 
        [~,leg] = legendflex([a1 a2],{maze_names{1},maze_names{2}},'anchor',{'nw','nw'},'ncol',1,'box','off','buffer',[10,-10],'xscale',.5,'fontsize',9); 
        leg(3).Children.FaceAlpha = a1.FaceAlpha;
        leg(4).Children.FaceAlpha = a2.FaceAlpha;

    % example ellipses
    ax2 = axes('Units','pixels','Position',[ax.Position(1) ynow-120 ax.Position(3) ax.Position(3)],'Color','none');
        ax2.XLim = ax.XLim;
        daspect([1 1 1])
        axis manual off
        ax2.Clipping = 'off';
    
        xvals = 0.9:-0.2:0.1;
        ye = 0;   
        cols = winter(length(xvals));
        ang = 90;
        for xx = 1:length(xvals)
            emin = 0.04;
            emax = emin/(1-xvals(xx)^2);
    
            e = ellipse(emax,emin,deg2rad(ang),xvals(xx),ye(1));
            x = e.XData;
            y = e.YData;
            y = y - max(y) + ax2.YLim(2);
            patch(x,y,'w','EdgeColor','k','FaceAlpha',1,'Clipping','off'); hold on;
            plot(mean(x(:)),mean(y(:)),'k+');
            delete(e);
            text(xvals(xx),max(y(:)),sprintf('%.1f',xvals(xx)),'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','bottom')
        end

%% >>>>>>>>>> orientation analysis results        
    xnow = xnow+320;
    ynow = ynow;

    ax = axes('Units','pixels','Position',[xnow ynow 230 125]);    
        ah = add_panel_title('B',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);  
    
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
    
                % get field anisotropy
                fanisotropy = (fdata.BoundingBox(:,4) - fdata.BoundingBox(:,3)) ./  (fdata.BoundingBox(:,4) + fdata.BoundingBox(:,3)); % height - width / height + width

                % get elongation/eccentricity
                e = sqrt(1 - (fdata.MinorAxisLength(:,1) ./ fdata.MajorAxisLength(:,1)));

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
                dat = [dat; fdata.MajorAxisLength(:,1) fdata.MinorAxisLength(:,1) fdata.Orientation(:,1) fdata.WeightedCentroid(:,1).*scale_x fdata.WeightedCentroid(:,2) fdata.Area(:,1) ds(:,1) repmat(uu,size(fdata,1),1) ds(:,2) ds(:,3:4) fanisotropy e];
            end
            datn(pp) = { dat };
        end

        % arena 1
        v1 = datn{1}(:,3);
        v2 = datn{2}(:,3);
        v3 = datn{3}(:,3);

        % get circular density of data
        ri = linspace(-90,90,360);
        xi = linspace(-180,180,360);        
        k = 0.25;
        f1 = circ_ksdensity(deg2rad(v1).*2,deg2rad(xi),[],k);
        f2 = circ_ksdensity(deg2rad(v2).*2,deg2rad(xi),[],k);
        f3 = circ_ksdensity(deg2rad(v3).*2,deg2rad(xi),[],k);

        % plot data
        a1 = area(ri,f1,'FaceColor',plot_set{1,1},'EdgeColor','none','FaceAlpha',alph); hold on;
        a2 = area(ri,f2,'FaceColor',plot_set{1,2},'EdgeColor','none','FaceAlpha',alph); hold on;

        % axis settings
        ax.XLim = [-90 90];
        ax.XTick = -180:45:180;
        box off
        xlabel(sprintf('Field orientation (%c)',176))
        ylabel('PDF')
        ytickformat('%.1f');
        ax.FontSize = 8;

        % text
        text(0.05,0.75,sprintf('N = %d fields',numel(v1)),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',ax.FontSize,'Units','normalized','Color',plot_set{1,1})
        text(0.05,0.65,sprintf('N = %d fields',numel(v2)),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',ax.FontSize,'Units','normalized','Color',plot_set{1,2})

        % stats
        rv1 = circ_r(deg2rad(v1).*2); % angle double
        rv2 = circ_r(deg2rad(v2).*2); % angle double
        % [~,leg] = legendflex([a1 a2],{sprintf('Arena 1: Rayleigh v = %.2f',rv1),sprintf('Hills: Rayleigh v = %.2f',rv2)},'anchor',{'sw','sw'},'ncol',1,'box','off','buffer',[0,-70],'xscale',0.5,'FontSize',8);   
        % leg(3).Children.FaceAlpha = alph;
        % leg(4).Children.FaceAlpha = alph;

    % example ellipses
    ax2 = axes('Units','pixels','Position',[ax.Position(1) ynow-95 ax.Position(3) ax.Position(3)],'Color','none');    
        ax2.XLim = ax.XLim;
        ax2.YLim = ax.XLim;
        daspect([1 1 1])
        axis manual off
        ax2.Clipping = 'off';
    
        xvals = -90:45:90;
        xvals = [-80 -45 0 45 80];
        ye = 90;   
        cols = winter(length(xvals));
        for xx = 1:length(xvals)
            emin = 8;
            emax = emin/(1-0.7^2);
    
            e = ellipse(emax,emin,deg2rad(xvals(xx)),xvals(xx),ye(1));
            x = e.XData;
            y = e.YData;
            patch(x,y,'w','EdgeColor','k','FaceAlpha',1,'Clipping','off'); hold on;
            plot(mean(x(:)),mean(y(:)),'k+');            
            delete(e);
            text(xvals(xx),max(y(:)),sprintf('%.f',xvals(xx)),'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','bottom')
        end

        % arrow([-93 70],[-53 70],'FaceColor',[.5 .5 .5],'Length',5,'EdgeColor',[.5 .5 .5]);
        % arrow([-93 70],[-93 110],'FaceColor',[.5 .5 .5],'Length',5,'EdgeColor',[.5 .5 .5]);
        % text(-93,110,'Y','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'Color',[.5 .5 .5])
        % text(-53,70,'X','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8,'Color',[.5 .5 .5])

%% >>>>>>>>>> All fields
    xnow = -5;
    ynow = ynow-325;
    xbuff = 210;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+50 220 90]);
        ah = add_panel_title('C',sprintf('Example fields'),'yoffset',-10,'xoffset',45,'width',400,'fontsize',fs);  
        % ellipses
        f = datn{1};
        alph = 0.3;
        nfields = 150;
        rng(99); % for reproducilibity
        rindx = randperm(size(f,1),nfields);
        for ff = 1:nfields
            e = ellipse(f(rindx(ff),1)./2,f(rindx(ff),2)./2,-deg2rad(f(rindx(ff),3)),f(rindx(ff),4),f(rindx(ff),5),'none');
            x = e.XData;
            y = e.YData;
            patch(x,y,(f(rindx(ff),12)),'EdgeColor','none','FaceAlpha',alph,'Clipping','off'); hold on;
            delete(e);
        end

        axis xy on
        box on
        daspect([1 1 1])
        % cmap_now = flipud(cmocean('thermal'));
        % cmap_now = flipud(cmocean('balance'));
        cmap_now = cmocean('balance');
        % cmap_now = cmocean('curl');  
        % cmap_now = inferno(128);
        % c = cmocean('balance'); 
        % c = flipud(c);
        % c = viridis(64);
        % colormap(ax1,[flipud(cmap_now);cmap_now])   
        % ax1.CLim = [-90 90];
        colormap(ax1,cmap_now)                   
        ax1.CLim = [-0.5 0.5];
        ax1.XLim = [0.5 size(posdata.anisotropy_map{1,1},2)+0.5];
        ax1.YLim = [0.5 size(posdata.anisotropy_map{1,1},1)+0.5];        
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,0,sprintf('%s, %d random fields',maze_names{1},nfields),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','top','rotation',0)

    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
        % ellipses
        f = datn{2};      
        for ff = 1:nfields
            e = ellipse(f(rindx(ff),1)./2,f(rindx(ff),2)./2,-deg2rad(f(rindx(ff),3)),f(rindx(ff),4),f(rindx(ff),5),'none');
            x = e.XData;
            y = e.YData;
            patch(x,y,(f(rindx(ff),12)),'EdgeColor','none','FaceAlpha',alph,'Clipping','off'); hold on;
            delete(e);
        end

        axis xy on
        box on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)   
        ax2.CLim = ax1.CLim;        
        ax2.XLim = [0.5 size(posdata.anisotropy_map{1,2},2)+0.5];
        ax2.YLim = [0.5 size(posdata.anisotropy_map{1,2},1)+0.5]; 
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(ax2.XLim(1),ax2.XLim(2),7);
        tops = ps(2:2:end);
        pt = plot([tops; tops],ax2.YLim,'k--');

        % text
        text(0,0,sprintf('%s, %d random fields',maze_names{2},nfields),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','top','rotation',0)

    % plot arena 2
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+xbuff ax2.Position(2) ax1.Position(3) ax1.Position(4)]);
        % ellipses
        f = datn{3};      
        for ff = 1:nfields
            e = ellipse(f(rindx(ff),1)./2,f(rindx(ff),2)./2,-deg2rad(f(rindx(ff),3)),f(rindx(ff),4),f(rindx(ff),5),'none');
            x = e.XData;
            y = e.YData;
            patch(x,y,(f(rindx(ff),12)),'EdgeColor','none','FaceAlpha',alph,'Clipping','off'); hold on;
            delete(e);
        end

        axis xy on
        box on
        daspect([1 1 1])
        colormap(ax3,ax1.Colormap)     
        ax3.CLim = ax1.CLim;
        ax3.XLim = [0.5 size(posdata.anisotropy_map{1,1},2)+0.5];
        ax3.YLim = [0.5 size(posdata.anisotropy_map{1,1},1)+0.5];        
        ax3.XTick = [];
        ax3.YTick = [];

        % text
        text(0,0,sprintf('%s, %d random fields',maze_names{3},nfields),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','top','rotation',0)

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+350 ynow+150 100 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 9;

        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(-0.3,1.5,sprintf('Anisotropy score:'),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')

    % example ellipses
    ax2e = axes('Units','pixels','Position',[axc.Position(1) axc.Position(2)+20 axc.Position(3) axc.Position(3)],'Color','none');
        ax2e.XLim = axc.XLim;
        ax2e.YLim = axc.XLim;
        daspect([1 1 1])
        axis manual off
        ax2e.Clipping = 'off';

        nf = 5;
        xvals = linspace(-0.5,0.5,nf); 
        ye = -0.5;
        alph = 0.5;
        scaler = 0.22;
        c = cmocean('balance',nf);
        for xx = 1:length(xvals)
            emin = 1;
            emax = -(emin + emin*xvals(xx))/(xvals(xx) - 1);
            v = sqrt(emin^2+emax^2);
            emin = emin/v*scaler;
            emax = emax/v*scaler;
            r = rectangle('Position',[xvals(xx)-(emin/2),ye(1)-(emax/2),emin,emax],'Curvature',[1 1],'EdgeColor','k','Clipping','off','FaceColor',[c(xx,:) 0.5]); hold on;
        end
        colormap(ax2e,c);

%% >>>>>>>>>> Field anisotropy map plot
    xnow = -5;
    ynow = ynow-145;
    xbuff = 210;

    % map data
    % arena 1
    f = datn{1};
    dist_cutoff = 160;
    mapset.binsize = 32; % (mm) firing rate map bin size
    
    dist_cutoff_bins = dist_cutoff ./ mapset.binsize;    
    anisotropy_map_1 = NaN(48,96);
    mapXY = f(:,4:5);
    for bb = 1:numel(anisotropy_map_1)
        [r,c] = ind2sub(size(anisotropy_map_1),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_1(bb) = mean(f(box_idx,12),1,'omitnan');
    end
    anisotropy_map_1 = imgaussfilt(anisotropy_map_1,1);

    % hills  
    f = datn{2};    
    anisotropy_map_2 = NaN(48,116);
    mapXY = f(:,4:5);
    for bb = 1:numel(anisotropy_map_2)
        [r,c] = ind2sub(size(anisotropy_map_2),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_2(bb) = mean(f(box_idx,12),1,'omitnan');
    end
    anisotropy_map_2 = imgaussfilt(anisotropy_map_2,1);

    % arena 2  
    f = datn{3};    
    anisotropy_map_3 = NaN(48,96);
    mapXY = f(:,4:5);
    for bb = 1:numel(anisotropy_map_3)
        [r,c] = ind2sub(size(anisotropy_map_3),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_3(bb) = mean(f(box_idx,12),1,'omitnan');
    end
    anisotropy_map_3 = imgaussfilt(anisotropy_map_3,1);

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+50 220 90]);
        ah = add_panel_title('D',sprintf('Average field anisotropy'),'yoffset',-10,'xoffset',45,'width',400,'fontsize',fs);     
    
        imagesc(anisotropy_map_1,'alphadata',~isnan(anisotropy_map_1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,cmap_now)                
        ax1.CLim = [-0.5 0.5];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,0,sprintf('%s, %d fields',maze_names{1},size(datn{1},1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
        imagesc(anisotropy_map_2,'alphadata',~isnan(anisotropy_map_2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)
        ax2.CLim = ax1.CLim;
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(anisotropy_map_2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k--')

        % text
        text(0,0,sprintf('%s, %d fields',maze_names{2},size(datn{2},1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')

    % plot arena 2
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+xbuff ax2.Position(2) ax1.Position(3) ax1.Position(4)]);
        imagesc(anisotropy_map_3,'alphadata',~isnan(anisotropy_map_3)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax3,ax1.Colormap)
        ax3.CLim = ax1.CLim;
        ax3.XTick = [];
        ax3.YTick = [];

        % text
        text(0,0,sprintf('%s, %d fields',maze_names{3},size(datn{3},1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')

%% >>>>>>>>>> Behaviour anisotropy example trajectory
    xnow = 20;
    ynow = ynow-85;
    xsiz = 125;
    ysiz = 70;
    xbuff = 10;
    xbuff2 = 310;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('E',sprintf('Example trajectory'),'yoffset',-10,'xoffset',20,'width',400,'fontsize',fs);  

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
        lw = 1;
        set(c,'LineWidth',lw);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        cmap = [cmap_now; flipud(cmap_now);cmap_now; flipud(cmap_now)];
        colormap(ax1,cmap)
        view(0,90)
        ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax1.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax1.CLim = [-180 180];

        % text
        text(0,0,sprintf('%s, example session',maze_names{1}),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
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

        % text
        text(0,0,sprintf('%s, example session',maze_names{2}),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+200 ynow+75 70 8]);
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

%% >>>>>>>>>> Behaviour anisotropy map
    xnow = xnow+xbuff2;
    
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
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
        ah = add_panel_title('F',sprintf('Elongation is not explained\nby behaviour'),'yoffset',0,'xoffset',20,'width',400,'fontsize',fs); 
    
        imagesc(a1,'alphadata',~isnan(a1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,cmap_now)        
        ax1.CLim = [-0.8 0.8];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,0,sprintf('%s, %d sessions',maze_names{1},size(all_amaps_a1,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')

    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
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
        plot([tops; tops],ax2.YLim,'k--')

        % text
        text(0,0,sprintf('%s, %d sessions',maze_names{2},size(all_amaps_a2,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+200 ynow+75 70 8]);
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

%% >>>>>>>>>> Field to wall angle, expected maps
    % collect data
    ucis = unique(clumaa.uci(pidx));
    datn = cell(1,3);
    for pp=1:3    
        idx = ismember(clumaa.uci,ucis{1}) & clumaa.partn==pp;
        if pp==2
            rmap = clumaa.ratemap_surficial{idx};       
        else
            rmap = clumaa.ratemap_planar{idx};                       
        end

        b = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0]+0.5; 
        a = [90 0 90 0];
        gpoly = [];
        res = 1000;
        for jj = 1:size(b,1)-1
            gpoly = [gpoly; linspace(b(jj,1),b(jj+1,1),res)' linspace(b(jj,2),b(jj+1,2),res)' repmat(a(jj),res,1)];
        end
        gpoly = unique(round(gpoly),'rows');

        ps = linspace(0,size(rmap,2),7);
        tops = ps(2:2:end);
        bottoms = ps(3:2:end-2);
        b2 = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0;tops(1) 0;tops(1) size(rmap,1);tops(2) size(rmap,1);tops(2) 0;tops(3) 0;tops(3) size(rmap,1);tops(3) 0;bottoms(1) 0;bottoms(1) size(rmap,1);bottoms(2) size(rmap,1);bottoms(2) 0]+0.5; 
        a2n = [90 0 90 0 0 90 0 90 0 90 90 0 90 0 90];
        gpoly2 = [];
        res = 1000;
        for jj = 1:size(b2,1)-1
            gpoly2 = [gpoly2; linspace(b2(jj,1),b2(jj+1,1),res)' linspace(b2(jj,2),b2(jj+1,2),res)' repmat(a2n(jj),res,1)];
        end
        gpoly2 = unique(round(gpoly2),'rows');

        dmap = NaN([size(rmap),3,2]);

        angle_sigma = 10; % bigger = distal walls have more of an effect
        dist_sigma = 16; % bigger = walls have more distal effects on anisotropy
        for jj = 1:numel(dmap(:,:,1))
            %% arena
            % current coordinate
            [y,x] = ind2sub(size(dmap),jj);

            % distance to all walls
            ds = pdist2([x,y],gpoly(:,[1:2]),'euclidean');

            % wall angles
            ws = gpoly(:,3);
            weights = normpdf(ds,0,angle_sigma); % Gaussian weight wall angles
            mean_angle = sum(weights(:) .* ws(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');
            mean_distance = sum(weights(:) .* ds(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');

            dmap(y,x,1,1) = (mean_angle./45)-1;
            dmap(y,x,2,1) = normpdf(mean_distance,0,dist_sigma);
            dmap(y,x,3,1) = dmap(y,x,1,1) .* dmap(y,x,2,1);

            %% hills
            % current coordinate
            [y,x] = ind2sub(size(dmap),jj);

            % distance to all walls
            ds = pdist2([x,y],gpoly2(:,[1:2]),'euclidean');

            % wall angles
            ws = gpoly2(:,3);
            weights = normpdf(ds,0,angle_sigma); % Gaussian weight wall angles
            mean_angle = sum(weights(:) .* ws(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');
            mean_distance = sum(weights(:) .* ds(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');

            dmap(y,x,1,2) = (mean_angle./45)-1;
            dmap(y,x,2,2) = normpdf(mean_distance,0,dist_sigma);
            dmap(y,x,3,2) = dmap(y,x,1,2) .* dmap(y,x,2,2);
        end
        datn(pp) = { dmap };
    end

%% >>>>>>>>>> Expected anisotropy map (walls only)
    xnow = 20;
    ynow = ynow-130;
    xsiz = 125;
    ysiz = 70;
    xbuff = 10;
    xbuff2 = 310;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
        ah = add_panel_title('G',sprintf('Elongation is not explained\nby geometry alone'),'yoffset',0,'xoffset',20,'width',400,'fontsize',fs); 
        % get data
        dmap = datn{1};
        mnow1a = dmap(:,:,3,1);

        % plot data
        imagesc(mnow1a,'alphadata',~isnan(mnow1a)); hold on;

        % axis settings
        axis xy on
        daspect([1 1 1])
        colormap(ax1,cmap_now)
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,0,sprintf('%s, simulation',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','top','rotation',0)

    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);    
        % get data
        dmap = datn{2};
        mnow1b = dmap(:,:,3,1);

        % plot data
        imagesc(mnow1b,'alphadata',~isnan(mnow1b)); hold on;

        % axis settings
        axis xy on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(mnow1b,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k--')

        % text
        text(0,0,sprintf('%s, simulation',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','top','rotation',0)

%% >>>>>>>>>> Behaviour anisotropy map
    xnow = xnow+xbuff2;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
        ah = add_panel_title('H',sprintf('Elongation is explained well\nby geometry & terrain'),'yoffset',0,'xoffset',20,'width',400,'fontsize',fs); 
        % get data
        dmap = datn{1};
        mnow2a = dmap(:,:,3,1);

        % plot data
        imagesc(mnow2a,'alphadata',~isnan(mnow2a)); hold on;

        % axis settings
        axis xy on
        daspect([1 1 1])
        colormap(ax1,cmap_now)
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,0,sprintf('%s, simulation',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','top','rotation',0)
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);    
        % get data
        dmap = datn{2};
        mnow2b = dmap(:,:,3,2);

        % plot data
        imagesc(mnow2b,'alphadata',~isnan(mnow2b)); hold on;

        % axis settings
        axis xy on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(mnow2b,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k--')

        % text
        text(0,0,sprintf('%s, simulation',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','top','rotation',0)

    % % compare each prediction to the actual data
    % [r1a,p1a] = corr(anisotropy_map_1(:),a1(:),'type','pearson','rows','pairwise'); % arena 1 field anisotropy vs arena 1 behaviour anisotropy
    % [r1b,p1b] = corr(anisotropy_map_1(:),mnow1a(:),'type','pearson','rows','pairwise'); % arena 1 field anisotropy vs arena 1 anisotropy (walls only)
    % [r1c,p1c] = corr(anisotropy_map_1(:),mnow2a(:),'type','pearson','rows','pairwise'); % arena 1 field anisotropy vs arena 1 anisotropy (walls only)
    % 
    % [r2a,p2a] = corr(anisotropy_map_2(:),a2(:),'type','pearson','rows','pairwise'); % hills field anisotropy vs hills behaviour anisotropy
    % [r2b,p2b] = corr(anisotropy_map_2(:),mnow1b(:),'type','pearson','rows','pairwise'); % hills field anisotropy vs hills anisotropy (walls only)
    % [r2c,p2c] = corr(anisotropy_map_2(:),mnow2b(:),'type','pearson','rows','pairwise'); % hills field anisotropy vs hills anisotropy (walls only)

keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig 3.png']; 
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





