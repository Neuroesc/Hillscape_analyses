%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  analysis function written for:
% Grieves, Duvelle and Taube (202X) 
%
% USAGE:
%       [out] = template(in) process with default settings
% 
%       [out] = template(in,optional1) process using optional argument 1
% 
%       [out] = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'param1'          -   (default = X) Scalar value, parameter to do something
% 
%       'param2'          -   (default = X) Scalar value, parameter to do something
% 
% INPUT:
%       in    - input as a vector
%
% OUTPUT:
%       out   - output as a vector
%
% EXAMPLES:
%       % run function using default values
%       out = template(in,varargin)
%
% See also: GIT_audit

% HISTORY:
% version 1.0.0, Release 16/02/23 Code conception
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white
    fsiz = 9;
    fs = [15 10];

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
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


%% >>>>>>>>>> Area vs distance from closest boundary
    xnow = 50;
    ynow = 700;
    k = 1.5;

    % create axis
    ax1 = axes('Units','pixels','Position',[xnow ynow 125 125]);
        ah = add_panel_title('A',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % % field-to-wall angle function
        % fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
        
        % get data
        pp = 1; % arena 1
        f1 = datn{pp};
        field_area = f1(:,6) .* 0.0001;
        wall_distance = f1(:,7) .* 32 ./ 1000; % distance to closest wall

        % plot data
        d = computeScatterDensity(wall_distance,field_area,'k',k);
        scatter(wall_distance,field_area,20,d,'filled','MarkerEdgeColor','none'); hold on;
        r = refline;
        set(r,'Color','k')

        % axis settings
        grid on;
        colormap(ax1,'turbo')
        ax1.CLim(1) = 0;
        ax1.YLim = [0 0.8];
        ax1.XLim = [0 0.75];
        ax1.YTick = 0:0.2:0.8;
        ax1.XTick = 0:0.25:0.75;        
        ylabel(sprintf('Field area (m^{2})'))
        xlabel(sprintf('Field-to-wall distance (m)'))
        % ytickformat('%.e');
        % ax1.YScale = 'log';

        % text
        text(0.1,1.02,sprintf('%s',maze_names{1}),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz)
        % r = corr(wall_distance,ftw,'rows','pairwise','type','Spearman');
        % text(0.1,1,sprintf('{\\itr} = %.2f',r),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',9)
        p = polyfit(wall_distance,field_area,1);
        text(0.1,1,sprintf('y = %.2fx + %.2f',p(1),p(2)),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','top','FontSize',fsiz);

    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax1.Position(1)+20 ax1.Position(2)-60 80 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('0 '),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' Max'),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Density'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

    ax1b = axes('Units','pixels','Position',[ax1.Position(1)+170 ynow 80 125]);
        cut_dist = 0.2;
    
        % get data
        f1 = datn{1}; % arena 1
        field_area1 = f1(:,6) .* 0.0001;
        wall_distance1 = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        idx1 = wall_distance1 < cut_dist; % fields close to wall
    
        f1 = datn{2}; % hills
        field_area2 = f1(:,6) .* 0.0001;
        wall_distance2 = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        idx2 = wall_distance2 < cut_dist; % fields close to wall
    
        % vectorise
        [ds,gs] = vectorDATAGROUP([],field_area1(idx1),field_area1(~idx1)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[1 1]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        ylabel(sprintf('Field area (m^{2})'))
        ax1b.XTick = 1:2;
        ax1b.YTick = ax1.YTick;
        ax1b.YLim = ax1.YLim;        
        ax1b.XLim = [0.5 2.5]; 
        ax1b.XTickLabel = {'Wall','Centre'};
        ax1b.FontSize = 8;
        % ytickformat('%.1e');
        % ax1b.YScale = 'log';

        % stats
        axt = axes('Units','pixels','Position',ax1b.Position,'Color','none');
            axis off
            axt.XLim = ax1b.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.05,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);
    
    % Hills
    % create axis
    ax2 = axes('Units','pixels','Position',[ax1b.Position(1)+140 ynow 125 125]);
        ah = add_panel_title('B',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % % field-to-wall angle function
        % fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
        
        % get data
        pp = 2; % arena 1
        f1 = datn{pp};
        field_area = f1(:,6) .* 0.0001;
        wall_distance = f1(:,7) .* 32 ./ 1000; % distance to closest wall

        % plot data
        d = computeScatterDensity(wall_distance,field_area,'k',k);
        scatter(wall_distance,field_area,20,d,'filled','MarkerEdgeColor','none'); hold on;
        r = refline;
        set(r,'Color','k')

        % axis settings
        grid on
        colormap(ax2,'turbo')
        ax2.CLim(1) = 0;
        ax2.YLim = ax1.YLim;
        ax2.XLim = [0 0.75];
        ax2.YTick = ax1.YTick;    
        ax2.XTick = 0:0.25:0.75;
        ylabel(sprintf('Field area (m^{2})'))
        xlabel(sprintf('Field-to-wall distance (m)'))
        % ytickformat('%.1e');
        % ax2.YScale = 'log';
        p = polyfit(wall_distance,field_area,1);
        text(0.1,1,sprintf('y = %.2fx + %.2f',p(1),p(2)),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','top','FontSize',fsiz);

        % text
        text(0.1,1.02,sprintf('%s',maze_names{2}),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz)

    ax2b = axes('Units','pixels','Position',[ax2.Position(1)+170 ax1b.Position(2) 80 125]);            
        % vectorise
        [ds,gs] = vectorDATAGROUP([],field_area2(idx2),field_area2(~idx2)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[2 2]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        ylabel(sprintf('Field area (m^{2})'))
        ax2b.XTick = 1:2;
        ax2b.YTick = ax2.YTick;
        ax2b.YLim = ax2.YLim; 
        ax2b.XLim = [0.5 2.5]; 
        ax2b.XTickLabel = {'Wall','Centre'};
        ax2b.FontSize = 8;
        % ytickformat('%.1e');
        % ax2b.YScale = 'log';

        % stats
        axt = axes('Units','pixels','Position',ax2b.Position,'Color','none');
            axis off
            axt.XLim = ax2b.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.05,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);


%% >>>>>>>>>> Eccentricity vs distance from closest boundary
    xnow = 50;
    ynow = ynow-210;

    % create axis
    ax1 = axes('Units','pixels','Position',[xnow ynow 125 125]);
        ah = add_panel_title('C',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % field-to-wall angle function
        fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
        
        % get data
        pp = 1; % arena 1
        f1 = datn{pp};
        field_elongation = f1(:,13);
        wall_distance = f1(:,7) .* 32 ./ 1000; % distance to closest wall

        % plot data
        d = computeScatterDensity(wall_distance,field_elongation,'k',k);
        scatter(wall_distance,field_elongation,20,d,'filled','MarkerEdgeColor','none'); hold on;
        r = refline;
        set(r,'Color','k')

        % axis settings
        grid on;
        colormap(ax1,'turbo')
        ax1.CLim(1) = 0;
        ax1.YLim = [0 1];
        ax1.XLim = [0 0.75];
        ax1.YTick = 0:0.25:1;
        ax1.XTick = 0:0.25:0.75;        
        ylabel(sprintf('Field elongation'))
        xlabel(sprintf('Field-to-wall distance (m)'))
        ytickformat('%.2f')
        p = polyfit(wall_distance,field_elongation,1);
        text(0.1,0,sprintf('y = %.2fx + %.2f',p(1),p(2)),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz);

        % text
        text(0.1,1.02,sprintf('%s',maze_names{1}),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz)
        % r = corr(wall_distance,ftw,'rows','pairwise','type','Spearman');
        % text(0.1,1,sprintf('{\\itr} = %.2f',r),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',9)

    ax1 = axes('Units','pixels','Position',[ax1.Position(1)+170 ynow 80 125]);
        cut_dist = 0.2;
    
        % get data
        f1 = datn{1}; % arena 1
        field_elongation1 = f1(:,13);
        wall_distance1 = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        idx1 = wall_distance1 < cut_dist; % fields close to wall
    
        f1 = datn{2}; % hills
        field_elongation2 = f1(:,13);
        wall_distance2 = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        idx2 = wall_distance2 < cut_dist; % fields close to wall
    
        % vectorise
        [ds,gs] = vectorDATAGROUP([],field_elongation1(idx1),field_elongation1(~idx1)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[1 1]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        ylabel(sprintf('Field elongation'))
        ax1.XTick = 1:2;
        ax1.YTick = 0:0.25:1;
        ax1.XLim = [0.5 2.5]; 
        ax1.XTickLabel = {'Wall','Centre'};
        ax1.FontSize = 8;
        ytickformat('%.1f');
    
        % stats
        axt = axes('Units','pixels','Position',ax1.Position,'Color','none');
            axis off
            axt.XLim = ax1.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.05,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);
    
    % Hills
    % create axis
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+140 ynow 125 125]);
        ah = add_panel_title('D',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % field-to-wall angle function
        fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
        
        % get data
        pp = 2; % arena 1
        f1 = datn{pp};
        field_elongation = f1(:,13);
        wall_distance = f1(:,7) .* 32 ./ 1000; % distance to closest wall

        % plot data
        d = computeScatterDensity(wall_distance,field_elongation,'k',k);
        scatter(wall_distance,field_elongation,20,d,'filled','MarkerEdgeColor','none'); hold on;
        r = refline;
        set(r,'Color','k')

        % axis settings
        grid on
        colormap(ax2,'turbo')
        ax2.CLim(1) = 0;
        ax2.YLim = [0 1];
        ax2.XLim = [0 0.75];
        ax2.YTick = 0:0.25:1;    
        ax2.XTick = 0:0.25:0.75;
        ylabel(sprintf('Field elongation'))
        xlabel(sprintf('Field-to-wall distance (m)'))
        ytickformat('%.2f')
        p = polyfit(wall_distance,field_elongation,1);
        text(0.1,0,sprintf('y = %.2fx + %.2f',p(1),p(2)),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz);

        % text
        text(0.1,1.02,sprintf('%s',maze_names{2}),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz)

    ax2 = axes('Units','pixels','Position',[ax2.Position(1)+170 ax1.Position(2) 80 125]);            
        % vectorise
        [ds,gs] = vectorDATAGROUP([],field_elongation2(idx2),field_elongation2(~idx2)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[2 2]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        ylabel(sprintf('Field elongation'))        
        ax2.XTick = 1:2;
        ax2.YTick = 0:0.25:1;        
        ax2.XLim = [0.5 2.5]; 
        ax2.XTickLabel = {'Wall','Centre'};
        ax2.FontSize = 8;
        ytickformat('%.1f');

        % stats
        axt = axes('Units','pixels','Position',ax2.Position,'Color','none');
            axis off
            axt.XLim = ax2.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.05,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);


%% >>>>>>>>>> Firing rate vs distance from closest boundary
    xnow = 50;
    ynow = ynow-210;

    % create axis
    ax1 = axes('Units','pixels','Position',[xnow ynow 125 125]);
        ah = add_panel_title('E',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % arena 1
        m1 = clumaa.ratemap_planar(pidx & clumaa.partn==1);
        m1 = cat(3,m1{:});

        dmap = zeros(size(m1(:,:,1)));
        dmap = padarray(dmap,[1 1],1,'both');
        dmap = bwdist(dmap);
        dmap = dmap(2:end-1,2:end-1);
        dmap = dmap .* 32 ./ 1000; % in m

        % linear
        edg = linspace(0,max(dmap(:)),9);   
        [~,~,bidx] = histcounts(dmap(:),edg);                
        fun_mean = @(x) mean(x(:),"all",'omitmissing');
        fun_sem = @(x) std(x(:),[],"all",'omitmissing') ./ sqrt(sum(~isnan(x(:)),'all'));       
        dat = NaN(length(edg)-1,size(m1,3));
        for jj = 1:size(m1,3)
            map_now = reshape(m1(:,:,jj),[],1);
            dat(:,jj) = accumarray(bidx(bidx>0),map_now(bidx>0),[],fun_mean);
        end
        a1 = sum(m1>1,3,'omitmissing') ./ sum(~isnan(m1),3,'omitmissing');
        ca = accumarray(bidx(bidx>0),a1(bidx>0),[],fun_mean);
        ma = mean(dat,2,'omitmissing');
        sa = nansem(dat,2);

            % centre vs walls
            edg2 = [0,cut_dist,max(dmap(:))];   
            [~,~,bidx] = histcounts(dmap(:),edg2);  
            dat_a = NaN(length(edg2)-1,size(m1,3));            
            for jj = 1:size(m1,3)
                map_now = reshape(m1(:,:,jj),[],1);
                dat_a(:,jj) = accumarray(bidx(bidx>0),map_now(bidx>0),[],fun_mean);
            end
            cab = accumarray(bidx(bidx>0),a1(bidx>0),[],fun_mean);

        % ridge
        m2 = clumaa.ratemap_planar(pidx & clumaa.partn==2);
        m2 = cat(3,m2{:});

        dmap = zeros(size(m2(:,:,1)));
        dmap = padarray(dmap,[1 1],1,'both');
        dmap = bwdist(dmap);
        dmap = dmap(2:end-1,2:end-1);
        dmap = dmap .* 32 ./ 1000; % in m

        [~,~,bidx] = histcounts(dmap(:),edg);                
        dat = NaN(length(edg)-1,size(m1,3));
        for jj = 1:size(m2,3)
            map_now = reshape(m2(:,:,jj),[],1);
            dat(:,jj) = accumarray(bidx(bidx>0),map_now(bidx>0),[],fun_mean);
        end
        a2 = sum(m2>1,3,'omitmissing') ./ sum(~isnan(m2),3,'omitmissing');
        cr = accumarray(bidx(bidx>0),a2(bidx>0),[],fun_mean);
        mr = mean(dat,2,'omitmissing');
        sr = nansem(dat,2);

            % centre vs walls
            [~,~,bidx] = histcounts(dmap(:),edg2);   
            dat_r = NaN(length(edg2)-1,size(m1,3));            
            for jj = 1:size(m2,3)
                map_now = reshape(m2(:,:,jj),[],1);
                dat_r(:,jj) = accumarray(bidx(bidx>0),map_now(bidx>0),[],fun_mean);
            end
            crb = accumarray(bidx(bidx>0),a2(bidx>0),[],fun_mean);

        % plot data
        xi = movmean(edg,2,'Endpoints','discard');
        errorbar(xi,ma,sa,'Color',plot_set{1,1},'Marker',plot_set{2,1}); hold on;
        errorbar(xi,mr,sr,'Color',plot_set{1,2},'Marker',plot_set{2,2});

        % axis settings
        grid off
        box off
        ax1.YLim = [0 1.21];        
        ax1.XLim = [0 0.75];
        ax1.YTick = 0:0.2:10;    
        ax1.XTick = 0:0.25:0.75;
        ylabel(sprintf('Mean firing rate (Hz)'))
        xlabel(sprintf('Field-to-wall distance (m)'))
        ytickformat('%.1f')

        % stats
        % [r,p] = corr(fs(:,1),fs(:,2),'rows','pairwise','type','Pearson');
        text(0.1,.12,sprintf('Arena 1'),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color',plot_set{1,1})
        % [r,p] = corr(fs(:,1),fs(:,3),'rows','pairwise','type','Pearson');    
        text(0.1,.02,sprintf('Ridge'),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color',plot_set{1,2})

    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+170 ax1.Position(2) 80 125]);            
        % plot data
        m = mean(dat_a,2,'omitmissing');
        e = nansem(dat_a,2);
        errorbar(1:2,m,e,'Color',plot_set{1,1},'Marker',plot_set{2,1}); hold on;
        m = mean(dat_r,2,'omitmissing');
        e = nansem(dat_r,2);        
        errorbar(1:2,m,e,'Color',plot_set{1,2},'Marker',plot_set{2,2});

        % axis settings
        ylabel(sprintf('Mean firing rate (Hz)'))
        ax2.YLim = ax1.YLim;
        ax2.XTick = 1:2;
        ax2.YTick = ax1.YTick;        
        ax2.XLim = [0.5 2.5]; 
        ax2.XTickLabel = {'Wall','Centre'};
        ax2.FontSize = 8;
        ytickformat('%.1f');
        grid off
        box off

        % stats
        d1 = dat_a(:);
        d2 = dat_r(:);
        d = [d1(:); d2(:)];
        g1 = [ones(size(d1)); ones(size(d2)).*2]; % maze group
        g2 = [reshape(ones(size(dat_a)).*[1;2],[],1); reshape(ones(size(dat_r)).*[1;2],[],1)];

        [P,T,STATS,TERMS] = anovan(d,[g1 g2],'display','off','varnames',{'maze','distance'},'model','full');
%   Source          Sum Sq.   d.f.   Mean Sq.     F     Prob>F
% ------------------------------------------------------------
%   maze              3.404      1    3.4038     6.81   0.0092
%   distance         11.507      1   11.5065    23.01   0     
%   maze:distance     1.03       1    1.0296     2.06   0.1515
%   Error           902.103   1804    0.5001                  
%   Total           917.776   1807                                                        
        % axt = axes('Units','pixels','Position',ax2.Position,'Color','none');
        %     axis off
        %     axt.XLim = ax2.XLim;
        %     axt.YLim = [0 1];        
        %     [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.05,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

%% >>>>>>>>>> Co-active cells vs distance from closest boundary
    xnow = xnow+310;
    ynow = ynow;
    ax3 = axes('Units','pixels','Position',[xnow ynow 125 125]);
        ah = add_panel_title('F',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % plot data
        xi = movmean(edg,2,'Endpoints','discard');
        plot(xi,ca.*100,'Color',plot_set{1,1},'Marker',plot_set{2,1}); hold on;
        plot(xi,cr.*100,'Color',plot_set{1,2},'Marker',plot_set{2,2});

        % axis settings
        grid off
        box off
        ax3.YLim = [0 40];
        ax3.XLim = [0 0.75];
        ax3.YTick = 0:10:100;    
        ax3.XTick = 0:0.25:0.75;
        ylabel(sprintf('Coactive place cells (%%)'))
        xlabel(sprintf('Field-to-wall distance (m)'))
        ytickformat('%d')

        % stats
        % [r,p] = corr(fs(:,1),fs(:,2),'rows','pairwise','type','Pearson');
        % text(0,1.12,sprintf('Arena 1'),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color',plot_set{1,1})
        % [r,p] = corr(fs(:,1),fs(:,3),'rows','pairwise','type','Pearson');    
        % text(0,1.02,sprintf('Ridge'),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color',plot_set{1,2})

    ax4 = axes('Units','pixels','Position',[ax3.Position(1)+170 ax3.Position(2) 80 125]);            
        % plot data
        plot(1:2,cab.*100,'Color',plot_set{1,1},'Marker',plot_set{2,1}); hold on;       
        plot(1:2,crb.*100,'Color',plot_set{1,2},'Marker',plot_set{2,2});

        % axis settings
        ylabel(sprintf('Coactive place cells (%%)'))
        ax4.YLim = ax3.YLim;
        ax4.XTick = 1:2;
        ax4.YTick = ax3.YTick;        
        ax4.XLim = [0.5 2.5]; 
        ax4.XTickLabel = {'Wall','Centre'};
        ax4.FontSize = 8;
        ytickformat('%d');
        grid off
        box off
% walls comparison
    % Difference	0.78 %
    % 95% CI	-4.1858% to 5.7912%
    % Chi-squared	0.094
    % DF 	1
    % Significance level	P = 0.7588
    % https://www.medcalc.org/calc/comparison_of_proportions.php
    % text(1,max([cab(1) crb(1)]).*1.1,sprintf('{\\itp} = 0.76'),'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','bottom');       
% center comparison
    % Difference	2.90 %
    % 95% CI	-2.4991% to 8.3216%
    % Chi-squared	1.104
    % DF 	1
    % Significance level	P = 0.2933
    % https://www.medcalc.org/calc/comparison_of_proportions.php
    % text(2,max([cab(2) crb(2)]).*1.1,sprintf('{\\itp} = 0.29'),'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','bottom');
% arena 1 wall vs center
    % Difference	5.44 %
    % 95% CI	0.0293% to 10.8179%
    % Chi-squared	3.886
    % DF 	1
    % Significance level	P = 0.0487
    % https://www.medcalc.org/calc/comparison_of_proportions.php    
% hills wall vs center
    % Difference	3.32 %
    % 95% CI	-1.6931% to 8.3188%
    % Chi-squared	1.688
    % DF 	1
    % Significance level	P = 0.1938
    % https://www.medcalc.org/calc/comparison_of_proportions.php

    % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S5.png']; 
        if fast_figs
            frame = getframe(gcf); % fig is the figure handle to save
            [raster, raster_map] = frame2im(frame); % raster is the rasterized image, raster_map is the colormap
            if isempty(raster_map)
                imwrite(raster, fname);
            else
                imwrite(raster, raster_map, fname); % fig_file is the path to the image
            end            
        else
            exportgraphics(gcf,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',res);  
        end
        close(gcf);   
    end
















