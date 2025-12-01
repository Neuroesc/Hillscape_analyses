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

%% >>>>>>>>>> Diagram explaining field to wall angle
    xnow = 50;
    ynow = 600;
    k = 1.5;

    % example ellipses
    ax2 = axes('Units','pixels','Position',[xnow-160 ynow 250 80],'Color','none');
        ah = add_panel_title('A',sprintf(''),'yoffset',0,'xoffset',150,'width',400,'fontsize',fs);  
    
        ax2.XLim = [-12 10];
        ax2.YLim = [-12 10];
        daspect([1 1 1])
        axis manual off
        ax2.Clipping = 'off';
    
        ovals = [0:15:90];
        xvals = (1:length(ovals)).*20;
        yvals = ones(size(xvals)).*4;
        ye = 0;   
        cols = winter(length(xvals));
        for xx = 1:length(xvals)
            emin = 5;
            emax = emin/(1-0.7^2);
    
            e = ellipse(emax,emin,deg2rad(ovals(xx)),xvals(xx),yvals(xx));
            x = e.XData;
            y = e.YData;
            y = y - min(y(:));
            patch(x,y,'k','EdgeColor','k','LineWidth',1,'FaceAlpha',0.3,'Clipping','off','EdgeAlpha',0); hold on;
            delete(e);
            text(mean(x(:)),max(y(:)),sprintf('%.f%c',ovals(xx),176),'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','bottom')

            hlen = [10 10];
            xCentre = mean(x(:));
            yCentre = mean(y(:));
            plot(xCentre,yCentre,'k+','MarkerSize',15)
            cosOrient = cosd(-ovals(xx));
            sinOrient = sind(-ovals(xx));
            xcoords = xCentre + hlen .* [cosOrient -cosOrient];
            ycoords = yCentre + hlen .* [-sinOrient sinOrient];
            patch(xcoords,ycoords,'k','LineWidth',1,'EdgeColor','k','FaceAlpha',0,'Clipping','off','EdgeAlpha',1,'LineStyle',':'); hold on;

            hlen = [10 14];            
            cosOrient = cosd(0);
            sinOrient = sind(0);
            xcoords = xCentre + hlen .* [cosOrient -cosOrient];
            ycoords = 0 + hlen .* [-sinOrient sinOrient];
            patch(xcoords,ycoords,'k','LineWidth',1,'EdgeColor','k','FaceAlpha',0,'Clipping','off','EdgeAlpha',1,'LineStyle','-'); hold on;

        end
        text(10,0,'Closest wall','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8)

%% >>>>>>>>>> field-to-wall angle vs distance from closest boundary
    xnow = 50;
    ynow = ynow-150;

    % create axis
    ax1 = axes('Units','pixels','Position',[xnow ynow 125 125]);
        ah = add_panel_title('B',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % field-to-wall angle function
        fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
        
        % get data
        pp = 1; % arena 1
        f1 = datn{pp};
        field_orientation = f1(:,3);
        idx = field_orientation==90;
        field_orientation(idx) = field_orientation(idx)-0.01; % fix 90 degree field orientation values
        wall_orientation = f1(:,9);
        idx = wall_orientation==90;
        wall_orientation(idx) = wall_orientation(idx)-0.01; % fix 90 degree field orientation values      
        wall_distance = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        ftw = fun([field_orientation wall_orientation]); % calculate field-to-wall angle

        % plot data
        d = computeScatterDensity(wall_distance,ftw,'k',k);
        scatter(wall_distance,ftw,20,d,'filled','MarkerEdgeColor','none'); hold on;
        r = refline;
        set(r,'Color','k')

        % axis settings
        grid on;
        colormap(ax1,'turbo')
        ax1.CLim(1) = 0;
        ax1.YLim = [0 90];
        ax1.XLim = [0 0.75];
        ax1.YTick = 0:30:90;
        ax1.XTick = 0:0.25:0.75;        
        ylabel(sprintf('Field-to-wall angle (%c)',176))
        xlabel(sprintf('Field-to-wall distance (m)'))
        ytickformat('%.f')

        % text
        text(0.1,1.02,sprintf('%s',maze_names{1}),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz)
        % r = corr(wall_distance,ftw,'rows','pairwise','type','Spearman');
        % text(0.1,1,sprintf('{\\itr} = %.2f',r),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',9)

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

    ax1 = axes('Units','pixels','Position',[ax1.Position(1)+170 ynow 80 125]);
        cut_dist = 0.2;
    
        % get data
        f1 = datn{1}; % arena 1
        wall_distance1 = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        idx1 = wall_distance1 < cut_dist; % fields close to wall
    
        f1 = datn{2}; % hills
        wall_distance2 = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        idx2 = wall_distance2 < cut_dist; % fields close to wall
    
        % vectorise
        [ds,gs] = vectorDATAGROUP([],ftw(idx1),ftw(~idx1)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[1 1]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        ylabel(sprintf('Field-to-wall angle (%c)',176))
        ax1.XTick = 1:2;
        ax1.YTick = 0:30:90;
        ax1.XLim = [0.5 2.5]; 
        ax1.YLim = [0 90];
        ax1.XTickLabel = {'Wall','Centre'};
        ax1.FontSize = 8;
        ytickformat('%.f');
    
        % stats
        axt = axes('Units','pixels','Position',ax1.Position,'Color','none');
            axis off
            axt.XLim = ax1.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.05,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

    % Hills
    % create axis
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+140 ynow 125 125]);
        ah = add_panel_title('C',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % field-to-wall angle function
        fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
        
        % get data
        pp = 2; % arena 1
        f1 = datn{pp};
        field_orientation = f1(:,3);
        idx = field_orientation==90;
        field_orientation(idx) = field_orientation(idx)-0.01; % fix 90 degree field orientation values
        wall_orientation = f1(:,9);
        idx = wall_orientation==90;
        wall_orientation(idx) = wall_orientation(idx)-0.01; % fix 90 degree field orientation values      
        wall_distance = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        ftw = fun([field_orientation wall_orientation]); % calculate field-to-wall angle

        % plot data
        d = computeScatterDensity(wall_distance,ftw,'k',k);
        scatter(wall_distance,ftw,20,d,'filled','MarkerEdgeColor','none'); hold on;
        r = refline;
        set(r,'Color','k')

        % axis settings
        grid on
        colormap(ax2,'turbo')
        ax2.CLim(1) = 0;
        ax2.YLim = [0 90];
        ax2.XLim = [0 0.75];
        ax2.YTick = 0:30:90;    
        ax2.XTick = 0:0.25:0.75;
        ylabel(sprintf('Field-to-wall angle (%c)',176))
        xlabel(sprintf('Field-to-wall distance (m)'))
        ytickformat('%.f')

        % text
        text(0.1,1.02,sprintf('%s',maze_names{2}),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz)

    ax2 = axes('Units','pixels','Position',[ax2.Position(1)+170 ax2.Position(2) 80 125]);            
        % vectorise
        [ds,gs] = vectorDATAGROUP([],ftw(idx2),ftw(~idx2)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[2 2]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        ylabel(sprintf('Field-to-wall angle (%c)',176))
        ax2.XTick = 1:2;
        ax2.YTick = 0:30:90;  
        ax2.YLim = [0 90];
        ax2.XLim = [0.5 2.5]; 
        ax2.XTickLabel = {'Wall','Centre'};
        ax2.FontSize = 8;
        ytickformat('%.f');

        % stats
        axt = axes('Units','pixels','Position',ax2.Position,'Color','none');
            axis off
            axt.XLim = ax2.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.05,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

    % Hills (including hilltops)
    % create axis
    ax2 = axes('Units','pixels','Position',[xnow ynow-210 125 125]);
        ah = add_panel_title('D',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % field-to-wall angle function
        fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
        
        % get data
        pp = 2; % arena 1
        f1 = datn{pp};
        field_orientation = f1(:,3);
        idx = field_orientation==90;
        field_orientation(idx) = field_orientation(idx)-0.01; % fix 90 degree field orientation values
        wall_orientation = f1(:,11);
        idx = wall_orientation==90;
        wall_orientation(idx) = wall_orientation(idx)-0.01; % fix 90 degree field orientation values      
        wall_distance = f1(:,7) .* 32 ./ 1000; % distance to closest wall
        ftw = fun([field_orientation wall_orientation]); % calculate field-to-wall angle

        % plot data
        d = computeScatterDensity(wall_distance,ftw,'k',k);
        scatter(wall_distance,ftw,20,d,'filled','MarkerEdgeColor','none'); hold on;
        r = refline;
        set(r,'Color','k')

        % axis settings
        grid on
        colormap(ax2,'turbo')
        ax2.CLim(1) = 0;
        ax2.YLim = [0 90];
        ax2.XLim = [0 0.75];
        ax2.YTick = 0:30:90;    
        ax2.XTick = 0:0.25:0.75;
        ylabel(sprintf('Field-to-wall angle (%c)',176))
        xlabel(sprintf('Field-to-wall distance (m)'))
        ytickformat('%.f')

        % text
        text(0.1,1.02,sprintf('%s (+ ridgetops)',maze_names{2}),'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','FontSize',fsiz)

    ax2 = axes('Units','pixels','Position',[ax2.Position(1)+170 ax2.Position(2) 80 125]);            
        % vectorise
        [ds,gs] = vectorDATAGROUP([],ftw(idx2),ftw(~idx2)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[2 2]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        ylabel(sprintf('Field-to-wall angle (%c)',176))
        ax2.XTick = 1:2;
        ax2.YTick = 0:30:90;  
        ax2.YLim = [0 90];
        ax2.XLim = [0.5 2.5]; 
        ax2.XTickLabel = {'Wall','Centre'};
        ax2.FontSize = 8;
        ytickformat('%.f');

        % stats
        axt = axes('Units','pixels','Position',ax2.Position,'Color','none');
            axis off
            axt.XLim = ax2.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.05,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

%% >>>>>>>>>> Field density vs distance from closest boundary
    xnow = xnow+310;
    ynow = ynow-210;

    % create axis
    ax2 = axes('Units','pixels','Position',[xnow ynow 240 125]);
        ah = add_panel_title('E',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);

        % get data
        f1 = datn{1}; % arena 1
        wall_distance1 = f1(:,7) .* 32 ./ 1000; % distance from field to closest wall
        f2 = datn{2}; % hills
        wall_distance2 = f2(:,7) .* 32 ./ 1000; % distance from field to closest wall

        xi = 0:0.01:0.75;
        bw = 0.02;
        f1 = ksdensity(wall_distance1(:),xi(:),"Bandwidth",bw);
        f2 = ksdensity(wall_distance2(:),xi(:),"Bandwidth",bw);

        % main plot    
        alph = 0.7;          
        a1 = area(xi,f1,'FaceColor',plot_set{1,1},'EdgeColor','k','FaceAlpha',alph); hold on;
        a2 = area(xi,f2,'FaceColor',plot_set{1,2},'EdgeColor','k','FaceAlpha',alph);       

        % axis settings
        ax2.XLim = [0 0.75]; 
        ax2.FontSize = 8;
        ytickformat('%.1f')
        box off
        xlabel(sprintf('Field-to-wall distance (m)'))    
        ylabel(sprintf('PDF'))  

        [~,leg] = legendflex([a1 a2],...
            {maze_names{1},maze_names{2}},...
            'anchor',{'ne','ne'},'ncol',1,'box','off','buffer',[-10,10],'xscale',.5,'fontsize',9); 
        
        leg(3).Children.FaceAlpha = a1.FaceAlpha;
        leg(4).Children.FaceAlpha = a2.FaceAlpha;
        leg(3).Children.FaceColor = a1.FaceColor;
        leg(4).Children.FaceColor = a2.FaceColor;

        % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S6.png']; 
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
















