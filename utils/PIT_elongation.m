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
    fig_now = figure('Units','pixels','Position',[100 100 1000 900],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % clumaa.planar_field_aspect_ratio = clumaa.planar_field_info(:,4) ./ clumaa.planar_field_info(:,5); % aspect ratio = field height / width
    % clumaa.surficial_field_aspect_ratio = clumaa.surficial_field_info(:,4) ./ clumaa.surficial_field_info(:,5); % aspect ratio = field height / width

    clumaa.planar_field_aspect_ratio = (clumaa.planar_field_info(:,4) - clumaa.planar_field_info(:,5)) ./ (clumaa.planar_field_info(:,4) + clumaa.planar_field_info(:,5)) ; % aspect ratio = field height / width
    clumaa.surficial_field_aspect_ratio = (clumaa.surficial_field_info(:,4) - clumaa.surficial_field_info(:,5)) ./ (clumaa.surficial_field_info(:,4) + clumaa.surficial_field_info(:,5)) ; % aspect ratio = field height / width

%% >>>>>>>>>> Example cells
    xnow = 60;
    ynow = 800;
    xsiz = [180, -10]; % size, buffer
    ysiz = [80, 15]; % size, buffer
    xvec = xnow : sum(xsiz) : 750;
    yvec = ynow : -sum(ysiz).*2-10 : 250;
    [x,y] = ndgrid(xvec,yvec);

    % find which cells to plot  
    idx = find(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5 & clumaa.repetition_score(:,1)<0.10 & clumaa.planar_spatial_info_shuffles(:,2)>30); % place cells in the hills
    [~,sidx] = sort(clumaa.planar_field_aspect_ratio(idx,1),'ascend');
    idx = idx(sidx);
    idx = idx([16 11 21 37 3 29 20 7 2 12 24 25 26 31 34]);

    for ii = 1:numel(x)
        % create axis
        ax = axes('Units','pixels','Position',[x(ii) y(ii) xsiz(1) ysiz(1)]);
            j = 0;
            m1 = clumaa.ratemap_planar{idx(ii+j)};
            m2 = clumaa.ratemap_planar{idx(ii+j)+1};

            % plot map
            imagesc(m1,'alphadata',~isnan(m1)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax,turbo)   
            ax.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
            ax.XTick = [];
            ax.YTick = [];
            ax.XColor = plot_set{1,1};
            ax.YColor = ax.XColor;
            ax.LineWidth = 3;

            % text
            text(0,0,sprintf('%.1f',ax.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','w','VerticalAlignment','bottom')
            text(0,1,sprintf('Cell %d',ii+j),'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color','k','VerticalAlignment','bottom')
            if ii==1 || ii==6 || ii==11
                text(-0.01,0.5,sprintf('Arena'),'Units','normalized','HorizontalAlignment','center','FontSize',12,'Color','k','VerticalAlignment','bottom','rotation',90,'Color',ax.XColor)
            end

        % create axis
        ax = axes('Units','pixels','Position',[x(ii) y(ii)-sum(ysiz) xsiz(1) ysiz(1)]);
            % plot map
            imagesc(m2,'alphadata',~isnan(m2)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax,turbo)   
            ax.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
            ax.XTick = [];
            ax.YTick = [];
            ax.XColor = plot_set{1,2};
            ax.YColor = ax.XColor;
            ax.LineWidth = 3;

            % additional plots
            ps = linspace(0,size(m2,2),7);
            tops = ps(2:2:end);
            valleys = ps(1:2:end);
            plot([tops; tops],ax.YLim,'w:')

            if ii==1 || ii==6 || ii==11
                text(-0.01,0.5,sprintf('Hills'),'Units','normalized','HorizontalAlignment','center','FontSize',12,'Color','k','VerticalAlignment','bottom','rotation',90,'Color',ax.XColor)
            end
    end
% return

%% >>>>>>>>>> aspect ratio
    xnow = 160;
    ynow = 100;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 200 150]);
        var = 'planar_field_aspect_ratio';
        v1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5 & clumaa.repetition_score(:,1)<0.5,1); % arena 1 data
        v2a = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5 & clumaa.repetition_score(:,1)<0.5,1); % hills data
        v2b = clumaa.surficial_field_aspect_ratio(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5 & clumaa.repetition_score(:,1)<0.5,1); % hills data
        v3 = clumaa.(var)(pidx & clumaa.partn==3 & clumaa.f_rate_hz>0.5 & clumaa.repetition_score(:,1)<0.5,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'meanmarker',{'o'},'meanfill',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',40,'dotalpha',0.5,'dotsmarker',plot_set(2,[1 2 2 3]));

        % axis settings
        ylabel(sprintf('Aspect ratio (Y/X)'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        ax.YLim = [-1 1];
        ax.XTickLabel = [];
        ax.FontSize = 10;
        ax.XTickLabelRotation = 25;
        % ax.YScale = 'log';
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s\n(flattened)',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);

%% >>>>>>>>>> orientation
    xnow = 460;
    ynow = 100;

    % create axis
    ax = polaraxes('Units','pixels','Position',[xnow ynow 200 150]);
        var = 'planar_field_info';
        v1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5 & clumaa.repetition_score(:,1)<0.5,6); % arena 1 data
        v2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5 & clumaa.repetition_score(:,1)<0.5,6); % hills data
        v3 = clumaa.(var)(pidx & clumaa.partn==3 & clumaa.f_rate_hz>0.5 & clumaa.repetition_score(:,1)<0.5,6); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2(:),v3(:)); % linearise data

        % get circular density of data
        ri = linspace(-pi,pi,360)';
        xi = movmean(ri,2,'Endpoints','discard');
        f1 = circ_ksdensity(deg2rad(v1).*2,xi,[],0.25);
        f2 = circ_ksdensity(deg2rad(v2).*2,xi,[],0.25);
        ri = ri(:)';
        f1 = f1(:)';
        f2 = f2(:)';

        % main plot
        p1 = polarhistogram('BinEdges',ri,'BinCounts',f1,'DisplayStyle','bar','EdgeColor','none','FaceColor',plot_set{1,1}); hold on;
        p2 = polarhistogram('BinEdges',ri,'BinCounts',f2,'DisplayStyle','bar','EdgeColor','none','FaceColor',plot_set{1,2});

        % axis settings
        ax.RTick = [];
        ax.RAxisLocation = 180;
        ax.ThetaZeroLocation = 'bottom';
        ax.LineWidth = 1.5;

        % additional plots
        m1 = circ_mean(deg2rad(v1).*2);
        m2 = circ_mean(deg2rad(v2).*2);
        p3 = polarplot([m1 m1],[0 1],'Color',plot_set{1,1},'LineWidth',1.5); hold on;
        p4 = polarplot([m2 m2],[0 1],'Color',plot_set{1,2},'LineWidth',1.5); hold on;
        text(0.5,-0.25,sprintf('Aligned to Y-axis'),'Units','normalized','HorizontalAlignment','center','FontSize',10,'Color','k','VerticalAlignment','bottom','rotation',0)
        text(0.5,1.15,sprintf('Aligned to X-axis'),'Units','normalized','HorizontalAlignment','center','FontSize',10,'Color','k','VerticalAlignment','bottom','rotation',0)

        % stats
        rv1 = circ_r(deg2rad(v1).*2);
        rv2 = circ_r(deg2rad(v2).*2);
        text(0.66,0.27,sprintf('r = %.2f',rv1),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color',plot_set{1,1},'VerticalAlignment','middle','rotation',0)
        text(0.66,0.73,sprintf('r = %.2f',rv2),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color',plot_set{1,2},'VerticalAlignment','middle','rotation',0)

        % legend
        axt = axes('Units','pixels','Position',[xnow ynow 200 150],'Color','none');
            axt.XLim = [0 1];
            axt.YLim = [0 1];
            axis off
            p1 = patch(repmat(-10,4,1),repmat(-10,4,1),p1.FaceColor,'FaceAlpha',p1.FaceAlpha,'EdgeColor',p1.EdgeColor); hold on;
            p2 = patch(repmat(-10,4,1),repmat(-10,4,1),p2.FaceColor,'FaceAlpha',p2.FaceAlpha,'EdgeColor',p2.EdgeColor);            
            p3 = plot(-10,-10,'Color',p3.Color,'LineWidth',p3.LineWidth);
            p4 = plot(-10,-10,'Color',p4.Color,'LineWidth',p4.LineWidth);
            [~,leg] = legendflex([p1,p2,p3,p4],{'Arena 1','Hills','Arena 1 mean','Hills mean'},'anchor',{'e','e'},'ncol',1,'box','off','buffer',[105,0],'xscale',.5,'fontsize',9); 
            leg(5).FaceAlpha = p1.FaceAlpha;
            leg(6).FaceAlpha = p2.FaceAlpha;

%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\PIT_elongation.png']; 
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


 


















