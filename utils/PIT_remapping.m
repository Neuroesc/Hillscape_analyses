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
%% >>>>>>>>>> Example cells
    xnow = 60;
    ynow = 800;
    xsiz = [180, -10]; % size, buffer
    ysiz = [80, 15]; % size, buffer
    xvec = xnow : sum(xsiz) : 750;
    yvec = ynow : -sum(ysiz).*2-10 : 250;
    [x,y] = ndgrid(xvec,yvec);

    % find which cells to plot  
    idx = find(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.1 & clumaa.planar_spatial_info_shuffles(:,2)>30); % place cells in the hills
    [~,sidx] = sort(clumaa.planar_spatial_info_shuffles(idx,1),'descend');
    idx = idx(sidx);
    idx = idx([1 2 3 5 7 8 9 10 13 12 14 15 16 17 18]);

    for ii = 1:numel(x)
        % hills ratemap
        ax1 = axes('Units','pixels','Position',[x(ii) y(ii)-sum(ysiz) xsiz(1) ysiz(1)]);
            % plot map
            imagesc(m2,'alphadata',~isnan(m2)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax1,turbo)   
            ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') max(m3(:),[],'omitnan') 1])];
            ax1.XTick = [];
            ax1.YTick = [];
            ax1.XColor = plot_set{1,2};
            ax1.YColor = ax1.XColor;
            ax1.LineWidth = 3;

            % additional plots
            ps = linspace(0,size(m1,2),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax1.YLim,'w:')

        % arena 1 ratemap
        ax2 = axes('Units','pixels','Position',[x(ii) y(ii)+ax1.Position(4)+ysiz(2) ax1.Position(3).*0.48 ax1.Position(4).*0.48]);
            % plot map
            imagesc(m1,'alphadata',~isnan(m1)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax2,ax1.Colormap)   
            ax2.CLim = ax1.CLim;
            ax2.XTick = [];
            ax2.YTick = [];

        % arena 2 ratemap
        ax3 = axes('Units','pixels','Position',[x(ii) y(ii)+ax1.Position(3)-(ax1.Position(3).*0.48) ax1.Position(3).*0.48 ax1.Position(4).*0.48]);
            % plot map
            imagesc(m3,'alphadata',~isnan(m3)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax3,ax1.Colormap)   
            ax3.CLim = ax1.CLim;
            ax3.XTick = [];
            ax3.YTick = [];

    end
return

%% >>>>>>>>>> correlations
    frate_cutoff = 1;
    if ~any(ismember(clumaa.Properties.VariableNames,'map_correlations'))
        ucis = unique(clumaa.uci);
        clumaa.map_correlations = NaN(size(clumaa,1),3); % preallocate
        
        for ii = 1:length(ucis)
            % normal correlations
            m1 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1}; % arena 1
            m2 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==2}; % hills
            m3 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==3}; % arena 2
    
            if max(m1(:),[],'omitnan')>=frate_cutoff || max(m2(:),[],'omitnan')>=frate_cutoff
                r1 = corr(m1(:),m2(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills
            end
            if max(m1(:),[],'omitnan')>=frate_cutoff || max(m3(:),[],'omitnan')>=frate_cutoff            
                r2 = corr(m1(:),m3(:),'rows','pairwise','type','Pearson'); % arena 1 vs arena 2
            end
    
            % correlations limited to overlapping bins
            ps = linspace(0,size(m1,2),7);
            tops = ps(2:2:end);
            valleys = ps(1:2:end);
            [cc,rr] = meshgrid(1:size(m1,1),1:size(m1,2));
            i = knnsearch(ps(:),(1:size(m1,2))');
    
            p_idx = 2:2:length(ps);
            p_log = ismember(repmat(i',size(m1,1),1),p_idx); % matrix same size as m1, 1s for points closes to peaks, 0s otherwise
            t_idx = 1:2:length(ps);
            t_log = ismember(repmat(i',size(m1,1),1),t_idx); % matrix same size as m1, 1s for points closes to troughs, 0s otherwise
    
            m1b = m1(t_log);
            m2b = m2(t_log);
            if max(m1b(:),[],'omitnan')>=frate_cutoff || max(m2b(:),[],'omitnan')>=frate_cutoff            
                r2b = corr(m1b(:),m2b(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills (troughs only)
            end

            clumaa.map_correlations(ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1,:) = [r1 r2 r2b];
        end
    end

%% >>>>>>>>>> correlation results
    xnow = 200;
    ynow = 100;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 200 150]);
        var = 'map_correlations';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,2); % arena 1 vs arena 2
        v2 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 vs hills
        v3 = clumaa.(var)(pidx & clumaa.partn==1,3); % arena 1 vs hills (troughs only)
        [ds,gs] = vectorDATAGROUP([],v1(:),v2(:),v3(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'meanmarker',{'o'},'meanfill',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',40,'dotalpha',0.5,'dotsmarker',plot_set(2,[1 2 2 3]));

        % axis settings
        ylabel(sprintf('Correlation (r)'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        % ax.XTickLabel = {[maze_names{1} ' vs ' maze_names{3}],[maze_names{1} ' vs ' maze_names{2}],[maze_names{1} ' vs ' maze_names{2} ' (overlapping)']};
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s vs %s',maze_names{1},maze_names{3}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s vs %s',maze_names{1},maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s vs %s\n(overlapping)',maze_names{1},maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')

        ax.FontSize = 10;
        ax.XTickLabelRotation = 25;
        % ax.YScale = 'log';
        % ax.YLim(1) = 0;
        % ax.YTick = [10 20 40 80 100];
        set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);
            N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];
        
% return
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\PIT_remapping.png']; 
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


 


















