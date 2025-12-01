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
    % Create figure
    fig_now = figure('Units','pixels','Position',[100 100 950 800],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

    var = 'repetition_score';
    r1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.1 & clumaa.repetition_score(:,1)>0.3,1); % arena 1 data
    r2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.1 & clumaa.repetition_score(:,1)>0.3,1); % hills data

    var = 'ratemap_planar';
    % m1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.1 & clumaa.repetition_score(:,1)>0.1); % arena 1 data
    % m1 = cellfun(@(x) imresize(x,[48 96],'bilinear'),m1,'UniformOutput',0);
    % m1 = cat(3,m1{:});   
    % m1 = mean(m1,1,'omitnan');

    m2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.1 & clumaa.repetition_score(:,1)>0.3,1); % arena 1 data
    m2 = cellfun(@(x) imresize(x,[48 96],'bilinear'),m2,'UniformOutput',0);
    m2 = cat(3,m2{:});   
    m2 = squeeze(sum(m2,1,'omitnan'));
    z2 = zscore(rot90(m2),[],2);

    u = ismember(clumaa.uci,clumaa.uci(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.1 & clumaa.repetition_score(:,1)>0.3));
    m1 = clumaa.(var)(u & clumaa.partn==1); % arena 1 data
    m1 = cellfun(@(x) imresize(x,[48 96],'bilinear'),m1,'UniformOutput',0);
    m1 = cat(3,m1{:});   
    m1 = squeeze(sum(m1,1,'omitnan'));
    z1 = zscore(rot90(m1),[],2);


    [~,midx1] = max(z1,[],2,'omitnan');
    [~,idx1] = sort(midx1,'ascend');
    z1 = z1(idx1,:);

    [~,midx2] = max(z2,[],2,'omitnan');
    [~,idx2] = sort(midx2,'ascend');
    z2 = z2(idx2,:);


ax = axes('Units','pixels','Position',[80 150 350 600]); % x-axis plot    
    imagesc(z1,'XData',1:size(z1,2)); hold on;
    axis ij on
    % daspect([1 1 1])
    ax = gca;
    colormap(ax,turbo)   
    ax.CLim = [-1 4];
    % ax.XTick = [];
    ylabel('Ranked cells')
    ax.YTick = unique([1 0:10:size(z1,1) size(z1,1)]);
    ax.YLabel.FontSize = 12;
    ax.FontSize = 12;
    xlabel('Position (bins)')

    % additional plots
    maze_sections = linspace(0,size(z1,2),7);
    hill_x = maze_sections(2:2:end);
    plot([hill_x; hill_x],ax.YLim,'w:','LineWidth',2)

% ax_f = axes('Units','pixels','Position',[ax.Position(1) ax.Position(2)-60 ax.Position(3) 55]); % x-axis plot
%     f = mean(z1,1,'omitnan');
%     s = nansem(z1,1);
% 
%     % errorbar(1:size(z1,2),f,s,'Color',plot_set{1,1},'Marker','o');
%     bb = boundedline(1:size(z1,2),f,s,'cmap',plot_set{1,1});
% 
%     xlabel('Position (bins)')
%     ax_f.XLim = ax.XLim;
%     ax_f.YLim = [-0.5 0.5];
%     ax_f.YTick = [];
%     % ax_f.XTick = ax.XTick;
%     % ax_f.YColor = 'none';        
%     box off
%     grid off
%     ytickformat('%.1f')
%     ax_f.FontSize = 12;
%     line(ax_f.XLim,[0 0],'Color',[.5 .5 .5])
% 
%     % text(ax_f,0,-1,sprintf('Elongated\nalong x-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','left')
%     % text(ax_f,1,-1,sprintf('Elongated\nalong y-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','right')



ax = axes('Units','pixels','Position',[ax.Position(1)+360 150 350 600]); % x-axis plot    
    imagesc(z2,'XData',1:size(z2,2)); hold on;
    axis ij on
    % daspect([1 1 1])
    ax = gca;
    colormap(ax,turbo)   
    ax.CLim = [-1 4];
    xlabel('Position (bins)')
    ax.YTick = [];
    % ylabel('Ranked cells')
    % ax.YTick = unique([1 0:10:size(z2,1) size(z2,1)]);
    ax.YLabel.FontSize = 12;
    ax.FontSize = 12;

    % additional plots
    maze_sections = linspace(0,size(z2,2),7);
    hill_x = maze_sections(2:2:end);
    plot([hill_x; hill_x],ax.YLim,'w:','LineWidth',2)

% ax_f = axes('Units','pixels','Position',[ax.Position(1) ax.Position(2)-60 ax.Position(3) 55]); % x-axis plot
%     f = mean(z2,1,'omitnan');
%     s = nansem(z2,1);
% 
%     % errorbar(1:size(z2,2),f,s,'Color',plot_set{1,2},'Marker','o');
%     bb = boundedline(1:size(z1,2),f,s,'cmap',plot_set{1,2});
% 
%     xlabel('Position (bins)')
%     ax_f.XLim = ax.XLim;
%     % ax_f.YTick = [];
%     % ax_f.YTick = [];
%     % ax_f.XTick = ax.XTick;
%     % ax_f.YColor = 'none';    
%     ax_f.YAxisLocation = 'right';
%     box off
%     grid off
%     ax_f.YLim = [-0.5 0.5];
%     ytickformat('%.1f')
%     ax_f.FontSize = 12;
%     line(ax_f.XLim,[0 0],'Color',[.5 .5 .5])
% 
%     % text(ax_f,0,-1,sprintf('Elongated\nalong x-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','left')
%     % text(ax_f,1,-1,sprintf('Elongated\nalong y-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','right')

%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\PIT_repetition_fieldpos.png']; 
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


return




















return
%% Parse inputs
    % Create figure
    fig_now = figure('Units','pixels','Position',[100 100 1000 800],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example cells
    xnow = 60;
    ynow = 650;
    xsiz = [180, -10]; % size, buffer
    ysiz = [80, 15]; % size, buffer
    xvec = xnow : sum(xsiz) : 750;
    yvec = ynow : -sum(ysiz)-5 : 300;
    [x,y] = ndgrid(xvec,yvec);

    % find which cells to plot  
    idx = find(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5 & clumaa.planar_spatial_info_shuffles(:,2)>30); % place cells in the hills
    [~,sidx] = sort(clumaa.repetition_score(idx,1),'descend');
    idx = idx(sidx);

    for ii = 1:numel(x)
        % create axis
        ax = axes('Units','pixels','Position',[x(ii) y(ii) xsiz(1) ysiz(1)]);
            m = clumaa.ratemap_planar{idx(ii)};

            % plot map
            imagesc(m,'alphadata',~isnan(m)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax,turbo)   
            caxis( [0 max([max(m(:),[],'omitnan') 1])] )
            ax.XTick = [];
            ax.YTick = [];
            ax.XColor = plot_set{1,2};
            ax.YColor = ax.XColor;
            ax.LineWidth = 3;

            % additional plots
            ps = linspace(0,size(m,2),7);
            tops = ps(2:2:end);
            valleys = ps(1:2:end);
            plot([tops; tops],ax.YLim,'w:')

            % text
            text(0,0,sprintf('%.1f',ax.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','w','VerticalAlignment','bottom')
            text(0,1,sprintf('Cell %d',ii),'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color','k','VerticalAlignment','bottom')
    end

%% >>>>>>>>>> repetition score
    xnow = 160;
    ynow = 100;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 200 150]);
        var = 'repetition_score';
        v1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5,1); % arena 1 data
        v2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5,1); % hills data
        v3 = clumaa.(var)(pidx & clumaa.partn==3 & clumaa.f_rate_hz>0.5,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2(:),v3(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'meanmarker',{'o'},'meanfill',1,'dotcolor',plot_set(1,:),'dotsize',40,'dotalpha',0.5,'dotsmarker',plot_set(2,:));

        % axis settings
        ylabel(sprintf('Repetition score'))    
        ax.XTick = 1:3;
        ax.XLim = [0.5 3.5]; 
        ax.XTickLabel = maze_names;
        ax.FontSize = 10;
        ax.XTickLabelRotation = 25;

        % stats
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.6);
        set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))
        N = [sum(~isnan(v1(:))) sum(~isnan(v2(:))) sum(~isnan(v3(:)))];

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])

        % legend
        % p1 = plot(-10,1,'ko','LineWidth',1.5);
        % p2 = plot(-10,1,'kv','LineWidth',1.5);
        % [~,leg] = legendflex([p1 p2],{'Balanced solution','Minimum error solution'},'anchor',{'nw','nw'},'ncol',1,'box','off','buffer',[0,20],'xscale',.5,'fontsize',9); 
        
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\PIT_repetition.png']; 
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


     % Create figure
    fig_now = figure('Units','pixels','Position',[100 100 900 900],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

    xnow = 60;
    ynow = 60;
    ax = axes('Units','pixels','Position',[xnow+280 ynow 250 500]);
        v2 = clumaa.amap_cent(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5); % hills data
        v2 = cell2mat(v2);
        v2idx = clumaa.repetition_score(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5,1); % hills data
        [~,sidx] = sort(v2idx,'descend');
    
        v2s = v2(sidx,:);
        imagesc(v2s,'XData',(1:size(v2s,2))-(size(v2s,2)./2)); hold on;
        axis ij on
        % daspect([1 1 1])
        ax = gca;
        colormap(ax,viridis)   
        ax.CLim = [-0.5 1];
        xlabel('Lag (bins)')
        ylabel('Ranked cells')
        ax.YTick = unique([1 0:25:size(v2s,1) size(v2s,1)]);
        ax.YLabel.FontSize = 12;
    
        % additional plots
        spacing_bins = 950/32;
        field_index = unique([0:spacing_bins:100 -(0:spacing_bins:100)]);
        plot([field_index; field_index],ax.YLim,'w:')

    ax = axes('Units','pixels','Position',[xnow ynow 250 500]);
        v2 = clumaa.amap_cent(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5); % hills data
        v2 = cell2mat(v2);
        v2idx = clumaa.repetition_score(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5,1); % hills data
        [~,sidx] = sort(v2idx,'descend');
    
        v2s = v2(sidx,:);
        imagesc(v2s,'XData',(1:size(v2s,2))-(size(v2s,2)./2)); hold on;
        axis ij on
        ax.YLabel.FontSize = 12;
        % daspect([1 1 1])
        ax = gca;
        colormap(ax,viridis)   
        ax.CLim = [-0.5 1];
        xlabel('Lag (bins)')
        % ylabel('Ranked cells')
        ax.YTick = unique([1 0:25:size(v2s,1) size(v2s,1)]);
    
        % additional plots
        spacing_bins = 950/32;
        field_index = unique([0:spacing_bins:100 -(0:spacing_bins:100)]);        
        plot([field_index; field_index],ax.YLim,'w:')

%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\PIT_repetition_all.png']; 
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


return













