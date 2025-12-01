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
    fig_now = figure('Units','pixels','Position',[100 100 1000 800],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    clumaa.planar_field_aspect_ratio = clumaa.planar_field_info(:,4) ./ clumaa.planar_field_info(:,5); % aspect ratio = field height / width

%% >>>>>>>>>> Example cells
    xnow = 60;
    ynow = 600;
    xsiz = [180, -10]; % size, buffer
    ysiz = [80, 15]; % size, buffer
    xvec = xnow : sum(xsiz) : 750;
    yvec = ynow : -sum(ysiz).*2-10 : 350;
    [x,y] = ndgrid(xvec,yvec);

    % find which cells to plot  
    idx = find(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.1 & clumaa.repetition_score(:,1)<0.20 & clumaa.planar_spatial_info_shuffles(:,2)>40); % place cells in the hills
    [~,sidx] = sort(clumaa.planar_spatial_info_shuffles(idx,1),'descend');
    idx = idx(sidx);
    idx = idx([1 2 3 5 6 7 8 9 10 12]);

    for ii = 1:numel(x)
        % create axis
        ax = axes('Units','pixels','Position',[x(ii) y(ii) xsiz(1) ysiz(1)]);
            m1 = clumaa.ratemap_planar{idx(ii)-1};
            m2 = clumaa.ratemap_planar{idx(ii)};

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
            text(0,1,sprintf('Cell %d',ii),'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color','k','VerticalAlignment','bottom')

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
            caxis( [0 max([max(m2(:),[],'omitnan') 1])] )
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

%% >>>>>>>>>> field radius
    xnow = 160;
    ynow = 100;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 200 150]);
        var = 'planar_field_info';
        v1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5,1); % arena 1 data
        v2a = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5,1); % hills data
        v2b = clumaa.surficial_field_info(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5,1); % hills data        
        v3 = clumaa.(var)(pidx & clumaa.partn==3 & clumaa.f_rate_hz>0.5,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'meanmarker',{'o'},'meanfill',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',40,'dotalpha',0.5,'dotsmarker',plot_set(2,[1 2 2 3]));

        % axis settings
        ylabel(sprintf('Field radius (cm)'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        % ax.XTickLabel = {maze_names{1} maze_names{2} sprintf([maze_names{2} ' (flattened)']) maze_names{3}};
        ax.FontSize = 10;
        % ax.XTickLabelRotation = 25;
        ax.YScale = 'log';
        % ax.YLim(1) = 0;
        ax.YTick = [10 20 40 80 100];
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s\n(flattened)',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')

        % additional plots
        % line(ax.XLim,[2 2],'Color',[.5 .5 .5]) 

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);
            N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];
        
% %% >>>>>>>>>> spatial info
%     xnow = 420;
%     ynow = 100;
% 
%     % create axis
%     ax = axes('Units','pixels','Position',[xnow ynow 200 150]);
%         var = 'planar_spatial_info_shuffles';
%         v1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5,2); % arena 1 data
%         v2a = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5,2); % hills data
%         v2b = clumaa.surficial_spatial_info_shuffles(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5,2); % hills data        
%         v3 = clumaa.(var)(pidx & clumaa.partn==3 & clumaa.f_rate_hz>0.5,2); % arena 2 data
%         [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data
% 
%         % main plot
%         meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'meanmarker',{'o'},'meanfill',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',40,'dotalpha',0.5,'dotsmarker',plot_set(2,[1 2 2 3]));
% 
%         % axis settings
%         ylabel(sprintf('Spatial info (z)'))    
%         ax.XTick = 1:4;
%         ax.XLim = [0.5 3.5]; 
%         % ax.XTickLabel = {maze_names{1} maze_names{2} sprintf([maze_names{2} ' (flattened)']) maze_names{3}};
%         ax.FontSize = 10;
%         % ax.XTickLabelRotation = 25;
%         % ax.YScale = 'log';
%         % ax.YLim(1) = 0;
%         ax.XTickLabel = {};
%         text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
%         text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
%         text(5/6,-0.01,sprintf('%s\n(flattened)',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
% 
%         % additional plots
%         line(ax.XLim,[2 2],'Color',[.5 .5 .5]) 
% 
%         % stats
%         axt = axes('Units','pixels','Position',ax.Position,'Color','none');
%             axis off
%             axt.XLim = ax.XLim;
%             axt.YLim = [0 1];        
%             [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);


%% >>>>>>>>>> fields
    xnow = 420;
    ynow = 100;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 200 150]);
        var = 'planar_fields';
        v1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.5,1); % arena 1 data
        v2a = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5,1); % hills data
        v2b = clumaa.surficial_fields(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5,1); % hills data        
        v3 = clumaa.(var)(pidx & clumaa.partn==3 & clumaa.f_rate_hz>0.5,1); % arena 2 data

        % convert to fields per m2
        v1 = v1 ./ (3*1.2);
        v2a = v2a ./ (3*1.2);
        v2b = v2b ./ (3.6*1.2);
        v3 = v3 ./ (3*1.2);

        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'meanmarker',{'o'},'meanfill',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',40,'dotalpha',0.5,'dotsmarker',plot_set(2,[1 2 2 3]));

        % axis settings
        ylabel(sprintf('Place fields / m^{2}'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        % ax.XTickLabel = {maze_names{1} maze_names{2} sprintf([maze_names{2} ' (flattened)']) maze_names{3}};
        ax.FontSize = 10;
        % ax.XTickLabelRotation = 25;
        % ax.YScale = 'log';
        % ax.YLim(1) = 0;
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s\n(flattened)',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        % set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))
        ytickformat('%.1f');

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);



%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\PIT_localisation.png']; 
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


 


















