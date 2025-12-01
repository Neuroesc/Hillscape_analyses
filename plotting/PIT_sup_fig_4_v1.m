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
    fs = [15 10];

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example planimetric cells
    xnow = 20;
    ynow = 700;
    xsiz = [130, 20]; % size, buffer
    ysiz = [80, 0]; % size, buffer
    xvec = xnow : sum(xsiz) : 10000;
    xvec = xvec(1:4);
    yvec = ynow : -sum(ysiz).*2.2 : -10000;
    yvec = yvec(1:4);
    [x,y] = ndgrid(xvec,yvec);

    % find which cells to plot  
    ucis = unique(clumaa.uci(pidx),'stable');

    % rank cells by correlation between arena 1 and hills
    v = clumaa.between_session_stability( ismember(clumaa.uci,ucis) & clumaa.partn==1 , 1);
    [~,ridx] = sort(v,'descend','MissingPlacement','last');
    ucis = ucis(ridx);
    ncells = numel(xvec).*numel(yvec); % number of cells per page
    ucis = ucis([3 4 5 7 8 9 11 13 14 15 16 17 18 19 20 21]);

    % plot figure
    count = 0;
    fig_count = 1;
    disp(sprintf('Fig %d...',fig_count))    
    for ii = 1:numel(ucis)
        uci = ucis{ii};
        disp(sprintf('\t...%s',uci))
        m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
        m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
        m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};
        vn = clumaa.between_session_stability( ismember(clumaa.uci,uci) & clumaa.partn==1 , 1);

        % hills ratemap
        jj = ii-((fig_count-1).*24);
        ax1 = axes('Units','pixels','Position',[x(jj) y(jj) xsiz(1) ysiz(1)]);
            % plot map
            imagesc(m2,'alphadata',~isnan(m2)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax1,turbo)   
            ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') max(m3(:),[],'omitnan') 1])];
            ax1.XTick = [];
            ax1.YTick = [];

            % additional plots
            ps = linspace(0,size(m1,2),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax1.YLim,'w:')

            % text
            text(0,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',9,'Color','k','VerticalAlignment','bottom','rotation',0)
            % text(0,0,sprintf('%s',uci),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',90,'Interpreter','none')
            text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
            text(1,1,sprintf('r = %.3f',vn),'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)

        % arena 1 ratemap
        ax2 = axes('Units','pixels','Position',[x(jj) y(jj)+ax1.Position(4)+ysiz(2) ax1.Position(3) ax1.Position(4)]);
            % plot map
            imagesc(m1,'alphadata',~isnan(m1)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax2,ax1.Colormap)   
            ax2.CLim = ax1.CLim;
            ax2.XTick = [];
            ax2.YTick = [];

            % text
            text(0,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',9,'Color','k','VerticalAlignment','bottom','rotation',0)
            text(1,1,sprintf('Cell %d',ii),'Units','normalized','HorizontalAlignment','right','FontSize',9,'Color','k','VerticalAlignment','bottom')

        % % arena 2 ratemap
        % ax3 = axes('Units','pixels','Position',[x(jj)+ax1.Position(3)-(ax1.Position(3).*0.48) y(jj)+ax1.Position(4)+ysiz(2) ax1.Position(3).*0.48 ax1.Position(4).*0.48]);
        %     % plot map
        %     imagesc(m3,'alphadata',~isnan(m3)); hold on;
        % 
        %     % axis settings
        %     axis xy on
        %     daspect([1 1 1])
        %     colormap(ax3,ax1.Colormap)   
        %     ax3.CLim = ax1.CLim;
        %     ax3.XTick = [];
        %     ax3.YTick = [];
        % 
        %     % text
        %     text(0,1,sprintf('%s',maze_names{3}),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)
    end

    % colorbar
    axc = axes('Units','pixels','Position',[ax1.Position(1)-430 ax1.Position(2)-20 80 8]);
        mat = (linspace(0,100,100));
        imagesc(ones(size(mat)),mat,mat);
        colormap(axc,turbo);
        axis xy off

        axc.YTick = [];
        axc.XTick = [];
        text(0.5,1.7,sprintf('Firing rate (Hz)'),'FontSize',7,'HorizontalAl','center','Units','normalized')
        text(1.02,0.5,sprintf('Max'),'FontSize',7,'HorizontalAl','left','Units','normalized')
        text(-0.02,0.5,sprintf('0'),'FontSize',7,'HorizontalAl','right','Units','normalized')

    % keyboard

%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S4.png']; 
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




% %% >>>>>>>>>> Example anti-planimetric cells
%     xnow = 20;
%     ynow = 330;
%     xsiz = [130, 20]; % size, buffer
%     ysiz = [80, 0]; % size, buffer
%     xvec = xnow : sum(xsiz) : 10000;
%     xvec = xvec(1:4);
%     yvec = ynow : -sum(ysiz).*1.75 : -10000;
%     yvec = yvec(1:3);
%     [x,y] = ndgrid(xvec,yvec);
% 
%     % find which cells to plot  
%     ucis = unique(clumaa.uci(pidx),'stable');
% 
%     % rank cells by correlation between arena 1 and hills
%     v = clumaa.between_session_stability( ismember(clumaa.uci,ucis) & clumaa.partn==1 , 1);
%     [~,ridx] = sort(v,'ascend','MissingPlacement','last');
%     ucis = ucis(ridx);
%     ncells = numel(xvec).*numel(yvec); % number of cells per page
%     % ucis = ucis([3 4 5 7 8 9 11 13 14 15 16 17]);
% 
%     % plot figure
%     count = 0;
%     fig_count = 1;
%     disp(sprintf('Fig %d...',fig_count))    
%     for ii = 1:ncells
%         uci = ucis{ii};
%         disp(sprintf('\t...%s',uci))
%         m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
%         m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
%         m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};
%         vn = clumaa.between_session_stability( ismember(clumaa.uci,uci) & clumaa.partn==1 , 1);
% 
%         % hills ratemap
%         jj = ii-((fig_count-1).*24);
%         ax1 = axes('Units','pixels','Position',[x(jj) y(jj) xsiz(1) ysiz(1)]);
%             % plot map
%             imagesc(m2,'alphadata',~isnan(m2)); hold on;
% 
%             % axis settings
%             axis xy on
%             daspect([1 1 1])
%             colormap(ax1,turbo)   
%             ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') max(m3(:),[],'omitnan') 1])];
%             ax1.XTick = [];
%             ax1.YTick = [];
% 
%             % additional plots
%             ps = linspace(0,size(m1,2),7);
%             tops = ps(2:2:end);
%             plot([tops; tops],ax1.YLim,'w:')
% 
%             % text
%             text(0,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)
%             text(0,0,sprintf('%s',uci),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',90,'Interpreter','none')
%             text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
%             text(1,1,sprintf('r = %.3f',vn),'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)
% 
%         % arena 1 ratemap
%         ax2 = axes('Units','pixels','Position',[x(jj) y(jj)+ax1.Position(4)+ysiz(2) ax1.Position(3).*0.48 ax1.Position(4).*0.48]);
%             % plot map
%             imagesc(m1,'alphadata',~isnan(m1)); hold on;
% 
%             % axis settings
%             axis xy on
%             daspect([1 1 1])
%             colormap(ax2,ax1.Colormap)   
%             ax2.CLim = ax1.CLim;
%             ax2.XTick = [];
%             ax2.YTick = [];
% 
%             % text
%             text(0,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)
% 
%         % arena 2 ratemap
%         ax3 = axes('Units','pixels','Position',[x(jj)+ax1.Position(3)-(ax1.Position(3).*0.48) y(jj)+ax1.Position(4)+ysiz(2) ax1.Position(3).*0.48 ax1.Position(4).*0.48]);
%             % plot map
%             imagesc(m3,'alphadata',~isnan(m3)); hold on;
% 
%             % axis settings
%             axis xy on
%             daspect([1 1 1])
%             colormap(ax3,ax1.Colormap)   
%             ax3.CLim = ax1.CLim;
%             ax3.XTick = [];
%             ax3.YTick = [];
% 
%             % text
%             text(0,1,sprintf('%s',maze_names{3}),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)
%     end










