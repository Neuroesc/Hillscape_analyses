%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  analysis function written for:
% Grieves, Duvelle and Taube (2024) 

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
%% >>>>>>>>>> Schematics
    xnow = 50;
    ynow = 650;
    xbuff = 50;

    % Arena
    ax1 = axes('Units','pixels','Position',[xnow ynow 200 120]);
        ah = add_panel_title('A',sprintf('%s',maze_names{7}),'yoffset',-20,'xoffset',25,'width',400,'fontsize',fs);
    
        PIT_plot_mazes(ax1,1,2)
        % title('Arena','FontSize',10,'FontWeight','normal')

        % additional text
        text(0.67,0.21,0,sprintf('3m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','center','rotation',20)
        text(0.12,0.05,0.5,sprintf('1.5m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','right','rotation',-34)
        text(0.7,0.85,0,sprintf('0.65m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','right')

    % Ridges
    ax2 = axes('Units','pixels','Position',ax1.Position+[ax1.Position(3)+xbuff 0 0 0]);
        ah = add_panel_title('B',sprintf('%s',maze_names{2}),'yoffset',-20,'xoffset',25,'width',400,'fontsize',fs);

        PIT_plot_mazes(ax2,2,2)
        % title(maze_names{2},'FontSize',10,'FontWeight','normal')

%% >>>>>>>>>> Photos
    xnow = 50;
    ynow = ynow-150;
    xbuff = 50;

    % Arena
    ax3 = axes('Units','pixels','Position',[xnow ynow ax1.Position(3) ax1.Position(4)]);
        % ah = add_panel_title('A',sprintf('Place cells recorded per rat'),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);   

        fname = [config.main_dir '\associated media\photos\Picture1.jpg'];
        im = imread(fname,'jpg');
        imshow(im);

    % Ridges
    ax4 = axes('Units','pixels','Position',ax3.Position+[ax3.Position(3)+xbuff 0 0 0]);
        % ah = add_panel_title('A',sprintf('Place cells recorded per rat'),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);   

        fname = [config.main_dir '\associated media\photos\Picture2.jpg'];
        im = imread(fname,'jpg');
        imshow(im);

        % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S1a.png']; 
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



























