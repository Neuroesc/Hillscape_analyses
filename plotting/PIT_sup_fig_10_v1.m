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
    xnow = 50;
    ynow = 680;
    xsiz = 200;
    ysiz = 120;
    xbuff = xsiz+30;
    cmap = cmocean('balance',128);
    % cmap = parula(128);

    % plot height
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('a',sprintf('Example trajectory'),'yoffset',0,'xoffset',0,'width',400);  

        % get data
        [r,c] = size(posdata.anisotropy_map{1,2});
        [rr,cc] = ndgrid(1:r,1:c);

        ncos = c;
        xi1 = linspace(-pi,5*pi,ncos);   
        height = (0.5 .* (cos(xi1) + 1)) .* 450;
        height_mat = repmat(height,r,1);

        % plot data
        imagesc('XData',[0 3],'YData',[0 1.5],'CData',height_mat); hold on;

        % axis settings
        daspect([1 1 1])
        axis on
        colormap(ax1,cmap)
        ax1.CLim = [0 450];
        xlabel('X (m)')
        ylabel('Y (m)')
        ytickformat('%.1f')
        xtickformat('%.1f')

        % additional plots
        ps = linspace(ax1.XLim(1),ax1.XLim(2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax1.YLim,'y--')  

    % plot slope
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+xbuff ynow xsiz ysiz]);
        ah = add_panel_title('a',sprintf('Example trajectory'),'yoffset',0,'xoffset',0,'width',400);  

        grad = atand(gradient(height) ./ gradient((1:c)./c.*3000));
        grad_mat = repmat(grad,r,1);

        % plot data
        imagesc('XData',[0 3],'YData',[0 1.5],'CData',grad_mat); hold on;

        % axis settings
        daspect([1 1 1])
        axis off
        colormap(ax2,cmap)
        ax2.CLim = [-55 55];
        xlabel('X (m)')
        ylabel('Y (m)')
        ytickformat('%.1f')
        xtickformat('%.1f')

        % additional plots
        ps = linspace(ax1.XLim(1),ax1.XLim(2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k--')  














        return
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig SXa.png']; 
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







