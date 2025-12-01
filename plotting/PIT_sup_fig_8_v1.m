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
    fs = [15 10];

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Behaviour anisotropy example trajectory
    xnow = -5;
    ynow = 720;
    xbuff = 210;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+50 220 90]);
        ah = add_panel_title('A',sprintf('Example trajectory'),'yoffset',-10,'xoffset',45,'width',400,'fontsize',fs);  

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
        lw = 1.5;
        set(c,'LineWidth',lw);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        cm = cmocean('balance');
        % cm = cmocean('thermal');        
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        colormap(ax1,cmap)
        view(0,90)
        ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax1.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax1.CLim = [-180 180];

        % text
        text(0,-0.02,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color','k','VerticalAlignment','top','rotation',0)
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
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
        text(0,-0.02,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color','k','VerticalAlignment','top','rotation',0)
        
    % plot arena 2
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+xbuff ax2.Position(2) ax1.Position(3) ax1.Position(4)]);
        % get data
        pp = 3;
        idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
        pox = pos.pox_planar(idx);
        poy = pos.poy_planar(idx);
        poz = pos.poz_planar(idx);
        poh = pos.yaw(idx);

        % plot data
        plot3(mf(:,1),mf(:,2),mf(:,3),'k'); hold on;
        c = cline(pox,poy,poz,poh);
        set(c,'LineWidth',lw);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        cm = cmocean('balance');
        % cm = cmocean('thermal');        
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        colormap(ax3,cmap)
        view(0,90)
        ax3.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax3.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax3.CLim = [-180 180];

        % text
        text(0,-0.02,sprintf('%s',maze_names{3}),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color','k','VerticalAlignment','top','rotation',0)
        
    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax2.Position(1)+20 ax1.Position(2)+100 80 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Azimuth (%c)',176),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

%% check 3D heading maps exist
    [clumaa,posdata] = PIT_3D_heading(config,clumaa,posdata,0);

%% >>>>>>>>>> Directional sampling in each maze, 3D HD heatmaps
    xnow = 50;
    ynow = ynow-100;
    ybuff = 150;
    xsiz = 150;
    ysiz = 90;
    xbuff = 50;

    % plot arena 3D HD map
    ax1b = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('B',sprintf('Movement direction - Example session'),'yoffset',-10,'xoffset',0,'width',400,'fontsize',fs);  
    
        % get data
        % dwellmap = posdata.hd_3d_dwell{eg_session,1}; % HD
        dwellmap = posdata.mov_3d_dwell{eg_session,1}; % movement 

        % plot data
        imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))

        % axis settings
        xlabel(sprintf('Azimuth (%c)',176))
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        view(0,90);
        ax1b.FontSize = 8;
        ax1b.CLim = [0,max(dwellmap(:))];
        colormap(gca,turbo)
        ax1b.XTick = -180:90:180;
        ax1b.YTick = -90:45:90;

    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax1b.Position(1)+260 ax1b.Position(2)+95 80 8]);
        x = linspace(ax1b.CLim(1),ax1b.CLim(2),100);
        imagesc(x,'XData',ax1b.CLim,'YData',[0 1]);
        colormap(axc,ax1b.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('0 '),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' Max'),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Dwell time'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

    % plot hills 3D HD map
    ax2b = axes('Units','pixels','Position',[ax1b.Position(1)+ax1b.Position(3)+xbuff ynow xsiz ysiz]);
        % get data
        % dwellmap = posdata.hd_3d_dwell{eg_session,2}; % HD
        dwellmap = posdata.mov_3d_dwell{eg_session,2}; % movement

        % plot data
        imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))

        % axis settings
        xlabel(sprintf('Azimuth (%c)',176))
        ylabel(sprintf('Tilt (%c)',176))        
        axis xy tight
        view(0,90);
        ax2b.FontSize = 8;
        ax2b.CLim = [0,max(dwellmap(:))];
        colormap(gca,turbo)
        ax2b.XTick = ax1b.XTick;
        ax2b.YTick = -90:45:90;       

    % plot hills 3D HD map
    ax3b = axes('Units','pixels','Position',[ax2b.Position(1)+ax2b.Position(3)+xbuff ynow xsiz ysiz]);
        % get data
        % dwellmap = posdata.hd_3d_dwell{eg_session,2}; % HD
        dwellmap = posdata.mov_3d_dwell{eg_session,3}; % movement

        % plot data
        imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))

        % axis settings
        xlabel(sprintf('Azimuth (%c)',176))
        ylabel(sprintf('Tilt (%c)',176))        
        axis xy tight
        view(0,90);
        ax3b.FontSize = 8;
        ax3b.CLim = [0,max(dwellmap(:))];
        colormap(gca,turbo)
        ax3b.XTick = ax1b.XTick;
        ax3b.YTick = -90:45:90;       

%% >>>>>>>>>> Directional sampling in each maze, averaged
    xnow = xnow;
    ynow = ynow-155;
    ysiz2 = 60;
    ybuff2 = 10;

    % plot arena 3D HD map (averaged)
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('C',sprintf('Averaged across sessions'),'yoffset',-10,'xoffset',0,'width',400,'fontsize',fs);  

        % get data
        % dwellmap = posdata.hd_3d_dwell(:,1); % HD
        dwellmap = posdata.mov_3d_dwell(:,1); % movement
        dwellmap = cat(3,dwellmap{:});
        for ii = 1:size(dwellmap,3)
            dwellmap(:,:,ii) = (dwellmap(:,:,ii) - mean(dwellmap(:,:,ii),'all','omitmissing')) ./ std(dwellmap(:,:,ii),[],'all','omitmissing');
        end
        dwellmap_a = mean(dwellmap,3,'omitmissing');

        % plot data
        imagesc([-180 180],[-90 90],dwellmap_a,'alphadata',~isnan(dwellmap_a))

        % axis settings
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        view(0,90);
        ax1.FontSize = 8;
        ax1.CLim = [-1,6];
        colormap(gca,turbo)
        ax1.XTick = -180:90:180;
        ax1.YTick = -90:45:90;
        xlabel(sprintf('Azimuth (%c)',176))

        % text
        text(0,1,sprintf('N = %d sessions',size(dwellmap,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','w')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax1.Position(1)+260 ax1.Position(2)+95 80 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Dwell time (z)'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

        % % plot yaw only
        % ax1c = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-ysiz2-ybuff2 ax1.Position(3) ysiz2]); 
        %     dwellmap = posdata.hd_3d_dwell(:,1); 
        %     v1 = cat(3,dwellmap{:});
        %     vm = squeeze(sum(v1,1,'omitmissing'))';
        %     vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
        % 
        %     x = linspace(-180,180,size(vm,2));
        %     m = mean(vm,1,'omitmissing');
        %     s = std(vm,[],1,'omitmissing');
        % 
        %     b = boundedline(x,m,s,'k','cmap',rgb('Black'));
        % 
        %     % axis settings
        %     ax1c.XTick = -180:90:180;
        %     ax1c.XLim = [-180 180];
        %     xlabel(sprintf('Azimuth (%c)',176))
        %     ax1c.YLim = [-2 3];
        %     ylabel(sprintf('Time (z)'))
        % 
        %     % additional plots
        %     line(ax1c.XLim,[0 0],'Color',[.5 .5 .5]);
        %     line([0 0],ax1c.YLim,'Color','k');
        %     line([90 90],ax1c.YLim,'Color','k','LineStyle',':');
        %     line([-90 -90],ax1c.YLim,'Color','k','LineStyle',':');   

    % plot hills 3D HD map (averaged)
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz ysiz]);
        % get data
        % dwellmap = posdata.hd_3d_dwell(:,2); % HD
        dwellmap = posdata.mov_3d_dwell(:,2); % movement        
        dwellmap = cat(3,dwellmap{:});
        for ii = 1:size(dwellmap,3)
            dwellmap(:,:,ii) = (dwellmap(:,:,ii) - mean(dwellmap(:,:,ii),'all','omitmissing')) ./ std(dwellmap(:,:,ii),[],'all','omitmissing');
        end
        dwellmap_a = mean(dwellmap,3,'omitmissing'); 

        % plot data
        imagesc([-180 180],[-90 90],dwellmap_a,'alphadata',~isnan(dwellmap_a))

        % axis settings
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        view(0,90);
        ax2.FontSize = 8;
        ax2.CLim = ax1.CLim;
        colormap(gca,turbo)
        ax2.XTick = -180:90:180;
        ax2.YTick = ax1.YTick;        
        xlabel(sprintf('Azimuth (%c)',176))

        % text
        text(0,1,sprintf('N = %d sessions',size(dwellmap,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','w')

        % % plot yaw only
        % ax2c = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)-ysiz2-ybuff2 ax2.Position(3) ysiz2]); 
        %     dwellmap = posdata.hd_3d_dwell(:,2); 
        %     v1 = cat(3,dwellmap{:});
        %     vm = squeeze(sum(v1,1,'omitmissing'))';
        %     vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
        % 
        %     x = linspace(-180,180,size(vm,2));
        %     m = mean(vm,1,'omitmissing');
        %     s = std(vm,[],1,'omitmissing');
        % 
        %     b = boundedline(x,m,s,'k','cmap',rgb('Black'));
        % 
        %     % axis settings
        %     ax2c.XTick = -180:90:180;
        %     ax2c.XLim = [-180 180];
        %     xlabel(sprintf('Azimuth (%c)',176))
        %     ax2c.YLim = ax1c.YLim;
        %     ax2c.YTick = [];
        %     % ylabel(sprintf('Time (z)'))
        % 
        %     % additional plots
        %     line(ax2c.XLim,[0 0],'Color',[.5 .5 .5]);
        %     line([0 0],ax2c.YLim,'Color','k');
        %     line([90 90],ax2c.YLim,'Color','k','LineStyle',':');
        %     line([-90 -90],ax2c.YLim,'Color','k','LineStyle',':'); 

    % plot arena 3D HD map (averaged)
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+xbuff ynow xsiz ysiz]);
        % get data
        % dwellmap = posdata.hd_3d_dwell(:,1); % HD
        dwellmap = posdata.mov_3d_dwell(:,3); % movement
        dwellmap = cat(3,dwellmap{:});
        for ii = 1:size(dwellmap,3)
            dwellmap(:,:,ii) = (dwellmap(:,:,ii) - mean(dwellmap(:,:,ii),'all','omitmissing')) ./ std(dwellmap(:,:,ii),[],'all','omitmissing');
        end
        dwellmap_a = mean(dwellmap,3,'omitmissing');

        % plot data
        imagesc([-180 180],[-90 90],dwellmap_a,'alphadata',~isnan(dwellmap_a))

        % axis settings
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        view(0,90);
        ax3.FontSize = 8;
        ax3.CLim = [-1,6];
        colormap(gca,turbo)
        ax3.XTick = -180:90:180;
        ax3.YTick = -90:45:90;
        xlabel(sprintf('Azimuth (%c)',176))

        % text
        text(0,1,sprintf('N = %d sessions',size(dwellmap,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','w')

        % % plot yaw only
        % ax1c = axes('Units','pixels','Position',[ax3.Position(1) ax3.Position(2)-ysiz2-ybuff2 ax3.Position(3) ysiz2]); 
        %     dwellmap = posdata.hd_3d_dwell(:,1); 
        %     v1 = cat(3,dwellmap{:});
        %     vm = squeeze(sum(v1,1,'omitmissing'))';
        %     vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
        % 
        %     x = linspace(-180,180,size(vm,2));
        %     m = mean(vm,1,'omitmissing');
        %     s = std(vm,[],1,'omitmissing');
        % 
        %     b = boundedline(x,m,s,'k','cmap',rgb('Black'));
        % 
        %     % axis settings
        %     ax1c.XTick = -180:90:180;
        %     ax1c.XLim = [-180 180];
        %     xlabel(sprintf('Azimuth (%c)',176))
        %     ax1c.YLim = [-2 3];
        %     ylabel(sprintf('Time (z)'))
        % 
        %     % additional plots
        %     line(ax1c.XLim,[0 0],'Color',[.5 .5 .5]);
        %     line([0 0],ax1c.YLim,'Color','k');
        %     line([90 90],ax1c.YLim,'Color','k','LineStyle',':');
        %     line([-90 -90],ax1c.YLim,'Color','k','LineStyle',':');  

%% >>>>>>>>>> Directional sampling in each maze, 3D HD heatmaps
    xnow = 50;
    ynow = ynow-170;
    ybuff = 150;
    xsiz = 150;
    ysiz = 90;
    xbuff = 50;

    % plot arena 3D HD map
    ax1b = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('D',sprintf('Head direction - Example session'),'yoffset',-10,'xoffset',0,'width',400,'fontsize',fs);  
    
        % get data
        dwellmap = posdata.hd_3d_dwell{eg_session,1}; % HD
        % dwellmap = posdata.mov_3d_dwell{eg_session,1}; % movement 

        % plot data
        imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))

        % axis settings
        xlabel(sprintf('Azimuth (%c)',176))
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        view(0,90);
        ax1b.FontSize = 8;
        ax1b.CLim = [0,max(dwellmap(:))];
        colormap(gca,turbo)
        ax1b.XTick = -180:90:180;
        ax1b.YTick = -90:45:90;

    % plot hills 3D HD map
    ax2b = axes('Units','pixels','Position',[ax1b.Position(1)+ax1b.Position(3)+xbuff ynow xsiz ysiz]);
        % get data
        dwellmap = posdata.hd_3d_dwell{eg_session,2}; % HD
        % dwellmap = posdata.mov_3d_dwell{eg_session,2}; % movement

        % plot data
        imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))

        % axis settings
        xlabel(sprintf('Azimuth (%c)',176))
        ylabel(sprintf('Tilt (%c)',176))        
        axis xy tight
        view(0,90);
        ax2b.FontSize = 8;
        ax2b.CLim = [0,max(dwellmap(:))];
        colormap(gca,turbo)
        ax2b.XTick = ax1b.XTick;
        ax2b.YTick = -90:45:90;       

    % plot hills 3D HD map
    ax3b = axes('Units','pixels','Position',[ax2b.Position(1)+ax2b.Position(3)+xbuff ynow xsiz ysiz]);
        % get data
        dwellmap = posdata.hd_3d_dwell{eg_session,3}; % HD
        % dwellmap = posdata.mov_3d_dwell{eg_session,3}; % movement

        % plot data
        imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))

        % axis settings
        xlabel(sprintf('Azimuth (%c)',176))
        ylabel(sprintf('Tilt (%c)',176))        
        axis xy tight
        view(0,90);
        ax3b.FontSize = 8;
        ax3b.CLim = [0,max(dwellmap(:))];
        colormap(gca,turbo)
        ax3b.XTick = ax1b.XTick;
        ax3b.YTick = -90:45:90;    

    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax2b.Position(1)+60 ax1b.Position(2)+100 80 8]);
        x = linspace(ax1b.CLim(1),ax1b.CLim(2),100);
        imagesc(x,'XData',ax1b.CLim,'YData',[0 1]);
        colormap(axc,ax1b.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.f ',ax1b.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.f',ax1b.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Dwell time'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> Directional sampling in each maze, averaged
    xnow = xnow;
    ynow = ynow-155;
    ysiz2 = 60;
    ybuff2 = 10;

    % plot arena 3D HD map (averaged)
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('E',sprintf('Averaged across sessions'),'yoffset',-10,'xoffset',0,'width',400,'fontsize',fs);  

        % get data
        dwellmap = posdata.hd_3d_dwell(:,1); % HD
        % dwellmap = posdata.mov_3d_dwell(:,1); % movement
        dwellmap = cat(3,dwellmap{:});
        for ii = 1:size(dwellmap,3)
            dwellmap(:,:,ii) = (dwellmap(:,:,ii) - mean(dwellmap(:,:,ii),'all','omitmissing')) ./ std(dwellmap(:,:,ii),[],'all','omitmissing');
        end
        dwellmap_a = mean(dwellmap,3,'omitmissing');

        % plot data
        imagesc([-180 180],[-90 90],dwellmap_a,'alphadata',~isnan(dwellmap_a))

        % axis settings
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        view(0,90);
        ax1.FontSize = 8;
        ax1.CLim = [-1,3];
        colormap(gca,turbo)
        ax1.XTick = -180:90:180;
        ax1.YTick = -90:45:90;
        xlabel(sprintf('Azimuth (%c)',176))

        % text
        text(0,1,sprintf('N = %d sessions',size(dwellmap,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','w')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax1.Position(1)+260 ax1.Position(2)+95 80 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Dwell time (z)'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

        % % plot yaw only
        % ax1c = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-ysiz2-ybuff2 ax1.Position(3) ysiz2]); 
        %     dwellmap = posdata.hd_3d_dwell(:,1); 
        %     v1 = cat(3,dwellmap{:});
        %     vm = squeeze(sum(v1,1,'omitmissing'))';
        %     vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
        % 
        %     x = linspace(-180,180,size(vm,2));
        %     m = mean(vm,1,'omitmissing');
        %     s = std(vm,[],1,'omitmissing');
        % 
        %     b = boundedline(x,m,s,'k','cmap',rgb('Black'));
        % 
        %     % axis settings
        %     ax1c.XTick = -180:90:180;
        %     ax1c.XLim = [-180 180];
        %     xlabel(sprintf('Azimuth (%c)',176))
        %     ax1c.YLim = [-2 3];
        %     ylabel(sprintf('Time (z)'))
        % 
        %     % additional plots
        %     line(ax1c.XLim,[0 0],'Color',[.5 .5 .5]);
        %     line([0 0],ax1c.YLim,'Color','k');
        %     line([90 90],ax1c.YLim,'Color','k','LineStyle',':');
        %     line([-90 -90],ax1c.YLim,'Color','k','LineStyle',':');   

    % plot hills 3D HD map (averaged)
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz ysiz]);
        % get data
        dwellmap = posdata.hd_3d_dwell(:,2); % HD
        % dwellmap = posdata.mov_3d_dwell(:,2); % movement        
        dwellmap = cat(3,dwellmap{:});
        for ii = 1:size(dwellmap,3)
            dwellmap(:,:,ii) = (dwellmap(:,:,ii) - mean(dwellmap(:,:,ii),'all','omitmissing')) ./ std(dwellmap(:,:,ii),[],'all','omitmissing');
        end
        dwellmap_a = mean(dwellmap,3,'omitmissing'); 

        % plot data
        imagesc([-180 180],[-90 90],dwellmap_a,'alphadata',~isnan(dwellmap_a))

        % axis settings
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        view(0,90);
        ax2.FontSize = 8;
        ax2.CLim = ax1.CLim;
        colormap(gca,turbo)
        ax2.XTick = -180:90:180;
        ax2.YTick = ax1.YTick;        
        xlabel(sprintf('Azimuth (%c)',176))

        % text
        text(0,1,sprintf('N = %d sessions',size(dwellmap,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','w')

        % % plot yaw only
        % ax2c = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)-ysiz2-ybuff2 ax2.Position(3) ysiz2]); 
        %     dwellmap = posdata.hd_3d_dwell(:,2); 
        %     v1 = cat(3,dwellmap{:});
        %     vm = squeeze(sum(v1,1,'omitmissing'))';
        %     vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
        % 
        %     x = linspace(-180,180,size(vm,2));
        %     m = mean(vm,1,'omitmissing');
        %     s = std(vm,[],1,'omitmissing');
        % 
        %     b = boundedline(x,m,s,'k','cmap',rgb('Black'));
        % 
        %     % axis settings
        %     ax2c.XTick = -180:90:180;
        %     ax2c.XLim = [-180 180];
        %     xlabel(sprintf('Azimuth (%c)',176))
        %     ax2c.YLim = ax1c.YLim;
        %     ax2c.YTick = [];
        %     % ylabel(sprintf('Time (z)'))
        % 
        %     % additional plots
        %     line(ax2c.XLim,[0 0],'Color',[.5 .5 .5]);
        %     line([0 0],ax2c.YLim,'Color','k');
        %     line([90 90],ax2c.YLim,'Color','k','LineStyle',':');
        %     line([-90 -90],ax2c.YLim,'Color','k','LineStyle',':'); 

    % plot arena 3D HD map (averaged)
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+xbuff ynow xsiz ysiz]);
        % get data
        dwellmap = posdata.hd_3d_dwell(:,3); % HD
        % dwellmap = posdata.mov_3d_dwell(:,3); % movement
        dwellmap = cat(3,dwellmap{:});
        for ii = 1:size(dwellmap,3)
            dwellmap(:,:,ii) = (dwellmap(:,:,ii) - mean(dwellmap(:,:,ii),'all','omitmissing')) ./ std(dwellmap(:,:,ii),[],'all','omitmissing');
        end
        dwellmap_a = mean(dwellmap,3,'omitmissing');

        % plot data
        imagesc([-180 180],[-90 90],dwellmap_a,'alphadata',~isnan(dwellmap_a))

        % axis settings
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        view(0,90);
        ax3.FontSize = 8;
        ax3.CLim = [-1,6];
        colormap(gca,turbo)
        ax3.XTick = -180:90:180;
        ax3.YTick = -90:45:90;
        xlabel(sprintf('Azimuth (%c)',176))

        % text
        text(0,1,sprintf('N = %d sessions',size(dwellmap,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','w')

        % % plot yaw only
        % ax1c = axes('Units','pixels','Position',[ax3.Position(1) ax3.Position(2)-ysiz2-ybuff2 ax3.Position(3) ysiz2]); 
        %     dwellmap = posdata.hd_3d_dwell(:,1); 
        %     v1 = cat(3,dwellmap{:});
        %     vm = squeeze(sum(v1,1,'omitmissing'))';
        %     vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
        % 
        %     x = linspace(-180,180,size(vm,2));
        %     m = mean(vm,1,'omitmissing');
        %     s = std(vm,[],1,'omitmissing');
        % 
        %     b = boundedline(x,m,s,'k','cmap',rgb('Black'));
        % 
        %     % axis settings
        %     ax1c.XTick = -180:90:180;
        %     ax1c.XLim = [-180 180];
        %     xlabel(sprintf('Azimuth (%c)',176))
        %     ax1c.YLim = [-2 3];
        %     ylabel(sprintf('Time (z)'))
        % 
        %     % additional plots
        %     line(ax1c.XLim,[0 0],'Color',[.5 .5 .5]);
        %     line([0 0],ax1c.YLim,'Color','k');
        %     line([90 90],ax1c.YLim,'Color','k','LineStyle',':');
        %     line([-90 -90],ax1c.YLim,'Color','k','LineStyle',':');  

        % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S7.png']; 
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







