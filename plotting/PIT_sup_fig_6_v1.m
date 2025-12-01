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

% distance covered vs session duration
% occupied bins
% velocity



%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

%% check 3D heading maps exist
    [clumaa,posdata] = PIT_3D_heading(config,clumaa,posdata,0);

%% >>>>>>>>>> Directional sampling in each maze
    xnow = 50;
    ynow = 730;
    xbuff = 50;
    ybuff = -40;
    xsiz = 200;
    ysiz = 150;

    % plot example arena session
    ax1a = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('a',sprintf('Behaviour'),'yoffset',-30,'xoffset',0,'width',400);  

        % get data
        eg_session = 3;
        pos = posdata.pos{eg_session};   
        stimes = posdata.session_times{eg_session};
        mf = posdata.maze_frame{eg_session};

        pp = 1;
        idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
        pox = pos.pox_surficial(idx);
        poy = pos.poy_surficial(idx);
        poz = pos.poz_surficial(idx);
        poh = pos.yaw(idx);
        mf = posdata.maze_frame{eg_session};

        % plot data
        plot3(mf(:,1),mf(:,2),mf(:,3),'k'); hold on;        
        c = cline(pox,poy,poz,poh);
        set(c,'LineWidth',1);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        cm = cmocean('balance');
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        colormap(ax1a,cmap)
        view(0,90)
        ax1a.XLim = [min(mf(:,4)) max(mf(:,4))];
        ax1a.YLim = [min(mf(:,2)) max(mf(:,2))];

        % additional plots
        ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax1a.YLim,'k:')

        % text
        text(0,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

    % plot arena
    ax1b = axes('Units','pixels','Position',[xnow ax1a.Position(2)-ysiz-ybuff xsiz ysiz]);
        % get data
        v1 = posdata.direction_3D(:,1); % arena 1
        f = @(x) (x - mean(x,"all",'omitmissing')) ./ std(x,[],"all",'omitmissing');
        % f = @(x) x ./ sum(x,'all','omitmissing');
        v1 = cellfun(f,v1,'UniformOutput',false);        
        v1 = cat(3,v1{:});
        mean_map = mean(v1,3,'omitnan');

        % plot data
        imagesc([-180 180],[-90 90],mean_map,'alphadata',~isnan(mean_map)); hold on;
        line(ax1b.XLim,[0 0],'Color',[.5 .5 .5]);
        line([0 0],ax1b.YLim,'Color','w');
        line([90 90],ax1b.YLim,'Color','w','LineStyle',':');
        line([-90 -90],ax1b.YLim,'Color','w','LineStyle',':');

        % axis settings
        ylabel(sprintf('Tilt (%c)',176))
        axis xy tight
        daspect([1 1 1])
        view(0,90);
        ax1b.CLim = [-1 5];
        colormap(gca,turbo)
        ax1b.XTick = [];
        ax1b.YTick = -90:45:90;    

    % plot yaw only
    ax1c = axes('Units','pixels','Position',[xnow ax1b.Position(2)-ysiz+115 xsiz 45]);
        v1 = posdata.direction_3D(:,1); % arena 1    
        v1 = cat(3,v1{:});        
        vm = squeeze(sum(v1,1,'omitmissing'))';
        vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');
        x = linspace(-180,180,size(vm,2));
        m = mean(vm,1,'omitmissing');
        s = std(vm,[],1,'omitmissing');

        b = boundedline(x,m,s,'k','cmap',rgb('Black'));

        % axis settings
        ax1c.XTick = -180:90:180;
        ax1c.XLim = [-180 180];
        xlabel(sprintf('Azimuth (%c)',176))
        ax1c.YLim = [-2 3];
        ylabel(sprintf('Time (z)'))

        % additional plots
        line(ax1c.XLim,[0 0],'Color',[.5 .5 .5]);
        line([0 0],ax1c.YLim,'Color','k');
        line([90 90],ax1c.YLim,'Color','k','LineStyle',':');
        line([-90 -90],ax1c.YLim,'Color','k','LineStyle',':');        

    % plot example hills session
    xnow = xnow+240;
    ax2a = axes('Units','pixels','Position',[ax1a.Position(1)+250 ax1a.Position(2) ax1a.Position(3) ax1a.Position(4)]);
        pp = 2;
        idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
        pox = pos.pox_planar(idx);
        poy = pos.poy_planar(idx);
        poz = pos.poz_planar(idx);
        poh = pos.yaw(idx);

        % plot data
        plot3(mf(:,1),mf(:,2),mf(:,3),'k'); hold on;        
        c = cline(pox,poy,poz,poh);
        set(c,'LineWidth',1);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        cm = cmocean('balance');
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        colormap(ax2a,cmap)
        view(0,90)
        ax2a.XLim = [min(mf(:,4)) max(mf(:,4))];
        ax2a.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax2a.YTick = [];

        % additional plots
        ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2a.YLim,'k:')

        % text
        text(0,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

    % colorbar
    axc = axes('Units','pixels','Position',[ax2a.Position(1)+ax2a.Position(3)-20 ax2a.Position(2)+40 10 70]);
        mat = (linspace(0,100,100))';
        imagesc(mat,ones(size(mat)),mat);
        colormap(axc,ax2a.Colormap);
        axis xy

        axc.YTick = [];
        axc.XTick = [];
        text(0,1.25,sprintf('Direction (%c)',176),'FontSize',7,'HorizontalAl','left','Units','normalized')
        text(0.5,0,sprintf('-180'),'FontSize',7,'HorizontalAl','center','Units','normalized','VerticalAl','top')
        text(0.5,1,sprintf('180'),'FontSize',7,'HorizontalAl','center','Units','normalized','VerticalAl','bottom')        

    % plot hills
    ax2b = axes('Units','pixels','Position',[ax2a.Position(1) ax2a.Position(2)-ysiz-ybuff ax1b.Position(3) ax1b.Position(4)]);
        % get data
        v1 = posdata.direction_3D(:,2); % arena 1
        v1 = cellfun(f,v1,'UniformOutput',false);
        v1 = cat(3,v1{:});
        mean_map = mean(v1,3,'omitnan');

        % plot data
        imagesc([-180 180],[-90 90],mean_map,'alphadata',~isnan(mean_map)); hold on;
        line(ax1b.XLim,[0 0],'Color',[.5 .5 .5]);
        line([0 0],ax1b.YLim,'Color','w');
        line([90 90],ax1b.YLim,'Color','w','LineStyle',':');
        line([-90 -90],ax1b.YLim,'Color','w','LineStyle',':');

        % axis settings
        axis xy tight
        daspect([1 1 1])
        view(0,90);
        ax2b.CLim = ax1b.CLim;
        colormap(gca,turbo)
        ax2b.XTick = [];
        ax2b.YTick = -90:45:90;    

    axc = axes('Units','pixels','Position',[ax2b.Position(1)+ax1a.Position(3)+10 ax2b.Position(2)+30 10 70]);
        mat = (linspace(0,100,100))';
        imagesc(mat,ones(size(mat)),mat);
        colormap(axc,ax1b.Colormap);
        axis xy

        axc.YTick = [];
        axc.XTick = [];
        text(0,1.25,sprintf('Time (z)'),'FontSize',7,'HorizontalAl','left','Units','normalized')
        text(0.5,1,sprintf('%.1f',ax1b.CLim(2)),'FontSize',7,'HorizontalAl','center','Units','normalized','VerticalAl','bottom')
        text(0.5,0,sprintf('%.1f',ax1b.CLim(1)),'FontSize',7,'HorizontalAl','center','Units','normalized','VerticalAl','top')

    % plot yaw only
    ax2c = axes('Units','pixels','Position',[ax2a.Position(1) ax1c.Position(2) ax1c.Position(3) ax1c.Position(4)]);
        v1 = posdata.direction_3D(:,2); % arena 1    
        v1 = cat(3,v1{:});        
        vm = squeeze(sum(v1,1,'omitmissing'))';
        vm = (vm - mean(vm,2,'omitmissing')) ./ std(vm,[],2,'omitmissing');        
        x = linspace(-180,180,size(vm,2));
        m = mean(vm,1,'omitmissing');
        s = std(vm,[],1,'omitmissing');

        b = boundedline(x,m,s,'k','cmap',rgb('Black'));
        line(ax1c.XLim,[0 0],'Color',[.5 .5 .5]);

        % axis settings
        ax2c.XTick = ax1c.XTick;
        ax2c.XLim = ax1c.XLim;        
        xlabel(sprintf('Azimuth (%c)',176))
        ax2c.YLim = ax1c.YLim;

        % additional plots
        line(ax1c.XLim,[0 0],'Color',[.5 .5 .5]);
        line([0 0],ax1c.YLim,'Color','k');
        line([90 90],ax1c.YLim,'Color','k','LineStyle',':');
        line([-90 -90],ax1c.YLim,'Color','k','LineStyle',':');  


        return





%% >>>>>>>>>> %% Quantify coverage
    mapset.ppm          = 1000;
    mapset.method       = 'histogram';
    mapset.binsize      = 100; % (mm) firing rate map bin size
    mapset.ssigma       = 64; % (mm) firing rate map smoothing sigma
    mapset.padding      = 0; % (mm) how much to pad the edges of the firing rate map
    mapset.mindwell     = 0; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    mapset.mindist      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.smethod      = 3; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing    
    mapset.frcut        = 0.2; % cutoff % for regions to be considered a place field
    mapset.arcut        = 36; % cutoff area for regions to be considered a place field
    mapset.drive_height_mm = 20;

    if 1%~any(ismember(posdata.Properties.VariableNames,'dwellmap')) % if the column(s) do not exist yet
        posdata.dwellmap = cell(size(posdata,1),3);

        for ss = 1:size(posdata,1) % for each recording session
            disp(sprintf('\tSession %d of %d (%.f%%)',ss,size(posdata,1),ss/size(posdata,1)*100))
            session_times = posdata.session_times{ss};
            pos = posdata.pos{ss};
    
            for pp = 1:size(session_times,1) % for each part
                pidxn = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part
                if pp==1 || pp==3
                    p3 = [pos.pox_planar(pidxn) pos.poy_planar(pidxn)];   
                    mf = posdata.maze_frame{eg_session};        
                    lx = mf(:,1); % in mm                
                    ly = mf(:,2); % in mm        
                elseif pp==2
                    p3 = [pos.pox_surficial(pidxn) pos.poy_surficial(pidxn)];   
                    mf = posdata.maze_frame{eg_session};        
                    lx = mf(:,4); % in mm                
                    ly = mf(:,2); % in mm  
                end
                rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
                
                % ratemap
                [~,dwellmap, ~, outputs,~] = rate_mapper(p3,[NaN NaN],rmset); 
                posdata.dwellmap(ss,pp) = { single(dwellmap) };

                % figure
                % imagesc(dwellmap,'alphadata',~isnan(dwellmap))
                % daspect([1 1 1])
                % keyboard
            end
        end
    end
  
%% >>>>>>>>>> Field anisotropy vs behaviour (case by case)
    xnow = 50;
    ynow = 500;

    % create axis
    ax1 = axes('Units','pixels','Position',[xnow ynow 200 160]);    
        ah = add_panel_title('a',sprintf(''),'yoffset',0,'xoffset',-70,'width',400);  

        pp = 1; % arena 1
        for ss = 1:size(posdata,1) % for each recording session
            m = posdata.dwellmap{ss,pp};
            [f,x] = ecdf(m(:));
            stairs(x,f,'Color',plot_set{1,1}); hold on;
        end

        ax1.XLim = [0 0.5];
        ax1.YLim = [0 1];


%% >>>>>>>>>> Save the overall figure
    if 0
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
            exportgraphics(gcf,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',350);  
        end
        close(gcf);   
    end







