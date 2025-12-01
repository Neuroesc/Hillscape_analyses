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
%% >>>>>>>>>> Hypothesis schematics
    xnow = 20;
    ynow = 580;
    ax = axes('Units','pixels','Position',[xnow ynow 590 400]); ax1 = ax;
        ah = add_panel_title('A',sprintf(''),'yoffset',-150,'xoffset',25,'width',400,'fontsize',fs);

        % plot image
        img = imread([config.main_dir '\associated media\hypotheses.jpg']);
        imshow(img);

        % text
        text((1/6)*1,0.95,sprintf('Volumetric'),'Units','normalized','FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','center','rotation',0)
        text((1/6)*2.9,0.95,sprintf('Earth-horizontal'),'Units','normalized','FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','center','rotation',0)
        text((1/6)*4.9,0.95,sprintf('Terrain/surface specific'),'Units','normalized','FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','center','rotation',0)

        yt = 0;
        text((1/6)*0+0.02,yt,sprintf('-Field surface area changes\n-Area change is dependent on field\n-New fields may be observed'),'Units','normalized','FontSize',7,'VerticalAlignment','top','HorizontalAlignment','left','rotation',0)
        text((1/6)*1.95+0.02,yt,sprintf('-Fields increase in surface area\n-From above fields appear unchanged\n-No new fields are observed'),'Units','normalized','FontSize',7,'VerticalAlignment','top','HorizontalAlignment','left','rotation',0)        
        text((1/6)*3.9+0.02,yt,sprintf('-Average field surface area is unchanged\n-Terrain changes result in remapping\n-Distal fields may also change'),'Units','normalized','FontSize',7,'VerticalAlignment','top','HorizontalAlignment','left','rotation',0)

        % axis settings
        daspect([1 1 1])
        axis off
% keyboard
%% >>>>>>>>>> Schematics of mazes
%% arena
    xnow = 30;
    ynow = 530;
    ax = axes('Units','pixels','Position',[xnow ynow 135 135]); ax1 = ax;
        ah = add_panel_title('B',sprintf(''),'yoffset',-50,'xoffset',15,'width',400,'fontsize',fs);
    
        PIT_plot_mazes(ax,1,2)
        title('Arena','FontSize',10,'FontWeight','normal')

        % additional text
        text(0.67,0.21,0,sprintf('3m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','center','rotation',20)
        text(0.12,0.05,0.5,sprintf('1.5m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','right','rotation',-34)
        text(0.7,0.85,0,sprintf('0.65m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','right')

    ax = axes('Units','pixels','Position',[xnow+ax1.Position(3)+10 ax1.Position(2) ax1.Position(3) ax1.Position(4)]); ax2 = ax;
        % indices
        pid = 2;
        p = posdata.pos{pid}; % [session, pot, pox, poy, poh, rx, ry, gx, gy, bx, by, poz, pox_planar, poy_planar, poz_planar, pox_surficial, poy_surficial, poz_surficial, poz_curve]
        sidx = p.session==1;

        % plot data
        c = cline(p.pox(sidx),p.poy(sidx),-p.poz(sidx),-p.poz(sidx));
        set(c,'LineWidth',2)

        % axis settings
        daspect([1 1 1])
        view(3)
        colormap(flipud(magma(128)))
        axis off tight
        ax.CLim = [100 500];
        camlight(ax)
        ax.Projection = 'perspective';

        % additional plots
        m = posdata.maze_frame{pid};
        pc = plotcube([range(m(:,1)) range(m(:,2)) 600],[min(m(:,1)) min(m(:,2)) min(m(:,3))],0,'k');
        ax.XLim = [min(m(:,1)) max(m(:,1))];
        ax.YLim = [min(m(:,2)) max(m(:,2))];
        ax.ZLim = [min(m(:,3)) max(m(:,3))];
        title('Arena','FontSize',10,'FontWeight','normal')

    axc = axes('Units','pixels','Position',[ax.Position(1)+70 ax.Position(2)+10 80 8]);
        mat = (linspace(0,100,100));
        imagesc(ones(size(mat)),mat,mat);
        colormap(axc,ax.Colormap);
        axis xy off

        axc.YTick = [];
        axc.XTick = [];
        text(0.5,1.7,sprintf('Height (mm)'),'FontSize',7,'HorizontalAl','center','Units','normalized')
        text(1.02,0.5,sprintf('500'),'FontSize',7,'HorizontalAl','left','Units','normalized')
        text(-0.02,0.5,sprintf('0'),'FontSize',7,'HorizontalAl','right','Units','normalized')

%% hills
    xnow = xnow+300;
    ax = axes('Units','pixels','Position',[xnow ax1.Position(2) ax1.Position(3) ax1.Position(4)]); ax1 = ax;
        PIT_plot_mazes(ax,2,2)
        title(maze_names{2},'FontSize',10,'FontWeight','normal')

        % additional text
        % text(0.67,0.18,0,sprintf('3m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','center','rotation',20)
        % text(1.08,0.95,0,sprintf('0.45m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','right')
        % text(0.12,0.05,0.5,sprintf('1.5m'),'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','right','rotation',-34)
% keyboard
    ax = axes('Units','pixels','Position',[xnow+ax1.Position(3)+10 ax1.Position(2) ax1.Position(3) ax1.Position(4)]); ax2 = ax;
        sidx = p.session==2;

        % plot data
        c = cline(p.pox(sidx),p.poy(sidx),-p.poz(sidx),-p.poz(sidx));
        set(c,'LineWidth',2)

        % axis settings
        daspect([1 1 1])
        view(3)
        colormap(flipud(magma(128)))
        axis off tight
        ax.CLim = [100 500];
        camlight(ax)
        ax.Projection = 'perspective';

        % additional plots
        m = posdata.maze_frame{pid};
        pc = plotcube([range(m(:,1)) range(m(:,2)) 600],[min(m(:,1)) min(m(:,2)) min(m(:,3))],0,'k');
        ax.XLim = [min(m(:,1)) max(m(:,1))];
        ax.YLim = [min(m(:,2)) max(m(:,2))];
        ax.ZLim = [min(m(:,3)) max(m(:,3))];
        title(maze_names{2},'FontSize',10,'FontWeight','normal')
% return
%% >>>>>>>>>> Schematic of experiment procedure
    xnow = 30;
    ynow = ynow-110;
    xbuff = 20;
    ax = axes('Units','pixels','Position',[xnow ynow 110 110],'Color','none'); ax1 = ax;
        ah = add_panel_title('C',sprintf(''),'yoffset',-50,'xoffset',15,'width',400,'fontsize',fs);
        PIT_plot_mazes(ax,1,2)

    ax = axes('Units','pixels','Position',ax1.Position+[ax1.Position(3)+xbuff 0 0 0],'Color','none'); ax2 = ax;    
        PIT_plot_mazes(ax,2,2)

    ax = axes('Units','pixels','Position',ax2.Position+[ax2.Position(3)+xbuff 0 0 0],'Color','none'); ax3 = ax;    
        PIT_plot_mazes(ax,1,2)

%% arrows etc
    ax = axes('Units','pixels','Position',[xnow ynow 500 ax1.Position(4)],'Color','none'); ax4 = ax;    
        ax.XLim = [0 1];
        ax.YLim = [0 1];
        ax.Color = 'none';
        axis manual off

        text(0.05,0.17,sprintf('%s\n(~25 mins)',maze_names{1}),'Units','normalized','FontSize',7,'VerticalAlignment','top','HorizontalAlignment','center')
        text(0.315,0.17,sprintf('%s\n(~45 mins)',maze_names{2}),'Units','normalized','FontSize',7,'VerticalAlignment','top','HorizontalAlignment','center')
        text(0.57,0.17,sprintf('%s\n(~25 mins)',maze_names{3}),'Units','normalized','FontSize',7,'VerticalAlignment','top','HorizontalAlignment','center')

        p1 = [0.12 0.2];
        p2 = [0.27 0.2];
        text(mean([p1(1) p2(1)]),0.17,sprintf('Home cage\n(~5 mins)'),'Units','normalized','FontSize',7,'VerticalAlignment','top','HorizontalAlignment','center')

        p1 = [0.38 0.2];
        p2 = [0.51 0.2];
        text(mean([p1(1) p2(1)]),0.17,sprintf('Home cage\n(~5 mins)'),'Units','normalized','FontSize',7,'VerticalAlignment','top','HorizontalAlignment','center')

        ya = 0.50;
        annotation('arrow', [0.15 0.25], [ya ya]); % x,y are (0–1)% x,y are normalized (0–1)
        annotation('arrow', [0.36 0.46], [ya ya]); % x,y are (0–1)% x,y are normalized (0–1)
% keyboard

%% >>>>>>>>>> Histology
xnow = xnow+400;
ynow = ynow-30;
    ax = axes('Units','pixels','Position',[xnow ynow 120 140]); ax1 = ax;
        ah = add_panel_title('D',sprintf(''),'yoffset',-30,'xoffset',0,'width',400,'fontsize',fs);

        iname = [config.data_out_dir 'Histology\RG40.tif'];
        [img1,m] = imread(iname,'tif');
        img1 = img1(:,:,1:3);
        imshow(img1(:,:,1:3)); hold on;

        xmin = 480;
        ymin = 880;
        ymax = 115;
        xmax = 850;
        plot([xmin xmax xmax xmin xmin],[ymax ymax ymin ymin ymax],'k','LineWidth',0.5);


    ax = axes('Units','pixels','Position',[xnow+ax1.Position(3)-20 ynow-30 70 ax1.Position(4)+20]); ax1 = ax;
        % plot image
        imshow(img1); hold on;

        plot(980,780,'k<','MarkerSize',5,'MarkerFaceColor','k');
        plot(1100,567,'k<','MarkerSize',5,'MarkerFaceColor','k');

        ax.XLim = [xmin xmax];
        ax.YLim = [ymax ymin];
        axis on
        box on
        ax.XTick = [];
        ax.YTick = [];

%% >>>>>>>>>> Example cells
    xnow = 20;
    ynow = ynow-115;
    xsiz = [90, 10]; % size, buffer
    ysiz = [70, 10]; % size, buffer
    xvec = xnow : sum(xsiz) : 1000;
    xvec = xvec(1:6);
    yvec = ynow : -sum(ysiz).*1.9 : -1000;
    yvec = yvec(1);
    [x,y] = ndgrid(xvec,yvec);
    ar_siz = [xsiz*0.48 ysiz*0.5];

    % find which cells to plot  
    % idx = find(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.1 & clumaa.repetition_score(:,1)<0.20 & clumaa.planar_spatial_info_shuffles(:,2)>30); % place cells in the hills
    % [~,sidx] = sort(clumaa.planar_spatial_info_shuffles(idx,1),'descend');
    % idx = idx(sidx);
    % idx = idx([1 2 5 10 13 14]);
    % ucis = {'RG19_221013_t1_c10','RG19_221026_t1_c4','RG26_230119_t1_c9','RG26_230119_t3_c1','RG11_220520_t2_c5','RG19_221014_t2_c9'};
    ucis = {'RG19_221013_t1_c10','RG19_221026_t1_c4','RG19_221013_t3_c13','RG26_230119_t3_c1','RG26_230119_t4_c5','RG19_221014_t2_c9'};

    for ii = 1:numel(x)
        % uci = clumaa.uci{idx(ii)};
        uci = ucis{ii};

        %% plot arena (spikes and positions)
        ax1 = axes('Units','pixels','Position',[x(ii) y(ii) xsiz(1) ysiz(1)]);
            if ii==1
                ah = add_panel_title('E',sprintf(''),'yoffset',0,'xoffset',20,'width',400,'fontsize',fs);                
                text(-0.02,0.97,maze_names{2},'Units','normalized','FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','left');
            end    

            % get data
            part_to_plot = 1;
            idxn = find( ismember(clumaa.uci,uci) & clumaa.partn==part_to_plot );
            pos_idx = clumaa.pos_idx(idxn);
            pos = posdata.pos{posdata.pos_idx==pos_idx};
            mf = posdata.maze_frame{posdata.pos_idx==pos_idx}; 

            session_times = posdata.session_times{posdata.pos_idx==pos_idx};
            pidxn = pos.pot > session_times(part_to_plot,1) & pos.pot < session_times(part_to_plot,2); % index for position data in this part

            spt = clumaa.spike_times_s{idxn};
            sidx = spt > session_times(part_to_plot,1) & spt < session_times(part_to_plot,2); % index for spike data in this part, first half
            spike_index = clumaa.spike_index{idxn};

            % plot data
            msiz = 4;
            lwidth = 0.5;
            pos_plot = plot(pos.pox_planar(pidxn),pos.poy_planar(pidxn),'Color',[.5 .5 .5],'LineWidth',lwidth); hold on;
            spk_plot = plot(pos.pox_planar(spike_index(sidx)),pos.poy_planar(spike_index(sidx)),'Marker','.','Color','r','MarkerSize',msiz,'LineStyle','none');
            plot(mf(:,1),mf(:,2),'Color','k');

            % axis settings
            daspect([1 1 1])
            axis off xy
            ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
            ax1.YLim = [min(mf(:,2)) max(mf(:,2))];

            % text
            text(0,1.05,sprintf('Cell %d',ii),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color','k','VerticalAlignment','bottom')
            if ii==1
                text(-0.02,0.5,sprintf(maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

                % legend
                [~,leg] = legendflex([pos_plot spk_plot],{'Positions','Spikes'},'anchor',{'e','e'},'ncol',2,'box','off','buffer',[90,55],'xscale',.5,'fontsize',9); 
                leg(3).LineWidth = 2;
                leg(6).MarkerSize = 18;
                leg(3).Color = pos_plot.Color;
                leg(6).Color = spk_plot.Color;

                % colorbar
                axc = axes('Units','pixels','Position',[ax1.Position(1)+250 ax1.Position(2)+85 80 8]);
                    mat = (linspace(0,100,100));
                    imagesc(ones(size(mat)),mat,mat);
                    colormap(axc,turbo);
                    axis xy off
            
                    axc.YTick = [];
                    axc.XTick = [];
                    text(0.5,1.7,sprintf('Firing rate (Hz)'),'FontSize',7,'HorizontalAl','center','Units','normalized')
                    text(1.02,0.5,sprintf('Max'),'FontSize',7,'HorizontalAl','left','Units','normalized')
                    text(-0.02,0.5,sprintf('0'),'FontSize',7,'HorizontalAl','right','Units','normalized')
            end

        %% plot hills (spikes and positions)
        ax2 = axes('Units','pixels','Position',[x(ii) ax1.Position(2)-ax1.Position(4)+20 xsiz(1) ysiz(1)]);
            % get data
            part_to_plot = 2;
            idxn = find( ismember(clumaa.uci,uci) & clumaa.partn==part_to_plot );
            pos_idx = clumaa.pos_idx(idxn);
            pos = posdata.pos{posdata.pos_idx==pos_idx};
            mf = posdata.maze_frame{posdata.pos_idx==pos_idx}; 

            session_times = posdata.session_times{posdata.pos_idx==pos_idx};
            pidxn = pos.pot > session_times(part_to_plot,1) & pos.pot < session_times(part_to_plot,2); % index for position data in this part

            spt = clumaa.spike_times_s{idxn};
            sidx = spt > session_times(part_to_plot,1) & spt < session_times(part_to_plot,2); % index for spike data in this part, first half
            spike_index = clumaa.spike_index{idxn};

            % plot data
            msiz = 4;
            lwidth = 0.5;
            pos_plot = plot(pos.pox_planar(pidxn),pos.poy_planar(pidxn),'Color',[.5 .5 .5],'LineWidth',lwidth); hold on;
            spk_plot = plot(pos.pox_planar(spike_index(sidx)),pos.poy_planar(spike_index(sidx)),'Marker','.','Color','r','MarkerSize',msiz,'LineStyle','none');
            plot(mf(:,1),mf(:,2),'Color','k');

            % axis settings
            daspect([1 1 1])
            axis off xy
            ax2.XLim = [min(mf(:,1)) max(mf(:,1))];
            ax2.YLim = [min(mf(:,2)) max(mf(:,2))];

            % additional plots
            ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax2.YLim,'k:')            

            % text
            if ii==1
                text(-0.02,0.5,sprintf(maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

        %% plot arena (ratemap)
        m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
        m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};
        m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};

        ax3 = axes('Units','pixels','Position',[x(ii) ax2.Position(2)-ax2.Position(4)+20 xsiz(1) ysiz(1)]);
            imagesc(m1,'alphadata',~isnan(m1)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax3,turbo)   
            ax3.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
            ax3.XTick = [];
            ax3.YTick = [];
            if ii==1
                text(-0.02,0.5,sprintf(maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

        %% plot hills (ratemap)
        ax4 = axes('Units','pixels','Position',[x(ii) ax3.Position(2)-ax3.Position(4)+20 xsiz(1) ysiz(1)]);
            imagesc(m2,'alphadata',~isnan(m2)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax4,turbo)   
            ax4.CLim = ax3.CLim;
            ax4.XTick = [];
            ax4.YTick = [];

            % additional plots
            ps = linspace(0,size(m2,2),7);
            tops = ps(2:2:end);
            valleys = ps(1:2:end);
            plot([tops; tops],ax4.YLim,'w:')
            
            % text
            text(0,-0.04,sprintf('%.1f',ax4.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','top')
            if ii==1
                text(-0.02,0.5,sprintf(maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end
    end

    keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig 1.png']; 
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


 


















