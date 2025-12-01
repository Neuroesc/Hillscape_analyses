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

    mapset.binsize = 32; % (mm) firing rate map bin siz

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example cells
    xnow = 20;
    ynow = 780;
    xsiz = [90, 10]; % size, buffer
    ysiz = [70, 10]; % size, buffer
    xvec = xnow : sum(xsiz) : 1000;
    xvec = xvec(1:6);
    yvec = ynow : -sum(ysiz).*1.9 : -1000;
    yvec = yvec(1);
    [x,y] = ndgrid(xvec,yvec);
    ar_siz = [xsiz*0.48 ysiz*0.5];

    % find which cells to plot  
    idx = find(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5 & clumaa.planar_spatial_info_shuffles(:,2)>30); % place cells in the hills
    [~,sidx] = sort(clumaa.repetition_score(idx,1),'ascend');
    idx = idx(sidx);
    idx = idx([2 3 5 7 8 9]);

    for ii = 1:numel(x)
        uci = clumaa.uci{idx(ii)};

        %% plot arena (spikes and positions)
        ax1 = axes('Units','pixels','Position',[x(ii) y(ii) xsiz(1) ysiz(1)]);
            if ii==1
                ah = add_panel_title('a',sprintf(''),'yoffset',0,'xoffset',20,'width',400,'fontsize',[18 13]);                
                text(0,0.97,maze_names{2},'Units','normalized','FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','left');
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
            pos_plot = plot(pos.pox_planar(pidxn),pos.poy_planar(pidxn),'Color','k','LineWidth',lwidth); hold on;
            spk_plot = plot(pos.pox_planar(spike_index(sidx)),pos.poy_planar(spike_index(sidx)),'Marker','.','Color','r','MarkerSize',msiz,'LineStyle','none');
            plot(mf(:,1),mf(:,2),'Color','k');

            % axis settings
            daspect([1 1 1])
            axis off xy
            ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
            ax1.YLim = [min(mf(:,2)) max(mf(:,2))];

            % text
            text(0,1.1,sprintf('Cell %d',ii),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color','k','VerticalAlignment','bottom')
            if ii==1
                text(0,0.5,sprintf('Arena 1'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

                % legend
                [~,leg] = legendflex([pos_plot spk_plot],{'Positions','Spikes'},'anchor',{'e','e'},'ncol',2,'box','off','buffer',[90,55],'xscale',.5,'fontsize',9); 
                leg(3).LineWidth = 2;
                leg(6).MarkerSize = 18;

                % colorbar
                axc = axes('Units','pixels','Position',[ax1.Position(1)+250 ax1.Position(2)+85 100 10]);
                    mat = (linspace(0,100,100));
                    imagesc(ones(size(mat)),mat,mat);
                    colormap(axc,turbo);
                    axis xy
            
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
            pos_plot = plot(pos.pox_planar(pidxn),pos.poy_planar(pidxn),'Color','k','LineWidth',lwidth); hold on;
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
            plot([tops; tops],ax2.YLim,'k')            

            % text
            if ii==1
                text(0,0.5,sprintf('Hills'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
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
                text(0,0.5,sprintf('Arena 1'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
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
            plot([tops; tops],ax4.YLim,'w')
            
            % text
            text(0,-0.02,sprintf('%.1f',ax4.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','top')
            if ii==1
                text(0,0.5,sprintf('Hills'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end
    end

%% >>>>>>>>>> repetition (field plots)
    % repetition scores
    var = 'repetition_score';
    rs_cutoff = -inf;
    r1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.repetition_score(:,1)>rs_cutoff,1); % arena 1 data
    r2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.repetition_score(:,1)>rs_cutoff,1); % hills data

    % ratemaps, averaged, zscored
    var = 'ratemap_planar';
    m2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.repetition_score(:,1)>rs_cutoff,1); % arena 1 data
    m2 = cellfun(@(x) imresize(x,[48 96],'bilinear'),m2,'UniformOutput',0);
    m2 = cat(3,m2{:});   
    m2 = squeeze(sum(m2,1,'omitnan'));
    z2 = zscore(rot90(m2),[],2);

    u = ismember(clumaa.uci,clumaa.uci(pidx & clumaa.partn==2 & clumaa.repetition_score(:,1)>rs_cutoff)); % hills data
    m1 = clumaa.(var)(u & clumaa.partn==1); % arena 1 data
    m1 = cellfun(@(x) imresize(x,[48 96],'bilinear'),m1,'UniformOutput',0);
    m1 = cat(3,m1{:});   
    m1 = squeeze(sum(m1,1,'omitnan'));
    z1 = zscore(rot90(m1),[],2);

    % ratemap peak vals
    [~,midx1] = max(z1,[],2,'omitnan');
    [~,idx1] = sort(midx1,'ascend');
    z1 = z1(idx1,:);

    [~,midx2] = max(z2,[],2,'omitnan');
    [~,idx2] = sort(midx2,'ascend');
    z2 = z2(idx2,:);

    % plotting settings
    xnow = 50;
    ynow = ynow-460;
    xbuff = 10;
    ysiz = 250;
    xsiz = 120;

    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]); % x-axis plot 
        ah = add_panel_title('b',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',[18 13]);                
    
        imagesc(z1,'XData',0:size(z1,2).*mapset.binsize./1000); hold on;
        axis ij on
        colormap(ax1,turbo)   
        ax1.CLim = [-1 4];
        ylabel('Ranked cells')
        ax1.YTick = unique([1 0:50:size(z1,1) size(z1,1)]);
        ax1.YLabel.FontSize = 12;
        ax1.XTickLabel = [];
    
        % additional plots
        maze_sections = 0:0.5:3;
        hill_x = maze_sections(2:2:end);
        plot([hill_x; hill_x],ax1.YLim,'w:','LineWidth',1)
    
        % text
        text(0.5,1,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','bottom')

    ax1i = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-70 ax1.Position(3) 60]); % x-axis plot 
        m = mean(z1,1,'omitnan');
        s = nansem(z1,1);
        xi = (1:size(z1,2)).*mapset.binsize./1000;

        bb = boundedline(xi,m,s,'cmap',plot_set{1,1});
        line(ax1i.XLim,[0 0],'Color',[.5 .5 .5])

        ax1i.XLim = ax1.XLim;
        ax1i.YLim = [-0.4 0.4];
        ax1i.XTick = ax1.XTick;
        xlabel('Position (m)')
        ax1i.YTick = -0.4:0.2:0.4;
        ylabel(sprintf('Mean %c SEM (z)',177))

    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]); % x-axis plot    
        imagesc(z2,'XData',0:size(z2,2).*mapset.binsize./1000); hold on;
        axis ij on
        colormap(ax2,turbo)   
        ax2.CLim = [-1 4];
        ax2.YTick = [];
        ax2.YLabel.FontSize = 12;
        ax2.XTickLabel = [];
    
        % additional plots
        maze_sections = 0:0.5:3;
        hill_x = maze_sections(2:2:end);
        plot([hill_x; hill_x],ax2.YLim,'w:','LineWidth',1)

        % text
        text(0.5,1,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','bottom')

    ax2i = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)-70 ax2.Position(3) 60]); % x-axis plot 
        m = mean(z2,1,'omitnan');
        s = nansem(z2,1);
        bb = boundedline(xi,m,s,'cmap',plot_set{1,2}); hold on;
        line(ax2i.XLim,[0 0],'Color',[.5 .5 .5])

        ax2i.XLim = ax2.XLim;
        ax2i.YLim = ax1i.YLim;        
        ax2i.XTick = ax2.XTick;
        xlabel('Position (m)')
        ax2i.YTick = [];
        ax2i.YColor = 'none';

    % colorbar
    axc = axes('Units','pixels','Position',[ax1.Position(1)+70 ax1.Position(2)+ysiz+20 100 10]);
        mat = (linspace(0,100,100));
        imagesc(ones(size(mat)),mat,mat);
        colormap(axc,ax1.Colormap);
        axis xy

        axc.YTick = [];
        axc.XTick = [];
        text(0.5,1.7,sprintf('Firing rate (z)'),'FontSize',7,'HorizontalAl','center','Units','normalized')
        text(1.02,0.5,sprintf('%.1f',ax1.CLim(2)),'FontSize',7,'HorizontalAl','left','Units','normalized')
        text(-0.02,0.5,sprintf('%.1f',ax1.CLim(1)),'FontSize',7,'HorizontalAl','right','Units','normalized')

%% >>>>>>>>>> Field distribution along hills
    ynow = ynow+120;
    xnow = xnow+310;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 200 120]);
        ah = add_panel_title('c',sprintf(''),'yoffset',0,'xoffset',-20,'width',400,'fontsize',[18 13]);                

        v1 = midx1 .* mapset.binsize ./ 1000; % convert bins to m
        v1 = wrapToPi(v1.*(2*pi)); % convert to phase and wrap to circle
        v2 = midx2 .* mapset.binsize ./ 1000; % convert bins to m
        v2 = wrapToPi(v2.*(2*pi)); % convert to phase and wrap to circle

        % get circular density of data
        ri = linspace(-pi,pi,360)';
        xi = movmean(ri,2,'Endpoints','discard');
        k = 0.25;
        f1 = circ_ksdensity(v1,xi,[],k);
        f2 = circ_ksdensity(v2,xi,[],k);
        f3 = circ_ksdensity(xi,xi,[],k);        
        ri = ri(:)';
        f1 = f1(:)';
        f2 = f2(:)';
        f3 = f3(:)';

        % main plot
        % p1 = polarhistogram('BinEdges',ri,'BinCounts',f1,'DisplayStyle','bar','EdgeColor','none','FaceColor',plot_set{1,1},'FaceAlpha',alph); hold on;
        % p2 = polarhistogram('BinEdges',ri,'BinCounts',f2,'DisplayStyle','bar','EdgeColor','none','FaceColor',plot_set{1,2},'FaceAlpha',alph);
        xip = rad2deg(xi) ./ 360;
        alph = 0.7;  
        a1 = area(xip,f1,'FaceColor',plot_set{1,1},'EdgeColor','none','FaceAlpha',alph); hold on;
        a2 = area(xip,f2,'FaceColor',plot_set{1,2},'EdgeColor','none','FaceAlpha',alph);
        a3 = plot(xip,(cos(xi)+1).*(ax.YLim(2)./2),'k','LineWidth',1.5);

        % axis settings
        ax.XLim = [-0.5 0.5]; 
        ax.FontSize = 10;
        box off
        xlabel(sprintf('Position relative to hilltop (m)'))    
        ylabel(sprintf('PDF'))    
        ytickformat('%.1f');

        % additional plots
        ch1 = line(xip,f3,'Color',[.5 .5 .5]);
        [~,leg] = legendflex([a1 a2 a3 ch1],{sprintf('%s',maze_names{1}),sprintf('%s',maze_names{2}),'Terrain','Chance'},'anchor',{'e','e'},'ncol',1,'box','off','buffer',[50,45],'xscale',.5,'fontsize',9); 
        leg(5).Children.FaceAlpha = a1.FaceAlpha;
        leg(6).Children.FaceAlpha = a2.FaceAlpha;

        % stats
        rv1 = circ_r(v1);
        rv2 = circ_r(v2);
        [ckp,ckf] = circ_ktest(v1(:),v2(:));
        text(0.03,1.15,sprintf('Rv = %.2f',rv1),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color',plot_set{1,1},'VerticalAlignment','bottom','rotation',0)
        text(0.03,1.05,sprintf('Rv = %.2f',rv2),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color',plot_set{1,2},'VerticalAlignment','bottom','rotation',0)
        text(0.42,1.05,sprintf('Arena vs Hills\n{\\itp} = %.2f',ckp),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S5.png']; 
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


 


















