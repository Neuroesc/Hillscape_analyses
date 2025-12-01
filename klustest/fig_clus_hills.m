function figCLUS(mtint,pdata,sdata,fig_vis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figCLUS  a basic figure function
% This function just stores the plotting for klustest. This figure gives all the information about a cluster
% in all parts (i.e. spikes, ratemap)
%
% USAGE:
%         figCLUS(mtint,pdata,sdata)
%
% INPUT:
%         mtint - mtint structure from klustest
%         pdata - pdata structure from klustest
%         sdata - sdatap structure from klustest
%
% See also: klustest figPART

% HISTORY:
% version 1.0.0, Release 12/04/19 Initial release, created to replace what was figPARTS
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
% deal with input variables
    inps = {'fig_vis','save_fig'};
    vals = {'''off''','0'};
    for ff = 1:length(inps)
        if ~exist(inps{ff},'var')
            eval([inps{ff} '=' vals{ff} ';']);
        end
    end

    if ~any(sdata.spikes) || all(isempty(sdata.spikes)) || all(isnan(sdata.spikes))
        return
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    fig_clust = figure('visible',fig_vis,'Units','pixels','Position',[50, 50, 1650, 920]); % open/create figure
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    fsiz = 10; % the fontsize for different texts
    flnw = 0.5; % the line width for different plots

    % add an annotation to the figure with some important info
    ann_str = sprintf('Cell: %s, Rat: %s, Date: %d, Tetrode: %d, Cluster: %d, Analysed: %s',sdata.uci{1},sdata.rat{1},sdata.date(1),sdata.tet(1),sdata.clu(1),datestr(now,'yyyy-mm-dd-HH-MM-SS'));
    annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',8,'LineStyle','none','interpreter','none');  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Position and spike plot (black lines & red dots)
    xnow = 50;
    ynow = 680;
    xwide = 350;
    ywide = 210;
    hbuff = xwide+40;
    vbuff = ywide; 
    axpos = [xnow,ynow,xwide,ywide; xnow+hbuff,ynow,xwide,ywide; xnow+hbuff*2,ynow,xwide,ywide; xnow+hbuff*3,ynow,xwide,ywide];
    prts = {1,1,'planar'; 2,1,'planar'; 2,2,'surficial'; 3,1,'planar'};
    for ii = 1:4
        % spikes and positions
        ax_ps = axes('Units','pixels','Position',axpos(ii,:));        
            % position and spike data (already cut to this part)
            part_now = prts{ii,1};            
            part_now_name = sdata.part{ prts{ii,1} };
            part_duration = sdata.duration( prts{ii,1} );  
            ppox = pdata.(part_now_name).(['pos_' prts{ii,3}])(:,1);
            ppoy = pdata.(part_now_name).(['pos_' prts{ii,3}])(:,2);
            ppot = pdata.(part_now_name).pot;
            pspx = ppox(sdata.spike_index{part_now});
            pspy = ppoy(sdata.spike_index{part_now});

            % plot position data, excluding pieces not included in this part
            % by inserting NaNs between intervals we can plot this as one line, which saves on memory
            dindax = abs([0; diff(ppot)])>0.1;
            pos_plot = [ppox ppoy];
            pos_plot(dindax,:) = NaN;
            plot(pos_plot(:,1),pos_plot(:,2),'Color',[.5 .5 .5]); hold on;

            % plot spikes after position data so they are all on top
            plot(pspx,pspy,'ro','MarkerFaceColor','r','MarkerSize',3)    

            % additional settings
            daspect([1 1 1])
            axis xy off tight
            text(.5,1.05,sprintf('%d spikes (%.2f Hz)',numel(pspx),numel(pspx)/part_duration),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')      
        
        % firing rate maps
        ax_ps = axes('Units','pixels','Position',axpos(ii,:)-[0 vbuff 0 0]);
            ratemap = sdata.( ['ratemap_' prts{ii,3}] ){ prts{ii,1} };
            spati = sdata.spatial_info_bsec(prts{ii,1},prts{ii,2});
            spars = sdata.sparsity(prts{ii,1},prts{ii,2});
            coher = sdata.spatial_coherence(prts{ii,1},prts{ii,2});

            im = imagesc(ratemap);
            set(im,'alphadata',~isnan(ratemap));
            daspect([1 1 1])
            colormap(gca,jet)  
            if ~any(ratemap(:)>0)
                caxis([0 1])                
            else
                caxis([0 nanmax(ratemap(:))])
            end     
            axis xy off tight
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz);   

            text(.5,1.05,sprintf('SI: %.2f, Sp: %.2f, Cohe: %.2f',spati,spars,coher),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
            tt = text(-0.05,0.5,sprintf('Max: %.2f, Mean: %.2f',nanmax(ratemap(:)),nanmean(ratemap(:))),'FontSize',fsiz,'Units','normalized','HorizontalAlignment','center');
            set(tt,'rotation',90);

            if ii==4
                axp = get(gca,'Position');
                cc = colorbar;
                set(gca,'Position',axp);
                set(cc,'Position',get(cc,'Position')+[-0.004 0 -0.004 0])
                title(cc,'Hz','FontSize',fsiz)
            end

        % grid autocorrelation
        ax_ps = axes('Units','pixels','Position',axpos(ii,:)-[0 vbuff*2 0 0]);
            amap = sdata.( ['amap_' prts{ii,3}] ){ prts{ii,1} };   
            msk = sdata.( ['amap_' prts{ii,3} '_mask'] ){ prts{ii,1} };
            if isnan(msk)
                msk = zeros(size(amap));
            end
            msk = ones(size(amap));

            imc = imagesc(amap);
            set(imc,'alphadata',double(msk)+0.3);    
            daspect([1 1 1])
            colormap(gca,inferno)       
            caxis([-0.2 1.0])
            axis xy off tight
            set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz); 
            text(.5,1.05,sprintf('G: %.2f, W: %.2f, O: %.2f',sdata.grid_score(prts{ii,1},prts{ii,2}),sdata.grid_wavelength(prts{ii,1},prts{ii,2}),sdata.grid_orientation(prts{ii,1},prts{ii,2})),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
            
            if ii==4
                axp = get(gca,'Position');
                cc = colorbar;
                set(gca,'Position',axp);
                set(cc,'Position',get(cc,'Position')+[-0.001 0 -0.004 0])
                title(cc,'Corr (r)','FontSize',fsiz)    
            end
            
        % 3D head direction
        ax_ps = axes('Units','pixels','Position',axpos(ii,:)-[0 vbuff*3-70 20 100]);            
            ratemap = sdata.hd_3d_ratemap{ prts{ii,1} };
            dwellmap = pdata.( part_now_name ).hd_3d_dwellmap;

            a = linspace(-180, 180, 65);
            e = linspace(-90, 90, 65);
            imagesc(a(:),e(:),ratemap,'alphadata',dwellmap>0.01);
            axis on tight xy
            view(0,90);    
            grid on
            colormap(gca,jet);
            if ~any(ratemap(:)>0)
                caxis([0 1])                
            else
                caxis([0 nanmax(ratemap(:))])
            end
            xlabel(sprintf('Yaw (%c)',176),'FontSize',fsiz)
            ylabel(sprintf('Pitch (%c)',176),'FontSize',fsiz)
            ax_ps.XTick = -180:90:180;
            ax_ps.YTick = -180:45:180;
            ax_ps.FontSize = fsiz;

            hold on;
            line(ax_ps.XLim,[0 0],'Color',[.5 .5 .5]);            
            
            
    end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Spike rasterplot running along bottom of figure
    axsr = axes('Units','pixels','Position',[30,35,1600,45]);
        % get the spike time data and bin it
        bsize = 5;
        edg = 0:bsize:pdata.duration; % vector of 1s time points at which we should calculate spike probability
        spt = mtint.tetrode(sdata.tet(1)).ts( mtint.tetrode(sdata.tet(1)).cut==sdata.clu(1) );
        f = histcounts(spt,edg);

        % plot this as an image
        im = imagesc(f);
        colormap(axsr,flipud(gray(256)))
        ax = gca;
        ax.YTick = [];
        ax.TickLength = [0.001, 0.001];
        ylabel('Spikes')
        hold on
        axis xy on
        xlabel('Time (s)')
        
        cols = jet(4);
        for ii = 1:4
            tvals = pdata.( sdata.part{ prts{ii,1} } ).times;
            for tt = 1:length(tvals(:,1))
                line([knnsearch(edg(:),tvals(tt,1)) knnsearch(edg(:),tvals(tt,2))],[0.5 0.5],'LineWidth',4,'Color',cols(ii,:));
            end
        end
        set(gca,'LineWidth',flnw,'layer','top','FontSize',fsiz,'XColor','k','YColor','k'); % set the axes to be on top, set their thickness to flnw, set all fonts to size fsiz and make sure the axis and text are black 
        set(ax,'xticklabel',num2str(get(gca,'xtick')'.*bsize,'%.f'))
    
        set(gcf,'visible','on')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    figfile = [pwd '\klustest\' pdata.combined_name '\figures\'];
    figname = [figfile sdata.uci{1} '_summary.png'];
    [~,~,~] = mkdir(figfile);
    %export_fig(figname,'-png','-r200') % t,r,b,l 
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',150);
    close(gcf);  
















































