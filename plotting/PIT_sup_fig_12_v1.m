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
%% >>>>>>>>>> Check that correlations and shuffles are complete
    clumaa = PIT_correlations(config,pidx,clumaa,posdata,0);

%% >>>>>>>>>> N fields in arena vs rep score in hills
    xnow = 50;
    ynow = 600;

    fname = [config.data_out_dir 'PIT_correlations_within_shuffles.mat'];
    load(fname,'shuffs_within');

    % create axis
    ax2 = axes('Units','pixels','Position',[xnow ynow 220 110]);  
        ah = add_panel_title('c',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400);         
        
        var = 'repetition_score';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
        unow = unique( clumaa.uci( clumaa.repetition_score(:,1) > prctile(v1,99) ) );
        pidx_now = ismember(clumaa.uci,unow) & pidx;
    
        pp = 2;
        var = 'within_session_stability';
        v1 = clumaa.(var)(pidx_now & clumaa.partn==pp,1); % within session stability, hills
        var = 'reflected_stability';
        v2 = clumaa.(var)(pidx_now & clumaa.partn==1,2); % reflected stability, hills
        % v3 = squeeze(shuffs_within(:,1,2)); % hills within session shuffles
        var = 'reflected_stability';
        v3 = clumaa.(var)(pidx & clumaa.partn==1,1); % reflected stability, arena 1

        xi = -1:0.01:1;
        bw = 0.05;
        f1 = ksdensity(v1(:),xi(:),"Bandwidth",bw);
        f2 = ksdensity(v2(:),xi(:),"Bandwidth",bw);
        f3 = ksdensity(v3(:),xi(:),"Bandwidth",bw);

        % main plot    
        alph = 0.7;          
        a1 = area(xi,f1,'FaceColor',plot_set{1,2},'EdgeColor','k','FaceAlpha',alph); hold on;
        a2 = area(xi,f2,'FaceColor','m','EdgeColor','k','FaceAlpha',alph);       
        a3 = area(xi,f3,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--'); hold on;       

        % axis settings
        ax2.XLim = [-0.5 1]; 
        ax2.FontSize = 8;
        ytickformat('%.1f')
        box off
        xlabel(sprintf('Correlation ({\\itr})'))    
        ylabel(sprintf('PDF'))  

        [~,leg] = legendflex([a1 a2 a3],...
            {sprintf('Session half 1 vs 2'),...
            sprintf('Hills reflected similarity'),sprintf('Arena 1 reflected similarity')},...
            'anchor',{'sw','sw'},'ncol',1,'box','off','buffer',[-10,-90],'xscale',.5,'fontsize',9); 
        
        leg(4).Children.FaceAlpha = a1.FaceAlpha;
        leg(5).Children.FaceAlpha = a2.FaceAlpha;
        leg(6).Children.FaceAlpha = a3.FaceAlpha;

%% >>>>>>>>>> Between session correlation results (Cohen's d values)
    % create axis
    ax = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+50 ynow 50 100]);
        % ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',-20,'width',400);
        e1 = stats_effect_size(v1,v2); % within vs reflect
        g1 = e1.cliff_delta;
        e2 = stats_effect_size(v2,v3); % reflect vs shuffle
        g2 = e2.cliff_delta;

        stem(1:2,[g1 g2],'filled','ko'); hold on;

        % axis settings
        ax.YLim = [-0.2 1];
        ax.XColor = 'none';
        ax.XLim = [0.5 2.5];
        ylabel(sprintf('Cliff''s Delta'))   
        ax.XTick = 1:2;
        ax.XTickLabel = {};
        ytickformat('%.1f');
        ax.FontSize = 8;
        box off

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]);

        % stats
        vcoef = 1;
        fsiz = 12;
        dat = {v1 v2;v2 v3};
        for ii = 1:2
            [z,p,~,~,~] = stats_perm_test(dat{ii,1},dat{ii,2},'method','cliff');
            if p<.001
                text(ii,ax.YLim(2).*vcoef,'***','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<0.01
                text(ii,ax.YLim(2).*vcoef,'**','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<.05
                text(ii,ax.YLim(2).*vcoef,'*','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            end
        end

        axi = axes('Units','pixels','Position',[ax.Position(1) ax.Position(2)-ax.Position(3)-5 ax.Position(3) ax.Position(3)]);
            axi.Color = 'none';
            axi.XLim = ax.XLim;
            axi.YLim = ax.XLim;
            rsiz = 0.6;
            rectangle('Position',[1-rsiz/2 1.75 rsiz rsiz],'FaceColor',[a1.FaceColor a1.FaceAlpha],'EdgeColor',a1.EdgeColor,'LineStyle',a1.LineStyle);
            rectangle('Position',[1-rsiz/2 0.6 rsiz rsiz],'FaceColor',[a2.FaceColor a2.FaceAlpha],'EdgeColor',a2.EdgeColor,'LineStyle',a2.LineStyle);
            text(1,1.2,'vs','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize)

            rectangle('Position',[2-rsiz/2 1.75 rsiz rsiz],'FaceColor',[a2.FaceColor a2.FaceAlpha],'EdgeColor',a2.EdgeColor,'LineStyle',a2.LineStyle);
            rectangle('Position',[2-rsiz/2 0.6 rsiz rsiz],'FaceColor',[a3.FaceColor a3.FaceAlpha],'EdgeColor',a3.EdgeColor,'LineStyle',a3.LineStyle);
            text(2,1.2,'vs','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize)

            axis off















    return
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig 6.png']; 
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







