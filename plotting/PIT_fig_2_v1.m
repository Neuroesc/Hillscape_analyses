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


% correlation cliffs delta
% field size and distance from boundaries
% percentage environment covered by fields
% range of field sizes (CV?) shuffle?
% nfields vs field size scatter
% nfields vs sum area of fields scatter

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
%% >>>>>>>>>> Check that correlations and shuffles are complete
    clumaa = PIT_correlations(config,pidx,clumaa,posdata,0);

%% >>>>>>>>>> Example cells (between session stability)
    xnow = 30;
    ynow = 780;
    xvec = [xnow xnow];
    yvec = [ynow ynow-60];
    ucis = {'RG19_221014_t3_c10','RG26_230119_t1_c9'};
    hsiz = 90;
    vsiz = 70;
    xbuff = 10;

    for uu = 1:length(ucis)
        m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==1)};        
        m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==2)};
        m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==3)};
    
        % plot arena 1
        ax1 = axes('Units','pixels','Position',[xvec(uu) yvec(uu) hsiz vsiz]);    
            imagesc(m1,'alphadata',~isnan(m1)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax1,turbo)   
            ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') max(m3(:),[],'omitnan') 1])];
            ax1.XTick = [];
            ax1.YTick = [];

            % text
            if uu==1
                text(0,1.05,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
                ah = add_panel_title('A',sprintf(''),'yoffset',-10,'xoffset',10,'width',400,'fontsize',fs);                                
            end
            text(-0.02,0.5,sprintf('Cell %d',uu),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')

        % plot hills
        ax2 = axes('Units','pixels','Position',ax1.Position+[ax1.Position(3)+xbuff 0 0 0]);
            imagesc(m2,'alphadata',~isnan(m2)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax2,turbo)   
            ax2.CLim = ax1.CLim;
            ax2.XTick = [];
            ax2.YTick = [];
    
            % additional plots
            ps = linspace(0,size(m2,2),7);
            tops = ps(2:2:end);
            valleys = ps(1:2:end);
            plot([tops; tops],ax2.YLim,'w:')

            % text
            if uu==1
                text(0,1.05,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
            end

        % plot arena 1
        ax3 = axes('Units','pixels','Position',ax2.Position+[ax2.Position(3)+xbuff 0 0 0]);    
            imagesc(m3,'alphadata',~isnan(m3)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax3,turbo)   
            ax3.CLim = ax1.CLim;
            ax3.XTick = [];
            ax3.YTick = [];

            % text
            if uu==1
                text(0,1.05,sprintf('%s',maze_names{3}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
            end
    end

    % vertical colormap
    axc = axes('Units','pixels','Position',[xnow+305 ynow-30 8 70]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100)';
        imagesc(x,'YData',ax1.CLim,'XData',[0 1]);
        colormap(ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0.5,1,sprintf('Max'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,0,sprintf('0 Hz'),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> Between session correlation results
    xnow = xnow+30;
    ynow = ynow-210;

    fname = [config.data_out_dir 'PIT_correlations_between_shuffles.mat'];
    load(fname,'shuffs');    

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 220 100]);
        % ah = add_panel_title('e',sprintf(''),'yoffset',0,'xoffset',-20,'width',400);
    
        var = 'between_session_stability';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,2); % arena 1 vs arena 2
        v2 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 vs hills
        v4 = shuffs(:,2); % arena 1 vs arena 2 shuffles
        v5 = shuffs(:,1); % arena 1 vs hills shuffles

        xi = -1:0.01:1;
        bw = 0.05;
        f1 = ksdensity(v1(:),xi(:),"Bandwidth",bw);
        f2 = ksdensity(v2(:),xi(:),"Bandwidth",bw);
        f4 = ksdensity(v4(:),xi(:),"Bandwidth",bw);
        f5 = ksdensity(v5(:),xi(:),"Bandwidth",bw);

        % main plot
        alph = 0.7;  
        a1 = area(xi,f1,'FaceColor','#012030','EdgeColor','k','FaceAlpha',alph); hold on;
        a2 = area(xi,f2,'FaceColor','#45C4B0','EdgeColor','k','FaceAlpha',alph);
        a4 = area(xi,f4,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--'); hold on;       
        a5 = area(xi,f5,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle',':'); hold on;       

        % axis settings
        ax.XLim = [-0.5 1]; 
        ax.FontSize = 10;
        set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))
        box off
        xlabel(sprintf('Correlation ({\\itr})'))    
        ylabel(sprintf('PDF'))    

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        [~,leg] = legendflex([a1 a2 a4 a5],...
            {sprintf('%s vs %s',maze_names{1},maze_names{3}),...
            sprintf('%s vs %s',maze_names{1},maze_names{2}),...
            sprintf('%s vs %s shuffle',maze_names{1},maze_names{3}),...
            sprintf('%s vs %s shuffle',maze_names{1},maze_names{2})},...
            'anchor',{'e','e'},'ncol',2,'box','off','buffer',[60,87],'xscale',.5,'fontsize',9); 
        
        leg(5).Children.FaceColor = a1.FaceColor;
        leg(6).Children.FaceColor = a2.FaceColor;
        leg(7).Children.FaceColor = a4.FaceColor;
        leg(8).Children.FaceColor = a5.FaceColor;        
        leg(5).Children.FaceAlpha = a1.FaceAlpha;
        leg(6).Children.FaceAlpha = a2.FaceAlpha;
        leg(7).Children.FaceAlpha = a4.FaceAlpha;
        leg(8).Children.FaceAlpha = a5.FaceAlpha;

        % % stats
        % [z,p,~,obs,~] = permutationZTest(v1,v4,'iti',1000,'method','cliff','rep',1); % arena 1&2 vs shuffle
        % disp(sprintf('arena 1&2 vs shuffle: z = %.3f, p = %.1e, Cliff''s delta = %.2f, permutation z-test',z,p,obs))
        % % median(v1,"all",'omitmissing')
        % % keyboard
        % 
        % [z,p,~,obs,~] = permutationZTest(v1,v2,'iti',1000,'method','cliff','rep',1); % arena 1&2 vs arena 1&hills
        % disp(sprintf('arena 1&2 vs arena 1&hills: z = %.3f, p = %.1e, Cliff''s delta = %.2f, permutation z-test',z,p,obs))
        % % keyboard
        % 
        % [z,p,~,obs,~] = permutationZTest(v2,v5,'iti',1000,'method','cliff','rep',1); % arena 1&hills vs shuffle
        % disp(sprintf('arena 1&hills vs shuffle: z = %.3f, p = %.1e, Cliff''s delta = %.2f, permutation z-test',z,p,obs))
        % % keyboard

%% >>>>>>>>>> Example cell (within session stability)
    xnow = xnow+310;
    ynow = ynow+150;
    
    uci = 'RG19_221021_t2_c1';
    m1a = clumaa.ratemap_planar_half{(ismember(clumaa.uci,uci) & clumaa.partn==1),1}; % arena half 1       
    m1b = clumaa.ratemap_planar_half{(ismember(clumaa.uci,uci) & clumaa.partn==1),2}; % arena half 2     
    m2a = clumaa.ratemap_planar_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),1}; % hills half 1     
    m2b = clumaa.ratemap_planar_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),2}; % hills half 2     

    % plot arena half 1
    ax1a = axes('Units','pixels','Position',[xnow ynow+60 90 70]);    
        imagesc(m1a,'alphadata',~isnan(m1a)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1a,turbo)   
        ax1a.CLim = [0 max([max(m1a(:),[],'omitnan') max(m1b(:),[],'omitnan') 1])];
        ax1a.XTick = [];
        ax1a.YTick = [];

        % text
        text(0,1.1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        ah = add_panel_title('B',sprintf(''),'yoffset',-10,'xoffset',10,'width',400,'fontsize',fs);                                
        text(-0.02,0.5,sprintf('1st half'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        text(0,0,sprintf('%.1f',ax1a.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')

    % plot arena half 2
    ax1b = axes('Units','pixels','Position',[xnow ynow 90 70]);    
        imagesc(m1b,'alphadata',~isnan(m1b)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1b,turbo)   
        ax1b.CLim = ax1a.CLim;
        ax1b.XTick = [];
        ax1b.YTick = [];

        % text
        text(-0.02,0.5,sprintf('2nd half'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

    % plot hills half 1
    ax2a = axes('Units','pixels','Position',[xnow+110 ynow+60 90 70]);
        imagesc(m2a,'alphadata',~isnan(m2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2a,turbo)   
        ax2a.CLim = [0 max([max(m2a(:),[],'omitnan') max(m2b(:),[],'omitnan') 1])];
        ax2a.XTick = [];
        ax2a.YTick = [];

        % additional plots
        ps = linspace(0,size(m2,2),7);
        tops = ps(2:2:end);
        valleys = ps(1:2:end);
        plot([tops; tops],ax2a.YLim,'w:')

        % text
        text(0,1.1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)        
        text(-0.02,0.5,sprintf('1st half'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        text(0,0,sprintf('%.1f',ax2a.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
        
    % plot hills half 2
    ax2b = axes('Units','pixels','Position',[xnow+110 ynow 90 70]);
        imagesc(m2b,'alphadata',~isnan(m2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2b,turbo)   
        ax2b.CLim = ax2a.CLim;
        ax2b.XTick = [];
        ax2b.YTick = [];

        % additional plots
        ps = linspace(0,size(m2,2),7);
        tops = ps(2:2:end);
        valleys = ps(1:2:end);
        plot([tops; tops],ax2b.YLim,'w:')

        % text
        text(-0.02,0.5,sprintf('2nd half'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

%% >>>>>>>>>> Within session correlation results
    ynow = ynow-150;
    % xnow = xnow;
    fname = [config.data_out_dir 'PIT_correlations_within_shuffles.mat'];
    load(fname,'shuffs_within');

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 220 100]);    
        var = 'within_session_stability';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1
        v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % hills
        v4 = squeeze(shuffs_within(:,1,1)); % arena shuffles
        v5 = squeeze(shuffs_within(:,1,2)); % hills shuffles

        xi = -1:0.01:1;
        bw = 0.05;
        f1 = ksdensity(v1(:),xi(:),"Bandwidth",bw);
        f2 = ksdensity(v2(:),xi(:),"Bandwidth",bw);
        f4 = ksdensity(v4(:),xi(:),"Bandwidth",bw);
        f5 = ksdensity(v5(:),xi(:),"Bandwidth",bw);

        % main plot    
        a1 = area(xi,f1,'FaceColor',plot_set{1,1},'EdgeColor','k','FaceAlpha',alph); hold on;
        a2 = area(xi,f2,'FaceColor',plot_set{1,2},'EdgeColor','k','FaceAlpha',alph);       
        a4 = area(xi,f4,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--'); hold on;       
        a5 = area(xi,f5,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle',':'); hold on;       

        % axis settings
        ax.XLim = [-0.5 1]; 
        ax.FontSize = 10;
        set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))
        box off
        xlabel(sprintf('Correlation ({\\itr})'))    
        ylabel(sprintf('PDF'))    

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        [~,leg] = legendflex([a1 a2 a4 a5],...
            {sprintf('%s',maze_names{1}),...
            sprintf('%s',maze_names{2}),...
            sprintf('%s shuffle',maze_names{1}),...
            sprintf('%s shuffle',maze_names{2})},...
            'anchor',{'e','e'},'ncol',2,'box','off','buffer',[-30,87],'xscale',.5,'fontsize',9); 

        leg(5).Children.FaceColor = a1.FaceColor;
        leg(6).Children.FaceColor = a2.FaceColor;
        leg(7).Children.FaceColor = a4.FaceColor;
        leg(8).Children.FaceColor = a5.FaceColor;
        leg(5).Children.FaceAlpha = a1.FaceAlpha;
        leg(6).Children.FaceAlpha = a2.FaceAlpha;
        leg(7).Children.FaceAlpha = a4.FaceAlpha;
        leg(8).Children.FaceAlpha = a5.FaceAlpha;

        % % % stats
        % [z,p,~,obs,~] = permutationZTest(v1,v4,'iti',1000,'method','cliff','rep',1); % arena 1 vs shuffles
        % disp(sprintf('arena 1 vs shuffles within: z = %.3f, p = %.1e, Cliff''s delta = %.2f, permutation z-test',z,p,obs))
        % % median(v1,"all",'omitmissing')
        % % keyboard
        % 
        % [z,p,~,obs,~] = permutationZTest(v2,v5,'iti',1000,'method','cliff','rep',1); % hills vs shuffle
        % disp(sprintf('hills vs shuffle within: z = %.3f, p = %.1e, Cliff''s delta = %.2f, permutation z-test',z,p,obs))
        % % median(v1,"all",'omitmissing')
        % % keyboard
        % 
        % [z,p,~,obs,~] = permutationZTest(v1,v2,'iti',1000,'method','cliff','rep',1); % arena 1 vs hills
        % disp(sprintf('arena 1 vs hills within: z = %.3f, p = %.1e, Cliff''s delta = %.2f, permutation z-test',z,p,obs))
        % % median(v1,"all",'omitmissing')
        % keyboard
% return
%% >>>>>>>>>> Firing rates
    xnow = 50;
    ynow = ynow-210;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 125 125]);
        ah = add_panel_title('C',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'f_rate_hz';
        ucis = unique(clumaa.uci(pidx));
        fst = NaN(length(ucis),3);
        for uu = 1:length(ucis)
            fst(uu,1) = clumaa.(var)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==1,1);% arena 1 data
            fst(uu,2) = clumaa.(var)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==2,1);% hills data
            fst(uu,3) = clumaa.(var)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==3,1);% arena 2 data
        end

        % main plot
        scatter(fst(:,1),fst(:,2),30,'k','filled','MarkerFaceColor','#45C4B0','MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;
        scatter(fst(:,1),fst(:,3),30,'k','filled','MarkerFaceColor','#012030','MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;

        % axis settings
        xlabel(sprintf('Firing rate 1'))            
        ylabel(sprintf('Firing rate 2'))    
        ax.XTick = [0.001 0.01 0.1 1 10 100];
        ax.YTick = [0.001 0.01 0.1 1 10 100];        
        ax.XLim = [0.001 5]; 
        ax.YLim = [0.001 5];         
        ax.FontSize = 8;
        ax.YScale = 'log';
        ax.XScale = 'log';

        % additional plots
        line([ax.XLim(1) 100],[ax.XLim(1) 100],'Color',[.5 .5 .5],'LineStyle','--')

        % stats
        [r,p] = corr(fst(:,1),fst(:,2),'rows','pairwise','type','Pearson');
        text(-0.02,1.12,sprintf('%s vs %s: r = %.1f',maze_names{1},maze_names{2},r),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','#45C4B0')
        [r,p] = corr(fst(:,1),fst(:,3),'rows','pairwise','type','Pearson');    
        text(-0.02,1.02,sprintf('%s vs %s: r = %.1f',maze_names{1},maze_names{3},r),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','#012030')

%% >>>>>>>>>> Fields per cell
    xnow = xnow+175;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 100 125]);
        ah = add_panel_title('D',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'planar_fields';
        pidx_temp = clumaa.cell_type(:,2)==3 & clumaa.planar_spatial_info_shuffles(:,2)>2;
        v1 = clumaa.(var)(pidx_temp & clumaa.partn==1,1); % arena 1 data
        v2a = clumaa.(var)(pidx_temp & clumaa.partn==2,1); % hills data
        v2b = clumaa.surficial_fields(pidx_temp & clumaa.partn==2,1); % hills data        
        v3 = clumaa.(var)(pidx_temp & clumaa.partn==3,1); % arena 2 data
        avg_fields = [median(v1(:),1,'omitnan') median(v2a(:),1,'omitnan') median(v2b(:),1,'omitnan') median(v3(:),1,'omitnan')];

        edg = -0.5:1:20;
        xi = movmean(edg,2,'omitnan','Endpoints','discard');
        [f1,e] = histcounts(v1,edg,'Normalization','probability');
        [f2a,e] = histcounts(v2a,edg,'Normalization','probability');
        [f2b,e] = histcounts(v2b,edg,'Normalization','probability');
        [f3,e] = histcounts(v3,edg,'Normalization','probability');
        f1 = [f1(1:8) sum(f1(9:end))];
        f2a = [f2a(1:8) sum(f2a(9:end))];
        f2b = [f2b(1:8) sum(f2b(9:end))];
        f3 = [f3(1:8) sum(f3(9:end))];
        xi = xi(1:9);

        p1 = plot(xi,f1,'Color',plot_set{1,1},'Marker',plot_set{2,1},'MarkerFaceColor',plot_set{1,1}); hold on;
        p2 = plot(xi,f2a,'Color',plot_set{1,2},'Marker',plot_set{2,2},'MarkerFaceColor',plot_set{1,2});
        % plot(xi,f3,'Color',plot_set{1,1},'Marker','o','MarkerFaceColor',plot_set(1,1));

        % axis settings
        ylabel(sprintf('Place cells (prop.)'))    
        ax.XTick = 0:8;
        ax.XLim = [0 8];
        ax.XTickLabel{end} = '>7';
        ax.XTickLabelRotation = 0;
        xlabel('Place fields')
        box off
        ax.FontSize = 8;
        ytickformat('%.1f');

        [~,leg] = legendflex([p1 p2],{maze_names{1},maze_names{2}},'anchor',{'nw','nw'},'ncol',1,'box','off','buffer',[25,10],'xscale',.5,'fontsize',9); 

        % stats
        [~,p,k] = kstest2(v1(:),v2a(:));
        text(-0.02,1.1,sprintf('D = %.1f, {\\itp} < .001',k),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom')
        % disp(sprintf('D = %.3f, p = %.1e, two-sample Kolmogorov-Smirnov test',k,p))
        % [median(v1(:),"all",'omitmissing') median(v2a(:),"all",'omitmissing')];

%% >>>>>>>>>> Fields per m2
    xnow = xnow+145;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 90 125]);
        ah = add_panel_title('E',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        % var = 'planar_fields';
        % v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
        % v2a = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data
        % v2b = clumaa.surficial_fields(pidx & clumaa.partn==2,1); % hills data        
        % v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 2 data

        var = 'planar_fields';
        v1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.repetition_score(:,1)<0.28,1); % arena 1 data
        v2a = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.repetition_score(:,1)<0.28,1); % hills data
        v2b = clumaa.surficial_fields(pidx & clumaa.partn==2 & clumaa.repetition_score(:,1)<0.28,1); % hills data        
        v3 = clumaa.(var)(pidx & clumaa.partn==3 & clumaa.repetition_score(:,1)<0.28,1); % arena 2 data

        % convert to fields per m2
        v1 = v1 ./ (3*1.5);
        v2a = v2a ./ (3*1.5);
        v2b = v2b ./ ((3*1.2)*1.5);
        v3 = v3 ./ (3*1.5);

        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dotcolor',plot_set(1,[1 2 2 3]),'dots',1,'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Place fields / m^{2}'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        % ax.XTickLabel = {maze_names{1} maze_names{2} sprintf([maze_names{2} ' (flattened)']) maze_names{3}};
        ax.FontSize = 8;
        % ax.XTickLabelRotation = 25;
        % ax.YScale = 'log';
        % ax.YLim(1) = 0;
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s          \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')
        % set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))
        ytickformat('%.1f');

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2));

% keyboard
%% >>>>>>>>>> Example cell (multiple fields)
    xnow = xnow+120;
    
    uci = 'RG11_220520_t2_c6';
    m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
    m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
    m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+65 110 70]);
        imagesc(m1,'alphadata',~isnan(m1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,turbo)   
        ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(-0.02,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        text(0,1.1,sprintf('Example cell'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')

    % plot hills
    ax2 = axes('Units','pixels','Position',[xnow ynow-15 110 70]);
        imagesc(m2,'alphadata',~isnan(m2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,turbo)   
        ax2.CLim = ax1.CLim;
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(m2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'w:')

        % text
        text(-0.02,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+10 ynow+55 70 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f Hz',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')

% %% >>>>>>>>>> N fields vs avg field area (scatter)
%     xnow = xnow+160;
% 
%     uci = 'RG11_220520_t2_c6';
%     m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
%     m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
%     m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};
% 
%     % create axis
%     ax = axes('Units','pixels','Position',[xnow ynow 125 125]);
%         ah = add_panel_title('c',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',[18 13]);                
% 
%         % collect data
%         ucis = unique(clumaa.uci(pidx));
%         dat = NaN(length(ucis),2,3);
%         for pp = 1:3    
%             for uu = 1:length(ucis) % for each cell
%                 idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==pp;
%                 if pp==2
%                     fdata = clumaa.surficial_field_data{idx}; % field data for this cell in this session in this part
%                     rmap = clumaa.ratemap_surficial{idx};       
%                 else
%                     fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
%                     rmap = clumaa.ratemap_planar{idx};                       
%                 end
%                 if isempty(fdata)
%                     continue
%                 end
% 
%                 dat(uu,:,pp) = [size(fdata,1) mean(fdata.Area(:),'all','omitmissing')./1e4];
%             end
%         end
% 
%         nindx = isnan(dat(:,1,1));
%         v1 = accumarray(dat(~nindx,1,1),dat(~nindx,2,1),[10,1],@mean,NaN);
%         v1s = accumarray(dat(~nindx,1,1),dat(~nindx,2,1),[10,1],@nansem,NaN);
%         nindx = isnan(dat(:,1,2));        
%         v2 = accumarray(dat(~nindx,1,2),dat(~nindx,2,2),[10,1],@mean,NaN);
%         v2s = accumarray(dat(~nindx,1,2),dat(~nindx,2,2),[10,1],@nansem,NaN);        
%         xi = 1:10;
% 
%         % main plot
%         errorbar(xi-0.1,v1,v1s,'Color',plot_set{1,1},'Marker','.'); hold on;
%         errorbar(xi+0.1,v2,v2s,'Color',plot_set{1,2},'Marker','.');         
%         % scatter(dat(:,1,1)-0.25,dat(:,2,1),20,'k','filled','MarkerFaceColor',plot_set{1,1},'MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;
%         % scatter(dat(:,1,2)+0.25,dat(:,2,2),20,'k','filled','MarkerFaceColor',plot_set{1,2},'MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;
%         % refline
% 
%         % axis settings
%         xlabel(sprintf('Total fields'))            
%         ylabel(sprintf('Avg. field area (m^{2})'))    
%         ax.XTick = 1:8;
%         % ax.YTick = [0.001 0.01 0.1 1 10 100];        
%         ax.XLim = [0 9]; 
%         ax.YLim(1) = 0;         
%         ax.FontSize = 8;
%         % ax.YScale = 'linear';
%         % ax.XScale = 'log';
%         ytickformat('%.1f')
%         box off
% 
%         % additional plots
%         % line([ax.XLim(1) 100],[ax.XLim(1) 100],'Color',[.5 .5 .5],'LineStyle','--')
% 
%         % % stats
%         % [r,p] = corr(fs(:,1),fs(:,2),'rows','pairwise','type','Pearson');
%         % text(0,1.12,sprintf('Arena 1 vs Hills: r = %.1f',r),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','#45C4B0')
%         % [r,p] = corr(fs(:,1),fs(:,3),'rows','pairwise','type','Pearson');    
%         % text(0,1.02,sprintf('Arena 1 vs Arena 2: r = %.1f',r),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','#012030')

%% >>>>>>>>>> field radius
    xnow = 50;
    ynow = ynow-215;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 100 125]);
        ah = add_panel_title('F',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'planar_field_info';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
        v2a = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data
        v2b = clumaa.surficial_field_info(pidx & clumaa.partn==2,1); % hills data        
        v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds.*10,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Field radius (mm)'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        ax.FontSize = 8;
        ax.YScale = 'log';
        ax.YTick = [100 200 400 800];
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s          \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);
            N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];

        % % spatial info stats
        % var = 'planar_spatial_info_shuffles';
        % v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
        % v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data
        % [z,p,~,obs,~] = stats_perm_test(v1,v2,'iti',1000,'method','cliff','rep',1); % arena 1&hills vs shuffle
        % disp(sprintf('z = %.3f, p = %.1e, Cliff''s delta = %.2f, permutation z-test',z,p,obs))
        % [median(v1(:),"all",'omitmissing') median(v2(:),"all",'omitmissing')];

%% >>>>>>>>>> field radius schematic
    xnow = xnow+135;
    ynow = ynow-15;
    mf = posdata.maze_frame{1}; 

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+80 110 70]);
        ah = add_panel_title('G',sprintf(''),'yoffset',-10,'xoffset',0,'width',400,'fontsize',fs);                
    
        plot(mf(:,1),mf(:,2),'Color','k'); hold on;

        x = linspace(min(mf(:,1)),max(mf(:,1)),7);
        xv = x([3 5]);
        y = linspace(min(mf(:,2)),max(mf(:,2)),5);
        yv = y([4 2]);
        for ii = 1:length(xv)
            rad = [mean(v1(:),1,'omitnan') mean(v1(:),1,'omitnan')].*10;
            cent = [xv(ii)-rad(1) yv(ii)-rad(1)];
            rectangle('Position',[cent rad.*2],'Curvature',1,'FaceColor',plot_set{1,1})
        end

        % axis settings
        daspect([1 1 1])
        axis off xy tight
        ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax1.YLim = [min(mf(:,2)) max(mf(:,2))];

        % text
        text(-0.02,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        text(0,1.1,sprintf('Exemplar firing activity'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

    ax1i = axes('Units','pixels','Position',[ax1.Position(1)+10 ax1.Position(2)+10 12 12]);
        img = imread([config.main_dir '\associated media\rat_silh_v2.png']);
        imshow(img);

        % axis settings
        daspect([1 1 1])
        axis off

    % plot hills
    ax2 = axes('Units','pixels','Position',[xnow ynow 110 70]);
        plot(mf(:,1),mf(:,2),'Color','k'); hold on;

        x = linspace(min(mf(:,1)),max(mf(:,1)),9);
        xv = x([3 5 7]);
        y = linspace(min(mf(:,2)),max(mf(:,2)),5);
        yv = y([4 3 2]);
        for ii = 1:length(xv)
            rad = [mean(v2a(:),1,'omitnan') mean(v2a(:),1,'omitnan')].*10;
            cent = [xv(ii)-rad(1) yv(ii)-rad(1)];
            rectangle('Position',[cent rad.*2],'Curvature',1,'FaceColor',plot_set{1,2})
        end

        % axis settings
        daspect([1 1 1])
        axis off xy tight
        ax2.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax2.YLim = [min(mf(:,2)) max(mf(:,2))];

        % text
        text(-0.02,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

    ax1i = axes('Units','pixels','Position',[ax2.Position(1)+10 ax2.Position(2)+10 12 12]);
        imshow(img);

        % axis settings
        daspect([1 1 1])
        axis off

% %% >>>>>>>>>> Example cell (field size)
%     xnow = xnow+140;
% 
%     uci = 'RG19_221013_t1_c10';
%     m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
%     m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
%     m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};
% 
%     % plot arena
%     ax1 = axes('Units','pixels','Position',[xnow ynow+80 110 70]);    
%         imagesc(m1,'alphadata',~isnan(m1)); hold on;
%         axis xy on
%         daspect([1 1 1])
%         colormap(ax1,turbo)   
%         ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
%         ax1.XTick = [];
%         ax1.YTick = [];
% 
%         % text
%         text(0.1,1.1,sprintf('Example cell'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
%         text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
%         text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
% 
%     % plot hills
%     ax2 = axes('Units','pixels','Position',[xnow ynow 110 70]);
%         imagesc(m2,'alphadata',~isnan(m2)); hold on;
%         axis xy on
%         daspect([1 1 1])
%         colormap(ax2,turbo)   
%         ax2.CLim = ax1.CLim;
%         ax2.XTick = [];
%         ax2.YTick = [];
% 
%         % additional plots
%         ps = linspace(0,size(m2,2),7);
%         tops = ps(2:2:end);
%         valleys = ps(1:2:end);
%         plot([tops; tops],ax2.YLim,'w')
% 
%         % text
%         text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
% 
%     % horizontal colormap
%     axc = axes('Units','pixels','Position',[xnow+10 ynow+70 70 8]);
%         x = linspace(ax1.CLim(1),ax1.CLim(2),100);
%         imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
%         colormap(ax1.Colormap);
%         axis on
%         axc.XTick = [];
%         axc.YTick = [];
%         axc.FontSize = 7;
%         text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
%         text(1,0.5,sprintf(' %.1f Hz',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> field coverage
    xnow = xnow+160;
    ynow = ynow+15;

    % field area analysis (if required)
    if 1%~any(ismember(clumaa.Properties.VariableNames,'planar_field_area')) % if the column(s) do not exist yet
        clumaa.planar_field_area = NaN(size(clumaa,1),2); % preallocate                
        for ii = 1:size(clumaa,1) 
            f = clumaa.planar_field_data{ii};
            m = clumaa.ratemap_planar{ii};
            ma = sum(~isnan(m(:))) .* ((32/10)^2);            
            if ~isempty(f)
                clumaa.planar_field_area(ii,:) = [sum(f.Area(:),'all','omitmissing') sum(f.Area(:),'all','omitmissing')./ma.*100];
            end

            f = clumaa.surficial_field_data{ii};
            m = clumaa.ratemap_surficial{ii};
            ma = sum(~isnan(m(:))) .* ((32/10)^2);
            if ~isempty(f)
                clumaa.surficial_field_area(ii,:) = [sum(f.Area(:),'all','omitmissing') sum(f.Area(:),'all','omitmissing')./ma.*100];    
            end
        end
    end

    % plot data
    ax1 = axes('Units','pixels','Position',[xnow ynow 100 125]);
        ah = add_panel_title('H',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);                

        var = 'planar_field_area';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,2); % arena 1 data
        v2a = clumaa.(var)(pidx & clumaa.partn==2,2); % hills data
        v2b = clumaa.surficial_field_area(pidx & clumaa.partn==2,2); % hills surficial data        
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotsize',20,'dotalpha',0.5,'dotsigma',0.1,'dotcolor',plot_set(1,[1 2 2 3]));

        % axis settings
        ylabel(sprintf('%% maze covered by fields'))    
        ax = gca;
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.YTick = [100 200 400 800];
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s          \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')
        ytickformat('%.f');

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'plot_omnibus',0);

%% >>>>>>>>>> spatial information content
    xnow = xnow+155;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 100 125]);
        ah = add_panel_title('I',sprintf(''),'yoffset',0,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'planar_spatial_info_shuffles';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,2); % arena 1 data
        v2a = clumaa.(var)(pidx & clumaa.partn==2,2); % hills data
        v2b = clumaa.surficial_spatial_info_shuffles(pidx & clumaa.partn==2,2); % arena 1 data        
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotsize',20,'dotalpha',0.5,'dotsigma',0.1,'dotcolor',plot_set(1,[1 2 2 3]));

        % axis settings
        ylabel(sprintf('Spatial information (z)'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.YTick = [100 200 400 800];
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s          \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(0.5,1.1,sprintf('n.s.'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize,'Units','normalized')

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'plot_omnibus',0);
            N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];

            % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig 2.png']; 
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

return












