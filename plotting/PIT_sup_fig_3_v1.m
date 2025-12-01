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
%% >>>>>>>>>> Check that correlations and shuffles are complete
    clumaa = PIT_correlations(config,pidx,clumaa,posdata,0);

%% >>>>>>>>>> Example cells
    xnow = 30;
    ynow = 780;
    xsiz = [120, 20]; % size, buffer
    ysiz = [70, 30]; % size, buffer
    xvec = xnow : sum(xsiz) : 1000;
    xvec = xvec(1:4);
    yvec = ynow : -sum(ysiz).*1.9 : -1000;
    yvec = yvec(1);
    [x,y] = ndgrid(xvec,yvec);
    ar_siz = [xsiz*0.48 ysiz*0.5];

    % find which cells to plot  
    idx = find(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5 & clumaa.planar_spatial_info_shuffles(:,2)>30); % place cells in the hills
    [~,sidx] = sort(clumaa.repetition_score(idx,1),'descend');
    idx = idx(sidx);
    idx = idx([35 35 35 35]);

    for ii = 1:numel(x)
        uci = clumaa.uci{idx(ii)};

        %% plot arena (ratemap)
        m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};
        if ii==4
            m2 = clumaa.ratemap_surficial{(ismember(clumaa.uci,uci) & clumaa.partn==2)};   
        else
            m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};   
        end
        m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};

        ax3 = axes('Units','pixels','Position',[x(ii) y(ii) xsiz(1) ysiz(1)]);
            if ii==1
                % correlations limited to edges
                msk = get_map_mask(m1a,'edges');
            elseif ii==2
                % correlations limited to corners
                msk = get_map_mask(m1a,'corners');
            elseif ii==3
                % correlations limited to overlapping bins
                msk = get_map_mask(m1a,'overlapping');
            elseif ii==4
                % correlations aligned
                msk = ones(size(m1));
            end
            msk = double(msk);        
            msk(msk==0) = 0.3;

            imagesc(m1,'alphadata',msk); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax3,turbo)   
            ax3.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
            ax3.XTick = [];
            ax3.YTick = [];
            if ii==1
                ah = add_panel_title('A',sprintf(''),'yoffset',-20,'xoffset',10,'width',400,'fontsize',fs);                
                text(-0.02,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
                text(0,1,sprintf('Edges'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)                        
            elseif ii==2
                text(0,1,sprintf('Corners'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)  
            elseif ii==3
                text(0,1,sprintf('Volumetric overlap'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)                                        
            elseif ii==4
                text(0,1,sprintf('Geodesic overlap'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)                        
            end

        %% plot hills (ratemap)
        if ii==4
            ax4 = axes('Units','pixels','Position',[x(ii) ax3.Position(2)-ax3.Position(4)-20 xsiz(1)+30 ysiz(1)]);
        else
            ax4 = axes('Units','pixels','Position',[x(ii) ax3.Position(2)-ax3.Position(4)-20 xsiz(1) ysiz(1)]);
        end
            if ii==1
                % correlations limited to edges
                msk = get_map_mask(m1a,'edges');
            elseif ii==2
                % correlations limited to corners
                msk = get_map_mask(m1a,'corners');
            elseif ii==3
                % correlations limited to overlapping bins
                msk = get_map_mask(m1a,'overlapping');
            elseif ii==4
                % correlations aligned
                msk = zeros(size(m2));
                msk(:,1:size(m1,2)) = 1;
            end
            msk = double(msk);
            msk(msk==0) = 0.3;

            imagesc(m2,'alphadata',msk); hold on;
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
                text(-0.02,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end
    end

%% arrows etc
    ax = axes('Units','pixels','Position',[xnow ynow-50 500 200],'Color','none');    
        ax.XLim = [0 1];
        ax.YLim = [0 1];
        ax.Color = 'none';
        axis manual off

        p1 = [0.12 0.265];
        p2 = [0.12 0.13];
        arrow(p1,p2,'Length',10,'TipAngle',20,'Width',1); hold off;

        p1 = [0.41 0.265];
        p2 = [0.41 0.13];
        arrow(p1,p2,'Length',10,'TipAngle',20,'Width',1); hold off;

        p1 = [0.68 0.265];
        p2 = [0.68 0.13];
        arrow(p1,p2,'Length',10,'TipAngle',20,'Width',1); hold off;

        p1 = [0.95 0.265];
        p2 = [0.95 0.13];
        arrow(p1,p2,'Length',10,'TipAngle',20,'Width',1); hold off;

%% >>>>>>>>>> Between session correlation results
    xnow = 60;
    ynow = ynow-310;
    fname = [config.data_out_dir 'PIT_correlations_between_shuffles.mat'];
    load(fname,'shuffs');

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 350 140]);
        ah = add_panel_title('B',sprintf(''),'yoffset',0,'xoffset',-20,'width',400,'fontsize',fs);
    
        var = 'between_session_stability';
        v1a = clumaa.(var)(pidx & clumaa.partn==1,5); % arena 1 vs hills (edges)
        v2a = clumaa.(var)(pidx & clumaa.partn==1,9); % arena 1 vs hills (corners)
        v3a = clumaa.(var)(pidx & clumaa.partn==1,3); % arena 1 vs hills (overlap)        
        v4a = max(clumaa.(var)(pidx & clumaa.partn==1,11:12),[],2,'omitnan'); % arena 1 vs hills (ends)
        % v1b = shuffs(:,5); % edges shuffled
        % v2b = shuffs(:,9); % corners shuffled
        % v3b = shuffs(:,3); % overlap shuffled        
        % v4b = max(shuffs(:,11:12),[],2,'omitnan'); % ends aligned shuffled
        v1b = clumaa.(var)(pidx & clumaa.partn==1,6); % arena 1 vs arena 2 (edges)
        v2b = clumaa.(var)(pidx & clumaa.partn==1,10); % arena 1 vs arena 2 (corners)
        v3b = clumaa.(var)(pidx & clumaa.partn==1,4); % arena 1 vs arena 2 (overlap)        
        v4b = clumaa.(var)(pidx & clumaa.partn==1,2); % arena 1 vs arena 2 (normal all map correlations)

        xi = -1:0.01:1;
        bw = 0.05;
        f1a = ksdensity(v1a(:),xi(:),"Bandwidth",bw);
        f2a = ksdensity(v2a(:),xi(:),"Bandwidth",bw);
        f3a = ksdensity(v3a(:),xi(:),"Bandwidth",bw);
        f4a = ksdensity(v4a(:),xi(:),"Bandwidth",bw);
        f1b = ksdensity(v1b(:),xi(:),"Bandwidth",bw);
        f2b = ksdensity(v2b(:),xi(:),"Bandwidth",bw);
        f3b = ksdensity(v3b(:),xi(:),"Bandwidth",bw);
        f4b = ksdensity(v4b(:),xi(:),"Bandwidth",bw);

        % main plot
        alph = 0.7;  
        cols = winter(4);
        a1a = area(xi,f1a,'FaceColor',cols(1,:),'EdgeColor','k','FaceAlpha',alph); hold on;
        a2a = area(xi,f2a,'FaceColor',cols(2,:),'EdgeColor','k','FaceAlpha',alph);
        a3a = area(xi,f3a,'FaceColor',cols(3,:),'EdgeColor','k','FaceAlpha',alph);
        a4a = area(xi,f4a,'FaceColor',cols(4,:),'EdgeColor','k','FaceAlpha',alph);   
        a1b = area(xi,f1b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--');
        a2b = area(xi,f2b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle',':');
        a3b = area(xi,f3b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','-');
        a4b = area(xi,f4b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','-.');      

        % axis settings
        ax.XLim = [-0.5 1]; 
        ax.FontSize = 10;
        ytickformat('%.1f');
        box off
        xlabel(sprintf('Correlation ({\\itr})'))    
        ylabel(sprintf('PDF'))    

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        [~,leg] = legendflex([a1a a2a a3a a4a a1b a2b a3b a4b],...
            {sprintf('Edges (%s\\times%s)',maze_names{4},maze_names{5}),...
            sprintf('Corners (%s\\times%s)',maze_names{4},maze_names{5}),...
            sprintf('Volumetric (%s\\times%s)',maze_names{4},maze_names{5}),...            
            sprintf('Geodesic (%s\\times%s)',maze_names{4},maze_names{5}),...
            sprintf('Edges (%s\\times%s)',maze_names{4},maze_names{6}),...
            sprintf('Corners (%s\\times%s)',maze_names{4},maze_names{6}),...
            sprintf('Volumetric (%s\\times%s)',maze_names{4},maze_names{6}),...            
            sprintf('Geodesic (%s\\times%s)',maze_names{4},maze_names{6})},...
            'anchor',{'e','e'},'ncol',2,'box','off','buffer',[-20,100],'xscale',.5,'fontsize',8); 
        leg(9).Children.FaceAlpha = a1a.FaceAlpha;
        leg(10).Children.FaceAlpha = a2a.FaceAlpha;
        leg(11).Children.FaceAlpha = a3a.FaceAlpha;
        leg(12).Children.FaceAlpha = a4a.FaceAlpha;
        leg(13).Children.FaceAlpha = a1b.FaceAlpha;
        leg(14).Children.FaceAlpha = a2b.FaceAlpha;
        leg(15).Children.FaceAlpha = a3b.FaceAlpha;
        leg(16).Children.FaceAlpha = a4b.FaceAlpha;

        leg(9).Children.FaceColor = a1a.FaceColor;
        leg(10).Children.FaceColor = a2a.FaceColor;
        leg(11).Children.FaceColor = a3a.FaceColor;
        leg(12).Children.FaceColor = a4a.FaceColor;
        leg(13).Children.FaceColor = a1b.FaceColor;
        leg(14).Children.FaceColor = a2b.FaceColor;
        leg(15).Children.FaceColor = a3b.FaceColor;
        leg(16).Children.FaceColor = a4b.FaceColor;

%% >>>>>>>>>> Between session correlation results (Cohen's d values)
    xnow = xnow+400;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 120 ax.Position(4)]);
        % ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',-20,'width',400);
        e1 = computeEffectSize(v1a,v1b);
        g1 = e1.cliff_delta;
        e2 = computeEffectSize(v2a,v2b);
        g2 = e2.cliff_delta;
        e3 = computeEffectSize(v3a,v3b);
        g3 = e3.cliff_delta;
        e4 = computeEffectSize(v4a,v4b);
        g4 = e4.cliff_delta;

        v = [g1 g2 g3 g4];
        for ii=1:4
            stem(ii,v(ii),'filled','MarkerFaceColor',cols(ii,:),'Color',cols(ii,:)); hold on;
        end

        % axis settings
        ax.YLim = [-1 1];
        ax.XLim = [0.5 4.5];
        ylabel(sprintf('Cliff''s Delta'))   
        ax.XTick = 1:4;
        ax.XTickLabel = {'Edges','Corners','Volumetric','Geodesic'};
        ytickformat('%.1f');
        ax.FontSize = 8;
        box off

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]);

        % stats
        vcoef = 1.1;
        fsiz = 12;
        dat = {v1a v1b;v2a v2b;v3a v3b;v4a v4b};
        for ii = 1:4
            [z,p,~,~,~] = permutationZTest(dat{ii,1},dat{ii,2},'method','cliff');
            if p<.001
                text(ii,ax.YLim(2).*vcoef,'***','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<0.01
                text(ii,ax.YLim(2).*vcoef,'**','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<.05
                text(ii,ax.YLim(2).*vcoef,'*','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            end
        end

%% >>>>>>>>>> Within session correlation results, Arena 1
    xnow = 60;
    ynow = ynow-220;
    fname = [config.data_out_dir 'PIT_correlations_within_shuffles.mat'];
    load(fname,'shuffs_within');

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 350 120]);
        ah = add_panel_title('C',sprintf(''),'yoffset',0,'xoffset',-20,'width',400,'fontsize',fs);

        var = 'within_session_stability';
        v1a = clumaa.(var)(pidx & clumaa.partn==1,4); % arena 1 (edges)
        v2a = clumaa.(var)(pidx & clumaa.partn==1,5); % arena 1 (center)   
        v3a = clumaa.(var)(pidx & clumaa.partn==1,6); % arena 1 (corners)        
        v1b = squeeze(shuffs_within(:,2,1)); % edges shuffled
        v2b = squeeze(shuffs_within(:,3,1)); % center shuffled        
        v3b = squeeze(shuffs_within(:,4,1)); % corners shuffled

        xi = -1:0.01:1;
        bw = 0.05;
        f1a = ksdensity(v1a(:),xi(:),"Bandwidth",bw);
        f2a = ksdensity(v2a(:),xi(:),"Bandwidth",bw);
        f3a = ksdensity(v3a(:),xi(:),"Bandwidth",bw);
        f1b = ksdensity(v1b(:),xi(:),"Bandwidth",bw);
        f2b = ksdensity(v2b(:),xi(:),"Bandwidth",bw);
        f3b = ksdensity(v3b(:),xi(:),"Bandwidth",bw);

        % main plot
        alph = 0.7;  
        cols = winter(3);
        a1a = area(xi,f1a,'FaceColor',cols(1,:),'EdgeColor','k','FaceAlpha',alph); hold on;
        a2a = area(xi,f2a,'FaceColor',cols(2,:),'EdgeColor','k','FaceAlpha',alph);
        a3a = area(xi,f3a,'FaceColor',cols(3,:),'EdgeColor','k','FaceAlpha',alph);
        a1b = area(xi,f1b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--');
        a2b = area(xi,f2b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle',':');
        a3b = area(xi,f3b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','-');

        % axis settings
        ax.XLim = [-0.5 1]; 
        ax.FontSize = 10;
        ytickformat('%.1f');
        box off
        ylabel(sprintf('PDF'))   
        ax.XTick = [];

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        [~,leg] = legendflex([a1a a2a a3a a1b a2b a3b],...
            {sprintf('Edges'),...
            sprintf('Centre'),...
            sprintf('Corners'),...            
            sprintf('Edges (shuffle)'),...
            sprintf('Centre (shuffle)'),...
            sprintf('Corners (shuffle)')},...            
            'anchor',{'e','e'},'ncol',2,'box','off','buffer',[0,50],'xscale',.5,'fontsize',9); 
        leg(7).Children.FaceAlpha = a1a.FaceAlpha;
        leg(8).Children.FaceAlpha = a2a.FaceAlpha;
        leg(9).Children.FaceAlpha = a3a.FaceAlpha;
        leg(10).Children.FaceAlpha = a4a.FaceAlpha;
        leg(11).Children.FaceAlpha = a1b.FaceAlpha;
        leg(12).Children.FaceAlpha = a2b.FaceAlpha;

        leg(7).Children.FaceColor = a1a.FaceColor;
        leg(8).Children.FaceColor = a2a.FaceColor;
        leg(9).Children.FaceColor = a3a.FaceColor;
        leg(10).Children.FaceColor = a4a.FaceColor;
        leg(11).Children.FaceColor = a1b.FaceColor;
        leg(12).Children.FaceColor = a2b.FaceColor;

        text(ax,0,1,maze_names{1},'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',10,'Color','k');

%% >>>>>>>>>> Within session correlation results (Cohen's d values) Arena 1
    xnow = xnow+400;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 120 ax.Position(4)]);
        % ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',-20,'width',400);
        e1 = computeEffectSize(v1a,v1b);
        g1 = e1.cliff_delta;
        e2 = computeEffectSize(v2a,v2b);
        g2 = e2.cliff_delta;
        e3 = computeEffectSize(v3a,v3b);
        g3 = e3.cliff_delta;

        v = [g1 g2 g3];
        for ii=1:3
            stem(ii,v(ii),'filled','MarkerFaceColor',cols(ii,:),'Color',cols(ii,:)); hold on;
        end

        % axis settings
        ax.YLim = [-1 1];
        ax.XLim = [0.5 3.5];
        ylabel(sprintf('Cliff''s Delta'))   
        ytickformat('%.1f');
        ax.FontSize = 8;
        ax.XTick = [];
        box off

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]);

        % stats
        vcoef = 1.1;
        fsiz = 12;
        dat = {v1a v1b;v2a v2b;v3a v3b};
        for ii = 1:3
            [z,p,~,~,~] = permutationZTest(dat{ii,1},dat{ii,2},'method','cliff');
            if p<.001
                text(ii,ax.YLim(2).*vcoef,'***','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<0.01
                text(ii,ax.YLim(2).*vcoef,'**','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<.05
                text(ii,ax.YLim(2).*vcoef,'*','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            end
        end

%% >>>>>>>>>> Within session correlation results, Hills
    xnow = 60;
    ynow = ynow-150;
    fname = [config.data_out_dir 'PIT_correlations_within_shuffles.mat'];
    load(fname,'shuffs');

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 350 120]);
        % ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',-20,'width',400);
    
        var = 'within_session_stability';
        v1a = clumaa.(var)(pidx & clumaa.partn==2,4); % arena 1 (edges)
        v2a = clumaa.(var)(pidx & clumaa.partn==2,5); % arena 1 (center)   
        v3a = clumaa.(var)(pidx & clumaa.partn==2,6); % arena 1 (corners)        
        v1b = squeeze(shuffs_within(:,2,2)); % edges shuffled
        v2b = squeeze(shuffs_within(:,3,2)); % center shuffled        
        v3b = squeeze(shuffs_within(:,4,2)); % corners shuffled

        xi = -1:0.01:1;
        bw = 0.05;
        f1a = ksdensity(v1a(:),xi(:),"Bandwidth",bw);
        f2a = ksdensity(v2a(:),xi(:),"Bandwidth",bw);
        f3a = ksdensity(v3a(:),xi(:),"Bandwidth",bw);
        f1b = ksdensity(v1b(:),xi(:),"Bandwidth",bw);
        f2b = ksdensity(v2b(:),xi(:),"Bandwidth",bw);
        f3b = ksdensity(v3b(:),xi(:),"Bandwidth",bw);

        % main plot
        alph = 0.7;  
        cols = winter(3);
        a1a = area(xi,f1a,'FaceColor',cols(1,:),'EdgeColor','k','FaceAlpha',alph); hold on;
        a2a = area(xi,f2a,'FaceColor',cols(2,:),'EdgeColor','k','FaceAlpha',alph);
        a3a = area(xi,f3a,'FaceColor',cols(3,:),'EdgeColor','k','FaceAlpha',alph);
        a1b = area(xi,f1b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--');
        a2b = area(xi,f2b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle',':');
        a3b = area(xi,f3b,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','-');

        % axis settings
        ax.XLim = [-0.5 1]; 
        ax.FontSize = 10;
        ytickformat('%.1f');
        box off
        xlabel(sprintf('Correlation ({\\itr})'))    
        ylabel(sprintf('PDF'))    

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        text(ax,0,1,maze_names{2},'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',10,'Color','k');

%% >>>>>>>>>> Within session correlation results (Cohen's d values) Hills
    xnow = xnow+400;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 120 ax.Position(4)]);
        % ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',-20,'width',400);
        e1 = computeEffectSize(v1a,v1b);
        g1 = e1.cliff_delta;
        e2 = computeEffectSize(v2a,v2b);
        g2 = e2.cliff_delta;
        e3 = computeEffectSize(v3a,v3b);
        g3 = e3.cliff_delta;

        v = [g1 g2 g3];
        for ii=1:3
            stem(ii,v(ii),'filled','MarkerFaceColor',cols(ii,:),'Color',cols(ii,:)); hold on;
        end

        % axis settings
        ax.YLim = [-1 1];
        ax.XLim = [0.5 3.5];
        ylabel(sprintf('Cliff''s Delta'))   
        ax.XTick = 1:4;
        ax.XTickLabel = {'Edges','Centre','Corners'};
        ytickformat('%.1f');
        ax.FontSize = 8;
        box off

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]);

        % stats
        vcoef = 1.1;
        fsiz = 12;
        dat = {v1a v1b;v2a v2b;v3a v3b};
        for ii = 1:3
            [z,p,~,~,~] = permutationZTest(dat{ii,1},dat{ii,2},'method','cliff');
            if p<.001
                text(ii,ax.YLim(2).*vcoef,'***','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<0.01
                text(ii,ax.YLim(2).*vcoef,'**','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<.05
                text(ii,ax.YLim(2).*vcoef,'*','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            end
        end

        % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S3.png']; 
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















