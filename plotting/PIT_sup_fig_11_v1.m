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
%% >>>>>>>>>> N fields in arena vs rep score in hills
    xnow = 50;
    ynow = 700;
    xsiz = 170;
    ysiz = 130;
    alph = 0.5;
    cmap = 'turbo';

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('A',sprintf(''),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);
    
        % collect data pairwise (this allows us to compare values from arena 1 and the hills)
        ucis = unique(clumaa.uci(pidx));
        v = NaN(length(ucis),2);
        var_1 = 'repetition_score'; % independent variable or x-axis
        var_2 = 'planar_fields'; % dependent variable or y-axis
        prts = [2 1]; % parts to collect data from
        cols = [1 1]; % columns to collect data from
        for uu = 1:length(ucis)
            v(uu,1) = clumaa.(var_1)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(1),cols(1));
            v(uu,2) = clumaa.(var_2)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(2),cols(2));
        end

        % main plot
        vplot = v;
        vplot(:,2) = vplot(:,2)+normrnd(0,0.3,size(v(:,2)));
        d = computeScatterDensity(vplot(:,1),vplot(:,2),'r1',64,'r2',64,'k',2);        
        scatter(vplot(:,1),vplot(:,2),30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',alph); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('N fields, %s',maze_names{1}))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];  
        ax.XLim = [-0.5 1.5]; 
        % ax2.YLim = [0 1];         
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.f')
        xtickformat('%.1f')        
        r = refline;
        set(r,'Color','k')
        ax.CLim = [min(d,[],'all','omitmissing') max(d,[],'all','omitmissing')];
        colormap(ax,cmap);

        % correlation
        [r,p] = corr(v(:,1),v(:,2),'rows','pairwise','type','Spearman');
        text(0,1.05,sprintf('{\\itr} = %.2f, {\\itp} = %.3f',r,p),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','k')
        % keyboard
    % inset bar graph
    axi = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+10 ax.Position(2) 60 ax.Position(4)]);
        var = 'repetition_score';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
        cutoff_rscore = prctile(v1,99);

        g1 = vplot(v(:,1)<=cutoff_rscore,2);
        g2 = vplot(v(:,1)>cutoff_rscore,2);

        % vectorise
        [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#0081A7'),hex2rgb('#F07167')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        axi.XTick = 1:2;
        axi.YTick = ax.YTick;
        axi.YLim = ax.YLim;        
        axi.XLim = [0.5 2.5]; 
        axi.FontSize = 8;
        box off
        axis off
        text(axi.XTick(1),axi.YLim(1),'Low','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)
        text(axi.XTick(2),axi.YLim(1),'High','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)

        % stats
        axt = axes('Units','pixels','Position',axi.Position,'Color','none');
            g1 = v(v(:,1)<=cutoff_rscore,2);
            g2 = v(v(:,1)>cutoff_rscore,2);
            [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data        
            axis off
            axt.XLim = axi.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'test','kw','plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

%% >>>>>>>>>> Field elongation in arena vs rep score in hills
    xnow = xnow+310;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('B',sprintf(''),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);
    
        % collect data pairwise (this allows us to compare values from arena 1 and the hills)
        ucis = unique(clumaa.uci(pidx));
        v = NaN(length(ucis),2);
        var_1 = 'repetition_score'; % independent variable or x-axis
        prts = [2 1]; % parts to collect data from
        cols = [1 1]; % columns to collect data from
        for uu = 1:length(ucis)
            v(uu,1) = clumaa.(var_1)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(1),cols(1));

            % get field data
            idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(2);
            fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
            if isempty(fdata)
                continue
            end
            e = sqrt(1 - (fdata.MinorAxisLength(:,1) ./ fdata.MajorAxisLength(:,1)));
            v(uu,2) = mean(e,'all','omitmissing');
        end

        % main plot
        d = computeScatterDensity(v(:,1),v(:,2),'r1',64,'r2',64,'k',2);
        scatter(v(:,1),v(:,2),30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',alph); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('Avg. elongation, %s',maze_names{1}))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];  
        ax.XLim = [-0.5 1.5]; 
        % ax2.YLim = [0 1];         
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.1f')
        xtickformat('%.1f')        
        r = refline;
        set(r,'Color','k')
        ax.CLim = [min(d,[],'all','omitmissing') max(d,[],'all','omitmissing')];
        colormap(ax,cmap);

        % correlation
        [r,p] = corr(v(:,1),v(:,2),'rows','pairwise','type','Spearman');
        text(0,1.05,sprintf('{\\itr} = %.2f, {\\itp} = %.3f',r,p),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','k')
        % return
    % inset bar graph
    axi = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+10 ax.Position(2) 60 ax.Position(4)]);
        g1 = v(v(:,1)<=cutoff_rscore,2);
        g2 = v(v(:,1)>cutoff_rscore,2);

        % vectorise
        [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#0081A7'),hex2rgb('#F07167')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        axi.XTick = 1:2;
        axi.YTick = ax.YTick;
        axi.YLim = ax.YLim;        
        axi.XLim = [0.5 2.5]; 
        axi.FontSize = 8;
        box off
        axis off
        text(axi.XTick(1),axi.YLim(1),'Low','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)
        text(axi.XTick(2),axi.YLim(1),'High','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)

        % stats
        axt = axes('Units','pixels','Position',axi.Position,'Color','none');
            axis off
            axt.XLim = axi.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'test','anova','plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

%% >>>>>>>>>> Anisotropy in arena vs rep score in hills
    xnow = 50;
    ynow = ynow-200;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('C',sprintf(''),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);
    
        % collect data pairwise (this allows us to compare values from arena 1 and the hills)
        ucis = unique(clumaa.uci(pidx));
        v = NaN(length(ucis),4);
        var_1 = 'repetition_score'; % independent variable or x-axis
        prts = [2 1]; % parts to collect data from
        cols = [1 1]; % columns to collect data from
        for uu = 1:length(ucis)
            v(uu,1) = clumaa.(var_1)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(1),cols(1));

            % get field data
            idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(2);
            fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
            if isempty(fdata)
                continue
            end
            rmap = clumaa.ratemap_planar{idx};                       

            % get field anisotropy
            fanisotropy = (fdata.BoundingBox(:,4) - fdata.BoundingBox(:,3)) ./  (fdata.BoundingBox(:,4) + fdata.BoundingBox(:,3)); % height - width / height + width

            % maze walls orientation
            b = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0]+0.5; 
            a = [90 0 90 0];
            gpoly = [];
            res = 100;
            for jj = 1:size(b,1)-1
                gpoly = [gpoly; linspace(b(jj,1),b(jj+1,1),res)' linspace(b(jj,2),b(jj+1,2),res)' repmat(a(jj),res,1)];
            end

            % distance to maze walls
            ds = NaN(size(fdata,1),2);
            for jj = 1:size(fdata,1)
                [i,d] = knnsearch(gpoly(:,[1:2]),[fdata.Centroid(jj,1),fdata.Centroid(jj,2)],'K',1);
                ds(jj,1:2) = [d gpoly(i,3)];              
            end

            % field to wall angle
            fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
            f1 = ds(:,2); % wall angles
            f2 = fdata.Orientation(:,1); % field orientations
            f1(f1==90) = f1(f1==90)-0.01;
            f2(f2==90) = f2(f2==90)-0.01;
            ftw = fun([f1 f2]);

            % accumulate
            v(uu,2) = mean(fanisotropy,'all','omitmissing');
            v(uu,3) = mean(ds(:,1),'all','omitmissing');
            v(uu,4) = mean(ftw,'all','omitmissing');            
        end

        % main plot
        d = computeScatterDensity(v(:,1),v(:,2),'r1',64,'r2',64,'k',2);
        scatter(v(:,1),v(:,2),30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',alph); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('Avg. anisotropy, %s',maze_names{1}))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];  
        ax.XLim = [-0.5 1.5]; 
        ax.YLim = [-1 1];         
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.1f')
        xtickformat('%.1f')        
        r = refline;
        set(r,'Color','k')
        ax.CLim = [min(d,[],'all','omitmissing') max(d,[],'all','omitmissing')];
        colormap(ax,cmap);

        % correlation
        [r,p] = corr(v(:,1),v(:,2),'rows','pairwise','type','Spearman');
        text(0,1.05,sprintf('{\\itr} = %.2f, {\\itp} = %.3f',r,p),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','k')
        
    % inset bar graph
    axi = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+10 ax.Position(2) 60 ax.Position(4)]);
        g1 = v(v(:,1)<=cutoff_rscore,2);
        g2 = v(v(:,1)>cutoff_rscore,2);

        % vectorise
        [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#0081A7'),hex2rgb('#F07167')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        axi.XTick = 1:2;
        axi.YTick = ax.YTick;
        axi.YLim = ax.YLim;        
        axi.XLim = [0.5 2.5]; 
        axi.FontSize = 8;
        box off
        axis off
        text(axi.XTick(1),axi.YLim(1),'Low','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)
        text(axi.XTick(2),axi.YLim(1),'High','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)

        % stats
        axt = axes('Units','pixels','Position',axi.Position,'Color','none');
            axis off
            axt.XLim = axi.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'test','anova','plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

%% >>>>>>>>>> Field area in arena vs rep score in hills
    xnow = xnow+310;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('D',sprintf(''),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);
    
        % collect data pairwise (this allows us to compare values from arena 1 and the hills)
        ucis = unique(clumaa.uci(pidx));
        vt = NaN(length(ucis),2);
        var_1 = 'repetition_score'; % independent variable or x-axis
        prts = [2 1]; % parts to collect data from
        cols = [1 1]; % columns to collect data from
        for uu = 1:length(ucis)
            vt(uu,1) = clumaa.(var_1)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(1),cols(1));

            % get field data
            idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(2);
            fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
            if isempty(fdata)
                continue
            end
            vt(uu,2) = mean(fdata.Area(:,1),'all','omitmissing');
        end

        % main plot
        vt(:,2) = vt(:,2) .* 0.0001; % convert to m2
        d = computeScatterDensity(vt(:,1),vt(:,2),'r1',64,'r2',64,'k',2);
        scatter(vt(:,1),vt(:,2),30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',alph); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('Avg. area, %s (m^{2})',maze_names{1}))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];  
        ax.XLim = [-0.5 1.5]; 
        % ax2.YLim = [0 1];         
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.1f')
        xtickformat('%.1f')        
        r = refline;
        set(r,'Color','k')
        ax.CLim = [min(d,[],'all','omitmissing') max(d,[],'all','omitmissing')];
        colormap(ax,cmap);

        % correlation
        [r,p] = corr(vt(:,1),vt(:,2),'rows','pairwise','type','Spearman');
        text(0,1.05,sprintf('{\\itr} = %.2f, {\\itp} = %.3f',r,p),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','k')
        
    % inset bar graph
    axi = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+10 ax.Position(2) 60 ax.Position(4)]);
        g1 = vt(vt(:,1)<=cutoff_rscore,2);
        g2 = vt(vt(:,1)>cutoff_rscore,2);

        % vectorise
        [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#0081A7'),hex2rgb('#F07167')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        axi.XTick = 1:2;
        axi.YTick = ax.YTick;
        axi.YLim = ax.YLim;        
        axi.XLim = [0.5 2.5]; 
        axi.FontSize = 8;
        box off
        axis off
        text(axi.XTick(1),axi.YLim(1),'Low','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)
        text(axi.XTick(2),axi.YLim(1),'High','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)

        % stats
        axt = axes('Units','pixels','Position',axi.Position,'Color','none');
            axis off
            axt.XLim = axi.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'test','anova','plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

%% >>>>>>>>>> Field distance from wall in arena vs rep score in hills
    xnow = 50;
    ynow = ynow-200;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('E',sprintf(''),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);
    
        % main plot
        d = computeScatterDensity(v(:,1),v(:,3),'r1',64,'r2',64,'k',2);
        scatter(v(:,1),v(:,3),30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',alph); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('Avg. distance, %s',maze_names{1}))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];  
        ax.XLim = [-0.5 1.5]; 
        % ax2.YLim = [0 1];         
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.f')
        xtickformat('%.1f')        
        r = refline;
        set(r,'Color','k')
        ax.CLim = [min(d,[],'all','omitmissing') max(d,[],'all','omitmissing')];
        colormap(ax,cmap);

        % correlation
        [r,p] = corr(v(:,1),v(:,3),'rows','pairwise','type','Spearman');
        text(0,1.05,sprintf('{\\itr} = %.2f, {\\itp} = %.3f',r,p),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','k')

    % inset bar graph
    axi = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+10 ax.Position(2) 60 ax.Position(4)]);
        g1 = v(v(:,1)<=cutoff_rscore,3);
        g2 = v(v(:,1)>cutoff_rscore,3);

        % vectorise
        [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#0081A7'),hex2rgb('#F07167')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        axi.XTick = 1:2;
        axi.YTick = ax.YTick;
        axi.YLim = ax.YLim;        
        axi.XLim = [0.5 2.5]; 
        axi.FontSize = 8;
        box off
        axis off
        text(axi.XTick(1),axi.YLim(1),'Low','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)
        text(axi.XTick(2),axi.YLim(1),'High','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)

        % stats
        axt = axes('Units','pixels','Position',axi.Position,'Color','none');
            axis off
            axt.XLim = axi.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'test','anova','plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

%% >>>>>>>>>> Field angle to wall in arena vs rep score in hills
    xnow = xnow+310;
    ynow = ynow;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('F',sprintf(''),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);
    
        % main plot
        d = computeScatterDensity(v(:,1),v(:,4),'r1',64,'r2',64,'k',2);
        scatter(v(:,1),v(:,4),30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',alph); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('Avg. angle, %s',maze_names{1}))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];  
        ax.XLim = [-0.5 1.5]; 
        ax.YLim = [0 90];   
        ax.YTick = 0:30:90;
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.f')
        xtickformat('%.1f')        
        r = refline;
        set(r,'Color','k')
        ax.CLim = [min(d,[],'all','omitmissing') max(d,[],'all','omitmissing')];
        colormap(ax,cmap);

        % correlation
        [r,p] = corr(v(:,1),v(:,4),'rows','pairwise','type','Spearman');
        text(0,1.05,sprintf('{\\itr} = %.2f, {\\itp} = %.3f',r,p),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','k')

    % inset bar graph
    axi = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+10 ax.Position(2) 60 ax.Position(4)]);
        g1 = v(v(:,1)<=cutoff_rscore,4);
        g2 = v(v(:,1)>cutoff_rscore,4);

        % vectorise
        [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#0081A7'),hex2rgb('#F07167')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        axi.XTick = 1:2;
        axi.YTick = ax.YTick;
        axi.YLim = ax.YLim;        
        axi.XLim = [0.5 2.5]; 
        axi.FontSize = 8;
        box off
        axis off
        text(axi.XTick(1),axi.YLim(1),'Low','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)
        text(axi.XTick(2),axi.YLim(1),'High','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)

        % stats
        axt = axes('Units','pixels','Position',axi.Position,'Color','none');
            axis off
            axt.XLim = axi.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'test','anova','plot_omnibus',1,'omnibus_text_y_gap_coeff',2);


%% >>>>>>>>>> Spatial info in arena vs rep score in hills
    xnow = 50;
    ynow = ynow-200;
    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('G',sprintf(''),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);
    
        % collect data pairwise (this allows us to compare values from arena 1 and the hills)
        ucis = unique(clumaa.uci(pidx));
        vt = NaN(length(ucis),2);
        var_1 = 'repetition_score'; % independent variable or x-axis
        var_2 = 'planar_spatial_info_shuffles';
        prts = [2 1]; % parts to collect data from
        cols = [1 1]; % columns to collect data from
        for uu = 1:length(ucis)
            vt(uu,1) = clumaa.(var_1)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(1),cols(1));
            vt(uu,2) = clumaa.(var_2)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(2),cols(2));
        end

        % main plot
        d = computeScatterDensity(vt(:,1),vt(:,2),'r1',64,'r2',64,'k',2);
        scatter(vt(:,1),vt(:,2),30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',alph); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('Spatial info, %s (b/s)',maze_names{1}))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];  
        ax.XLim = [-0.5 1.5]; 
        % ax2.YLim = [0 1];         
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.f')
        xtickformat('%.1f')        
        r = refline;
        set(r,'Color','k')
        ax.CLim = [min(d,[],'all','omitmissing') max(d,[],'all','omitmissing')];
        colormap(ax,cmap);

        % correlation
        [r,p] = corr(vt(:,1),vt(:,2),'rows','pairwise','type','Spearman');
        text(0,1.05,sprintf('{\\itr} = %.2f, {\\itp} = %.3f',r,p),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','k')
        
    % inset bar graph
    axi = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+10 ax.Position(2) 60 ax.Position(4)]);
        g1 = vt(vt(:,1)<=cutoff_rscore,2);
        g2 = vt(vt(:,1)>cutoff_rscore,2);

        % vectorise
        [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#0081A7'),hex2rgb('#F07167')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        axi.XTick = 1:2;
        axi.YTick = ax.YTick;
        axi.YLim = ax.YLim;        
        axi.XLim = [0.5 2.5]; 
        axi.FontSize = 8;
        box off
        axis off
        text(axi.XTick(1),axi.YLim(1),'Low','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)
        text(axi.XTick(2),axi.YLim(1),'High','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)

        % stats
        axt = axes('Units','pixels','Position',axi.Position,'Color','none');
            axis off
            axt.XLim = axi.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'test','anova','plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

%% >>>>>>>>>> Stability vs repetition hills
    xnow = xnow+310;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('H',sprintf(''),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);
    
        ucis = unique(clumaa.uci(pidx));
        v = NaN(length(ucis),2);
        var_1 = 'repetition_score'; % independent variable or x-axis
        var_2 = 'within_session_stability'; % dependent variable or y-axis
        prts = [2 2]; % parts to collect data from
        cols = [1 1]; % columns to collect data from
        for uu = 1:length(ucis)
            v(uu,1) = clumaa.(var_1)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(1),cols(1));
            v(uu,2) = clumaa.(var_2)(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==prts(2),cols(2));
        end

        % main plot
        d = computeScatterDensity(v(:,1),v(:,2),'r1',64,'r2',64,'k',2);
        scatter(v(:,1),v(:,2),30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',alph); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('Within session stability, %s',maze_names{2}))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];  
        % ax.XLim = [-0.8 0.8]; 
        % ax.YLim = [-1 0.8];         
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.1f')
        xtickformat('%.1f')        
        r = refline;
        set(r,'Color','k')
        ax.CLim = [min(d,[],'all','omitmissing') max(d,[],'all','omitmissing')];
        colormap(ax,cmap);

        % correlation
        [r,p] = corr(v(:,1),v(:,2),'rows','pairwise','type','Spearman');
        text(0,1.05,sprintf('{\\itr} = %.2f, {\\itp} = %.3f',r,p),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','k')
        
    % inset bar graph
    axi = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+10 ax.Position(2) 60 ax.Position(4)]);
        g1 = v(v(:,1)<=cutoff_rscore,2);
        g2 = v(v(:,1)>cutoff_rscore,2);

        % vectorise
        [ds,gs] = vectorDATAGROUP([],g1,g2); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#0081A7'),hex2rgb('#F07167')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        axi.XTick = 1:2;
        axi.YTick = ax.YTick;
        axi.YLim = ax.YLim;        
        axi.XLim = [0.5 2.5]; 
        axi.FontSize = 8;
        box off
        axis off
        text(axi.XTick(1),axi.YLim(1),'Low','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)
        text(axi.XTick(2),axi.YLim(1),'High','HorizontalAlignment','center','VerticalAlignment','top','FontSize',axi.FontSize)

        % stats
        axt = axes('Units','pixels','Position',axi.Position,'Color','none');
            axis off
            axt.XLim = axi.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'test','anova','plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

            % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S11.png']; 
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







