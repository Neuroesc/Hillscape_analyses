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
    warning('off','MATLAB:legend:IgnoringExtraEntries')

    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example cell (fields aligned to maze)
    xnow = 30;
    ynow = 700;
    
    uci = 'RG26_230119_t3_c5';
    m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
    m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
    m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+80 120 70]);
        ah = add_panel_title('a',sprintf('Example cell'),'yoffset',-10,'xoffset',10,'width',400);  
    
        imagesc(m1,'alphadata',~isnan(m1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,turbo)   
        ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')

    % plot hills
    ax2 = axes('Units','pixels','Position',[xnow ynow 120 70]);
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
        text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+10 ynow+70 70 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f Hz',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')

    % aspect ratio analysis (if needed)
    if ~any(ismember(clumaa.Properties.VariableNames,'planar_field_aspect_ratio')) % if the column(s) do not exist yet
        clumaa.planar_field_aspect_ratio = (clumaa.planar_field_info(:,4) - clumaa.planar_field_info(:,5)) ./ (clumaa.planar_field_info(:,4) + clumaa.planar_field_info(:,5)) ; % aspect ratio = field height / width
        clumaa.surficial_field_aspect_ratio = (clumaa.surficial_field_info(:,4) - clumaa.surficial_field_info(:,5)) ./ (clumaa.surficial_field_info(:,4) + clumaa.surficial_field_info(:,5)) ; % aspect ratio = field height / width
    end  

%% >>>>>>>>>> elongation
    xnow = xnow+170;
    ynow = ynow;

    planar_field_elongation = sqrt(1 - (clumaa.planar_field_info(:,3) ./ clumaa.planar_field_info(:,2)));
    surficial_field_elongation = sqrt(1 - (clumaa.surficial_field_info(:,3) ./ clumaa.surficial_field_info(:,2)));

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow+15 90 125]);
        ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',0,'width',400);  
        v1 = planar_field_elongation(pidx & clumaa.partn==1,1); % arena 1 data
        v2a = planar_field_elongation(pidx & clumaa.partn==2,1); % hills data
        v2b = surficial_field_elongation(pidx & clumaa.partn==2,1); % hills data
        v3 = planar_field_elongation(pidx & clumaa.partn==3,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'meanmarker',{'o'},'meanfill',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',20,'dotalpha',0.5,'dotsmarker',plot_set(2,[1 2 2 3]),'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Eccentricity'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        % ax.YLim = [-1 1];
        ax.FontSize = 8;
        ax.XTickLabelRotation = 25;
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s      \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);

            for ii = 1:max(gs(:))
                [~,p,~] = ttest(ds(gs==ii),0); % compare group to zero
                if p<=.001
                    text(ii,0,sprintf('***'),'FontSize',10,'HorizontalAlignment','center','Rotation',0,'VerticalAlignment','bottom','Units','data')
                elseif p<=0.01
                    text(ii,0,sprintf('**'),'FontSize',10,'HorizontalAlignment','center','Rotation',0,'VerticalAlignment','bottom','Units','data')
                elseif p<.05
                    text(ii,0,sprintf('*'),'FontSize',10,'HorizontalAlignment','center','Rotation',0,'VerticalAlignment','bottom','Units','data')
                end
            end

%% >>>>>>>>>> aspect ratio
    xnow = xnow+140;
    ynow = ynow;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow+15 90 125]);
        ah = add_panel_title('c',sprintf(''),'yoffset',0,'xoffset',0,'width',400);  
        var = 'planar_field_aspect_ratio';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
        v2a = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data
        v2b = clumaa.surficial_field_aspect_ratio(pidx & clumaa.partn==2,1); % hills data
        v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'meanmarker',{'o'},'meanfill',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',20,'dotalpha',0.5,'dotsmarker',plot_set(2,[1 2 2 3]),'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Aspect ratio'))    
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        ax.YLim = [-1 1];
        ax.FontSize = 8;
        ax.XTickLabelRotation = 25;
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s      \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);

            for ii = 1:max(gs(:))
                [~,p,~] = ttest(ds(gs==ii),0); % compare group to zero
                if p<=.001
                    text(ii,0,sprintf('***'),'FontSize',10,'HorizontalAlignment','center','Rotation',0,'VerticalAlignment','bottom','Units','data')
                elseif p<=0.01
                    text(ii,0,sprintf('**'),'FontSize',10,'HorizontalAlignment','center','Rotation',0,'VerticalAlignment','bottom','Units','data')
                elseif p<.05
                    text(ii,0,sprintf('*'),'FontSize',10,'HorizontalAlignment','center','Rotation',0,'VerticalAlignment','bottom','Units','data')
                end
            end

%% >>>>>>>>>> Average fields
    xnow = xnow + 145;

    var = 'planar_amap';
    v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
    all_amaps_a1 = cat(3,v1{:});
    a1 = mean(all_amaps_a1,3,'omitnan');

    v2a = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data (planar)
    all_amaps_a2a = cat(3,v2a{:});
    a2a = mean(all_amaps_a2a,3,'omitnan');

    v2b = clumaa.surficial_amap(pidx & clumaa.partn==2,1); % hills data (surficial)
    v2b = cellfun(@(x) imresize(x,[95 230],'bilinear'),v2b,'UniformOutput',false);
    all_amaps_a2b = cat(3,v2b{:});
    a2b = mean(all_amaps_a2b,3,'omitnan');

    v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 2 data
    all_amaps_a3 = cat(3,v3{:});
    a3 = mean(all_amaps_a3,3,'omitnan');

    mapset.binsize = 32; % (mm) firing rate map bin size

    ds = size(a1).*mapset.binsize./2;

    % create axis
    siz = 60;
    xbuff = 10;
    ax1 = axes('Units','pixels','Position',[xnow ynow+80 siz siz]); % arena 1
        ah = add_panel_title('d',sprintf('Avg. field shape'),'yoffset',0,'xoffset',0,'width',400);  
    
        imagesc([-ds(2) ds(2)],[-ds(1) ds(1)],a1);
        daspect([1 1 1])
        colormap(turbo)
        ax1.CLim = [0.1 1];

        b = 800;
        ax1.XLim = [-b b];
        ax1.YLim = [-b b];
        ax1.XTick = [];
        ax1.FontSize = 7;
        ylabel('mm')
        ax1.YTick = [-800 0 800];
        text(0,1,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom')

    ax2 = axes('Units','pixels','Position',[xnow ynow siz siz]); % arena 2
        imagesc([-ds(2) ds(2)],[-ds(1) ds(1)],a3);
        daspect([1 1 1])
        colormap(turbo)
        ax2.CLim = ax1.CLim;

        ax2.XLim = [-b b];
        ax2.YLim = [-b b];
        ylabel('mm')
        xlabel('mm')
        ax2.XTick = [-600 0 600];
        ax2.YTick = [-600 0 600];
        ax2.XTickLabelRotation = 0;        
        ax2.FontSize = 7;
        text(0,1,sprintf('%s',maze_names{3}),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom')

    ax3 = axes('Units','pixels','Position',[xnow+ax1.Position(3)+xbuff ynow+80 siz siz]); % hills
        imagesc([-ds(2) ds(2)],[-ds(1) ds(1)],a2a);
        daspect([1 1 1])
        colormap(turbo)
        ax3.CLim = ax1.CLim;

        b = 800;
        ax3.XLim = [-b b];
        ax3.YLim = [-b b];
        ax3.XTick = [];
        ax3.YTick = [];
        ax3.FontSize = 7;
        text(0,1,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom')

    ax4 = axes('Units','pixels','Position',[xnow+ax1.Position(3)+xbuff ynow siz siz]); % hills flattened
        imagesc([-ds(2) ds(2)],[-ds(1) ds(1)],a2b);
        daspect([1 1 1])
        colormap(turbo)
        ax4.CLim = ax1.CLim;

        ax4.XLim = [-b b];
        ax4.YLim = [-b b];
        ax4.YTick = [];
        xlabel('mm')
        ax4.XTick = [-600 0 600];
        ax4.XTickLabelRotation = 0;
        ax4.FontSize = 7;
        text(0,1,sprintf('%s (flattened)',maze_names{2}),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+30 ynow-40 70 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.f ',0),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.f',1),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,0,sprintf('Correlation'),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> Example cells (fields aligned to maze)
% cell 1
    xnow = 30;
    ynow = ynow-220;
    
    uci = 'RG26_230119_t2_c15';
    m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
    m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
    m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+110 120 70]);
        ah = add_panel_title('e',sprintf('Example cell'),'yoffset',-15,'xoffset',10,'width',400);  
    
        imagesc(m1,'alphadata',~isnan(m1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,turbo)   
        ax1.CLim = [0 max([max(m1(:),[],'omitnan') 1])];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
        % text(0,1,sprintf('Cell 1'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

        % ellipses
        f = clumaa.planar_field_data{(ismember(clumaa.uci,uci) & clumaa.partn==1)};
        for ff = 1:size(f,1)
            e = ellipse(f.MajorAxisLength(ff,1)./2,f.MinorAxisLength(ff,1)./2,-deg2rad(f.Orientation(ff,1)),f.Centroid(ff,1),f.Centroid(ff,2));
            e.LineWidth = 1;
            e.Color = 'w';
            plot(f.Centroid(ff,1),f.Centroid(ff,2),'w+','MarkerSize',10);
        end

    % plot hills
    ax2 = axes('Units','pixels','Position',[xnow ynow+30 120 70]);
        imagesc(m2,'alphadata',~isnan(m2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,turbo)   
        ax2.CLim = [0 max([max(m2(:),[],'omitnan') 1])];
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(m2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'w:')

        % text
        text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        text(0,0,sprintf('%.1f',ax2.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')

        % ellipses
        f = clumaa.planar_field_data{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
        for ff = 1:size(f,1)
            e = ellipse(f.MajorAxisLength(ff,1)./2,f.MinorAxisLength(ff,1)./2,-deg2rad(f.Orientation(ff,1)),f.Centroid(ff,1),f.Centroid(ff,2));
            e.LineWidth = 1;
            e.Color = 'w';
            plot(f.Centroid(ff,1),f.Centroid(ff,2),'w+','MarkerSize',10);
        end

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+10 ynow+100 70 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f Hz',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> All fields
    % collect data
    ucis = unique(clumaa.uci(pidx));
    datn = cell(1,3);
    for pp=1:3    
        dat = [];
        for uu = 1:length(ucis) % for each cell
            idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==pp;
            if pp==2
                fdata = clumaa.surficial_field_data{idx}; % field data for this cell in this session in this part
                rmap = clumaa.ratemap_surficial{idx};       
            else
                fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
                rmap = clumaa.ratemap_planar{idx};                       
            end
            if isempty(fdata)
                continue
            end

            if pp==2
                scale_x = 116 ./ size(rmap,2); % we need to scale field centroids to match a common map size
            else
                scale_x = 96 ./ size(rmap,2); % we need to scale field centroids to match a common map size
            end

            b = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0]+0.5; 
            a = [90 0 90 0];
            gpoly = [];
            res = 100;
            for jj = 1:size(b,1)-1
                gpoly = [gpoly; linspace(b(jj,1),b(jj+1,1),res)' linspace(b(jj,2),b(jj+1,2),res)' repmat(a(jj),res,1)];
            end

            ps = linspace(0,size(rmap,2),7);
            tops = ps(2:2:end);
            b2 = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0;tops(1) 0;tops(1) size(rmap,1);tops(2) size(rmap,1);tops(2) 0;tops(3) 0;tops(3) size(rmap,1)]+0.5; 
            a2 = [90 0 90 0 0 90 0 90 0 90];
            gpoly2 = [];
            res = 100;
            for jj = 1:size(b2,1)-1
                gpoly2 = [gpoly2; linspace(b2(jj,1),b2(jj+1,1),res)' linspace(b2(jj,2),b2(jj+1,2),res)' repmat(a2(jj),res,1)];
            end

            ds = NaN(size(fdata,1),4);
            for jj = 1:size(fdata,1)
                [i,d] = knnsearch(gpoly(:,[1:2]),[fdata.Centroid(jj,1),fdata.Centroid(jj,2)],'K',1);
                ds(jj,1:2) = [d gpoly(i,3)];

                [i,d] = knnsearch(gpoly2(:,[1:2]),[fdata.Centroid(jj,1),fdata.Centroid(jj,2)],'K',1);
                ds(jj,3:4) = [d gpoly2(i,3)];                
            end
            dat = [dat; fdata.MajorAxisLength(:,1) fdata.MinorAxisLength(:,1) fdata.Orientation(:,1) fdata.WeightedCentroid(:,1).*scale_x fdata.WeightedCentroid(:,2) fdata.Area(:,1) ds(:,1) repmat(uu,size(fdata,1),1) ds(:,2) ds(:,3:4)];
        end
        datn(pp) = { dat };
    end

    xnow = xnow+135;
    ynow = ynow-10;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+50 240 90]);
        ah = add_panel_title('f',sprintf('Example fields'),'yoffset',10,'xoffset',35,'width',400);  
        % ellipses
        f = datn{1};
        alph = 0.3;
        nfields = 150;
        rng(999); % for reproducilibity
        rindx = randperm(size(f,1),nfields);
        for ff = 1:nfields
            e = ellipse(f(rindx(ff),1)./2,f(rindx(ff),2)./2,-deg2rad(f(rindx(ff),3)),f(rindx(ff),4),f(rindx(ff),5),'none');
            x = e.XData;
            y = e.YData;
            patch(x,y,-f(ff,3),'EdgeColor','none','FaceAlpha',alph,'Clipping','off'); hold on;
            delete(e);
        end

        axis xy on
        box on
        daspect([1 1 1])
        c = cmocean('thermal');
        % c = viridis(64);
        colormap(ax1,[c;flipud(c)])   
        ax1.CLim = [-90 90];
        ax1.XLim = [0.5 size(posdata.anisotropy_map{1,1},2)+0.5];
        ax1.YLim = [0.5 size(posdata.anisotropy_map{1,1},1)+0.5];        
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,1.1,sprintf('%s, %d random fields',maze_names{1},nfields),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+220 ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
        % ellipses
        f = datn{2};
        for ff = 1:nfields
            e = ellipse(f(rindx(ff),1)./2,f(rindx(ff),2)./2,-deg2rad(f(rindx(ff),3)),f(rindx(ff),4),f(rindx(ff),5),'none');
            x = e.XData;
            y = e.YData;
            patch(x,y,-f(ff,3),'EdgeColor','none','FaceAlpha',alph,'Clipping','off'); hold on;
            delete(e);
        end

        axis xy on
        box on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)   
        ax2.CLim = ax1.CLim;        
        ax2.XLim = [0.5 size(posdata.anisotropy_map{1,2},2)+0.5];
        ax2.YLim = [0.5 size(posdata.anisotropy_map{1,2},1)+0.5]; 
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(ax2.XLim(1),ax2.XLim(2),7);
        tops = ps(2:2:end);
        pt = plot([tops; tops],ax2.YLim,'k--');

        % text
        text(0,1.1,sprintf('%s, %d random fields',maze_names{2},nfields),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+50 ynow+35 90 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,0,sprintf('Field orientation (%c)',176),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',axc.FontSize,'Units','normalized')

        [~,leg] = legendflex(pt,{'Hilltops'},'anchor',{'s','s'},'ncol',1,'box','off','buffer',[0,-20],'xscale',1,'FontSize',8);   

%% >>>>>>>>>> orientation
    xnow = 40;
    ynow = ynow-155;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow+75 60 60]);
        ah = add_panel_title('g',sprintf(''),'yoffset',0,'xoffset',0,'width',400); 
        v1 = [-90 90];
        v2 = [0];

        % get circular density of data
        ri = linspace(-90,90,360);
        xi = linspace(-180,180,360); 
        k = 0.8;
        f1 = circ_ksdensity(deg2rad(v1).*2,deg2rad(xi),[],k);
        f2 = circ_ksdensity(deg2rad(v2).*2,deg2rad(xi),[],k);

        % plot data
        alph = 0.5;
        a1 = area(ri,f1,'FaceColor','k','EdgeColor','none','FaceAlpha',1); hold on;

        % axis settings
        ax.XLim = [-90 90];
        ax.XTick = -180:90:180;
        ax.YTick = [];
        ax.FontSize = 7;        
        box off
        ylabel('PDF')
        ytickformat('%.1f');
        text(0.5,1.0,sprintf('Fields aligned to y-axis'),'Units','normalized','HorizontalAlignment','center','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)
        text(0.5,1.18,sprintf('Prediction'),'Units','normalized','HorizontalAlignment','center','FontSize',9,'Color','k','VerticalAlignment','bottom','rotation',0)

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow-15 60 60]);
        % plot data
        alph = 0.5;
        a2 = area(ri,f2,'FaceColor','k','EdgeColor','none','FaceAlpha',1); hold on;

        % axis settings
        ax.XLim = [-90 90];
        ax.XTick = -180:90:180;
        ax.YTick = [];  
        ax.FontSize = 7;
        box off
        xlabel(sprintf('Field orientation (%c)',176))
        ylabel('PDF')
        ytickformat('%.1f');
        text(0.5,1.0,sprintf('Fields aligned to x-axis'),'Units','normalized','HorizontalAlignment','center','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)

%% >>>>>>>>>> orientation analysis results        
    ax = axes('Units','pixels','Position',[xnow+115 ynow+40 150 100]);
        % ah = add_panel_title('g',sprintf(''),'yoffset',0,'xoffset',-10,'width',400);  
    
        % var = 'planar_field_info';
        % v1 = clumaa.(var)(pidx & clumaa.partn==1,6); % arena 1 data
        % v2 = clumaa.(var)(pidx & clumaa.partn==2,6); % hills data
        % v3 = clumaa.(var)(pidx & clumaa.partn==3,6); % arena 2 data

        % arena 1
        v1 = datn{1}(:,3);
        v2 = datn{2}(:,3);
        v3 = datn{3}(:,3);

        % get circular density of data
        k = 0.25;
        f1 = circ_ksdensity(deg2rad(v1).*2,deg2rad(xi),[],k);
        f2 = circ_ksdensity(deg2rad(v2).*2,deg2rad(xi),[],k);
        f3 = circ_ksdensity(deg2rad(v3).*2,deg2rad(xi),[],k);

        % plot data
        a1 = area(ri,f1,'FaceColor',plot_set{1,1},'EdgeColor','none','FaceAlpha',alph); hold on;
        a2 = area(ri,f2,'FaceColor',plot_set{1,2},'EdgeColor','none','FaceAlpha',alph); hold on;

        % axis settings
        ax.XLim = [-90 90];
        ax.XTick = -180:45:180;
        box off
        xlabel(sprintf('Field orientation (%c)',176))
        ylabel('PDF')
        ytickformat('%.1f');
        text(0.5,1,sprintf('Real data'),'Units','normalized','HorizontalAlignment','center','FontSize',9,'Color','k','VerticalAlignment','bottom','rotation',0)

        % stats
        rv1 = circ_r(deg2rad(v1).*2); % angle double
        rv2 = circ_r(deg2rad(v2).*2); % angle double
        [~,leg] = legendflex([a1 a2],{sprintf('Arena 1: Rayleigh v = %.2f',rv1),sprintf('Hills: Rayleigh v = %.2f',rv2)},'anchor',{'sw','sw'},'ncol',1,'box','off','buffer',[0,-70],'xscale',0.5,'FontSize',8);   
        leg(3).Children.FaceAlpha = alph;
        leg(4).Children.FaceAlpha = alph;
% 
% %% >>>>>>>>>> orientation analysis schematic
%     xnow = xnow+200;
%     ynow = ynow;
% 
%     uci = 'RG26_230119_t2_c15';
%     m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
%     m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
%     m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};
% 
%     % plot arena
%     ax1 = axes('Units','pixels','Position',[xnow ynow+100 120 70]);
%         ah = add_panel_title('g',sprintf(''),'yoffset',-10,'xoffset',0,'width',400);  
% 
%         patch([0 0 size(m1,2) size(m1,2)],[0 size(m1,1) size(m1,1) 0],'k','FaceColor','none','EdgeColor','k'); hold on;
%         axis xy off tight
%         daspect([1 1 1])
%         ax1.XTick = [];
%         ax1.YTick = [];
% 
%         % text
%         text(0.5,1,sprintf('Aligned to closest wall'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
% 
%         % ellipse 1
%         a = [15 5];
%         c = [40 35];
%         o = deg2rad(0);
%         e = ellipse(a(1),a(2),o,c(1),c(2));
%         e.LineWidth = 1;
%         e.Color = 'k';
%         line([c(1)-(a(1)) c(1)+(a(1))],[c(2) c(2)],'Color','k'); % major axis
%         line([c(1)-(a(1)) c(1)+(a(1))],[size(m1,1) size(m1,1)],'Color','r'); % wall segment
%         line([c(1)-(a(1)) c(1)-(a(1))],[c(2) size(m1,1)],'Color','k','LineStyle',':');
%         line([c(1)+(a(1)) c(1)+(a(1))],[c(2) size(m1,1)],'Color','k','LineStyle',':');
% 
%         % ellipse 2
%         a = [15 5];
%         c = [85 25];
%         o = deg2rad(90);
%         e = ellipse(a(1),a(2),o,c(1),c(2));
%         e.LineWidth = 1;
%         e.Color = 'k';
%         line([c(1) c(1)],[c(2)-(a(1)) c(2)+(a(1))],'Color','k'); % major axis
%         line([size(m1,2) size(m1,2)],[c(2)-(a(1)) c(2)+(a(1))],'Color','r'); % wall segment
%         line([c(1) size(m1,2)],[c(2)-(a(1)) c(2)-(a(1))],'Color','k','LineStyle',':');
%         line([c(1) size(m1,2)],[c(2)+(a(1)) c(2)+(a(1))],'Color','k','LineStyle',':');
% 
%     % plot hills
%     ax2 = axes('Units','pixels','Position',[xnow ynow+10 120 70]);
%         patch([0 0 size(m2,2) size(m2,2) 0],[0 size(m2,1) size(m2,1) 0 0],'k','FaceColor','none','EdgeColor','k'); hold on;
%         axis xy off tight
%         daspect([1 1 1])
%         ax2.XTick = [];
%         ax2.YTick = [];
% 
%         % text
%         text(0.5,1,sprintf('Aligned to hills'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
% 
%         % ellipse 1
%         a = [15 5];
%         c = [25 25];
%         o = deg2rad(90);
%         e = ellipse(a(1),a(2),o,c(1),c(2));
%         e.LineWidth = 1;
%         e.Color = 'k';
%         line([c(1) c(1)],[c(2)-(a(1)) c(2)+(a(1))],'Color','k');
%         scatter(c(1),size(m1,1),5,'r');
%         line([c(1) c(1)],[c(2) size(m1,1)],'Color','k','LineStyle',':');
% 
%         % ellipse 2
%         a = [15 5];
%         c = [85 25];
%         o = deg2rad(90);
%         e = ellipse(a(1),a(2),o,c(1),c(2));
%         e.LineWidth = 1;
%         e.Color = 'k';
%         line([c(1) c(1)],[c(2)-(a(1)) c(2)+(a(1))],'Color','k'); % major axis
%         line([size(m1,2) size(m1,2)],[c(2)-(a(1)) c(2)+(a(1))],'Color','r'); % wall segment
%         line([c(1) size(m1,2)],[c(2)-(a(1)) c(2)-(a(1))],'Color','k','LineStyle',':');
%         line([c(1) size(m1,2)],[c(2)+(a(1)) c(2)+(a(1))],'Color','k','LineStyle',':');
% 
%         % additional plots
%         ps = linspace(ax2.XLim(1),ax2.XLim(2),7);
%         tops = ps(2:2:end);
%         plot([tops; tops],ax2.YLim,'Color',[.5 .5 .5],'LineStyle','-')

%% >>>>>>>>>> orientation analysis predictions
    xnow = xnow+310;
    ynow = ynow;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow+75 60 60]);
        ah = add_panel_title('h',sprintf(''),'yoffset',0,'xoffset',0,'width',400); 
    
        fun = @(x) rad2deg( atan( abs( (tand(x(:,1))-tand(x(:,2))) ./ (1+tand(x(:,1)).*tand(x(:,2))) ) ) ); % angle between two lines given angle from x-axis, output in degs
        v1 = fun([90-0.01 90-0.01;0 0;0 0]); % aligned to walls
        v2 = fun([90-0.01 90-0.01;90-0.01 0]); % aligned to hills

        % get circular density of data
        k = 0.20;
        xi = linspace(0,90,360);
        k1 = circ_ksdensity(deg2rad(v1),deg2rad(xi),[],k);
        k2 = circ_ksdensity(deg2rad(v2),deg2rad(xi),[],k);

        % plot data
        a1 = area(xi,k1,'FaceColor','k','EdgeColor','none','FaceAlpha',1); hold on;

        % axis settings
        ax.XLim = [min(xi,[],'all') max(xi,[],'all')];
        ax.XTick = min(xi,[],'all'):45:max(xi,[],'all');
        ax.YTick = [];
        ax.FontSize = 7;        
        box off
        ylabel('PDF')
        ytickformat('%.1f');
        text(0.5,1.0,sprintf('Fields aligned to walls'),'Units','normalized','HorizontalAlignment','center','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)
        text(0.5,1.18,sprintf('Prediction'),'Units','normalized','HorizontalAlignment','center','FontSize',9,'Color','k','VerticalAlignment','bottom','rotation',0)

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow-15 60 60]);
        % plot data
        a2 = area(xi,k2,'FaceColor','k','EdgeColor','none','FaceAlpha',1); hold on;

        % axis settings
        ax.XLim = [min(xi,[],'all') max(xi,[],'all')];
        ax.XTick = min(xi,[],'all'):45:max(xi,[],'all');
        ax.YTick = [];  
        ax.FontSize = 7;
        box off
        xlabel(sprintf('Field to wall angle (%c)',176))
        ylabel('PDF')
        ytickformat('%.1f');
        text(0.5,1.0,sprintf('Fields aligned to hills'),'Units','normalized','HorizontalAlignment','center','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',0)

%% >>>>>>>>>> wall alignment analysis results        
    % create axis
    ax = axes('Units','pixels','Position',[xnow+115 ynow+40 150 100]);            
        % arena 1
        f1 = datn{1};
        f1(f1(:,3)==90,3) = f1(f1(:,3)==90,3)-0.01;
        f1(f1(:,9)==90,9) = f1(f1(:,9)==90,9)-0.01;        
        v1 = fun(f1(:,[3 9]));

        % hills
        f2 = datn{2};
        f2(f2(:,3)==90,3) = f2(f2(:,3)==90,3)-0.01;
        f2(f2(:,9)==90,9) = f2(f2(:,9)==90,9)-0.01;  
        f2(f2(:,11)==90,11) = f2(f2(:,11)==90,11)-0.01;                 
        v2 = fun(f2(:,[3 9]));

        % hills (if hill tops are boundaries)
        v3 = fun(f2(:,[3 11]));

        % get circular density of data
        k = 0.05;
        xi = linspace(0,90,360);         
        k1 = circ_ksdensity(deg2rad(v1),deg2rad(xi),[],k);
        k2 = circ_ksdensity(deg2rad(v2),deg2rad(xi),[],k);
        k3 = circ_ksdensity(deg2rad(v3),deg2rad(xi),[],k);
        
% keyboard
        % plot data
        alph = 0.5;
        a1 = area(xi,k1,'FaceColor',plot_set{1,1},'EdgeColor','none','FaceAlpha',alph); hold on;
        a2 = area(xi,k2,'FaceColor',plot_set{1,2},'EdgeColor','none','FaceAlpha',alph); hold on;
        a3 = plot(xi,k3,'k','LineStyle','--');

        % axis settings
        ax.XLim = [min(xi,[],'all') max(xi,[],'all')];
        ax.XTick = min(xi,[],'all'):45:max(xi,[],'all');
        box off
        xlabel(sprintf('Field to wall angle (%c)',176))
        ylabel('PDF')
        ytickformat('%.1f');
        text(0.5,1,sprintf('Real data'),'Units','normalized','HorizontalAlignment','center','FontSize',9,'Color','k','VerticalAlignment','bottom','rotation',0)

        % stats
        % rv1 = circ_r(deg2rad(v1).*4); % angle quadruple
        % rv2 = circ_r(deg2rad(v2).*4); % angle quadruple
        [~,leg] = legendflex([a1 a2 a3],{maze_names{1},maze_names{2},'Hills (hilltops = walls)'},'anchor',{'sw','sw'},'ncol',1,'box','off','buffer',[0,-85],'xscale',0.5,'FontSize',8);   
        leg(4).Children.FaceAlpha = alph;
        leg(5).Children.FaceAlpha = alph;
% keyboard

%% >>>>>>>>>> Field anisotropy analysis (if needed)
    if ~any(ismember(clumaa.Properties.VariableNames,'field_anisotropy')) % if the column(s) do not exist yet
        clumaa.field_anisotropy = cell(size(clumaa,1),1);

        for ss = 1:size(posdata,1) % for each recording session
            disp(sprintf('\tSession %d of %d (%.f%%)',ss,size(posdata,1),ss/size(posdata,1)*100))
            sidx = ismember(clumaa.pos_idx(:),posdata.pos_idx(ss));
            session_times = posdata.session_times{ss};
    
            ucis = unique(clumaa.uci(sidx));
    
            for pp = 1:size(session_times,1) % for each part
                for uu = 1:length(ucis) % for each cell recorded in this session
                    idx = ismember(clumaa.uci,ucis{uu}) & clumaa.partn==pp;
                    if pp==2
                        fdata = clumaa.surficial_field_data{idx}; % field data for this cell in this session in this part
                        rmap = clumaa.ratemap_surficial{idx};
                        bmap = posdata.anisotropy_map{ss,pp};
                    else
                        fdata = clumaa.planar_field_data{idx}; % field data for this cell in this session in this part
                        rmap = clumaa.ratemap_planar{idx};       
                        bmap = posdata.anisotropy_map{ss,pp};
                    end
                    if isempty(fdata)
                        continue
                    end
    
                    % get field anisotropy
                    fanisotropy = (fdata.BoundingBox(:,4) - fdata.BoundingBox(:,3)) ./  (fdata.BoundingBox(:,4) + fdata.BoundingBox(:,3)); % height - width / height + width
    
                    % get behaviour anisotropy
                    % fcents = fdata(:,2:3); % centroids
                    fcents = fdata.WeightedCentroid(:,1:2); % weighted centroids
                    % banisotropy = interp2(bmap,fcents(:,1),fcents(:,2),'nearest');
    
                    banisotropy = NaN(size(fcents,1),1);
                    for ff = 1:size(fdata,1) % for every field
                        pxls = fdata.PixelIdxList{ff,1};
                        banisotropy(ff) = median(reshape(bmap(pxls),[],1),'omitnan');
                    end
                    
                    % accumulate
                    if pp==2
                        scale_x = 116 ./ size(bmap,2); % we need to scale field centroids to match a common map size
                    else
                        scale_x = 1;
                    end
                    clumaa.field_anisotropy(idx) = { single([fanisotropy banisotropy fdata.WeightedCentroid(:,1).*scale_x fdata.WeightedCentroid(:,2)]) };
    
                    if 0
                        figure
                        subplot(1,2,1)
                        imagesc(rmap); hold on;
                        plot(fcents(:,1),fcents(:,2),'ko')
                        daspect([1 1 1])
                        axis xy
        
                        subplot(1,2,2)
                        imagesc(bmap); hold on;
                        plot(fcents(:,1),fcents(:,2),'ko')
                        daspect([1 1 1])
                        axis xy
                        keyboard
                    end
                end
            end
        end
    end

%% >>>>>>>>>> Behaviour anisotropy example trajectory
    xnow = 20;
    ynow = ynow-160;
    xsiz = 125;
    ysiz = 80;
    xbuff = 10;
    yplot_buff = 100;

    % plot arena
    ax1a = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('i',sprintf('Example behavioural anisotropy'),'yoffset',-20,'xoffset',20,'width',400);  

        % get data
        eg_session = 2;
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
        set(c,'LineWidth',1);

        % axis settings
        ax1a = gca;
        axis xy off tight
        daspect([1 1 1])
        cm = cmocean('balance');
        % cm = cmocean('thermal');        
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        colormap(ax1a,cmap)
        view(0,90)
        ax1a.XLim = [min(mf(:,1)) max(mf(:,1))];
        ax1a.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax1a.CLim = [-180 180];

        % text
        % text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        % text(0,1,sprintf('Example session'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1a.Position(1)+ax1a.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
        % get data
        pp = 2;
        idx = pos.pot > stimes(pp,1) & pos.pot < stimes(pp,2);
        pox = pos.pox_surficial(idx);
        poy = pos.poy_surficial(idx);
        poz = pos.poz_surficial(idx);
        poh = pos.yaw(idx);
        mf = posdata.maze_frame{eg_session};

        % plot data
        plot3(mf(:,4),mf(:,2),mf(:,3),'k'); hold on;        
        c = cline(pox,poy,poz,poh);
        set(c,'LineWidth',1);

        % axis settings
        axis xy off tight
        daspect([1 1 1])
        colormap(ax2,ax1a.Colormap)
        view(0,90)
        ax2.XLim = [min(mf(:,4)) max(mf(:,4))];
        ax2.YLim = [min(mf(:,2)) max(mf(:,2))];
        ax2.CLim = [-180 180];

        % additional plots
        ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

%% >>>>>>>>>> Field anisotropy example cell
    ynow = ynow-yplot_buff;

    uci = 'RG26_230119_t3_c5';
    m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
    m2 = clumaa.ratemap_surficial{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
    m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('k',sprintf('Example place field anisotropy'),'yoffset',-20,'xoffset',20,'width',400);  
    
        imagesc(m1,'alphadata',ones(size(m1)).*0.4); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,turbo)   
        ax1.CLim = [0 max([max(m1(:),[],'omitnan') 1])];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        % text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
        % text(0,1,sprintf('Example fields'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        
        % ellipses
        cm = cmocean('balance');
        % cm = cmocean('thermal');                
        cm = cm(32:224,:);
        cmap = [cm; flipud(cm);cm; flipud(cm)];
        cidx = linspace(-180,180,size(cmap,1));
        f = clumaa.planar_field_data{(ismember(clumaa.uci,uci) & clumaa.partn==1)};
        for ff = 1:size(f,1)
            e = ellipse(f.MajorAxisLength(ff,1)./2,f.MinorAxisLength(ff,1)./2,-deg2rad(f.Orientation(ff,1)),f.Centroid(ff,1),f.Centroid(ff,2),'Clipping','off');
            e.LineWidth = 2.5;
            c_idx = knnsearch(cidx(:),f.Orientation(ff,1));
            e.Color = cmap(c_idx,:);
        end

    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
        imagesc(m2,'alphadata',ones(size(m2)).*0.4); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,turbo)   
        ax2.CLim = [0 max([max(m2(:),[],'omitnan') 1])];
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(m2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

        % text
        % text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

        % ellipses
        f = clumaa.surficial_field_data{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
        for ff = 1:size(f,1)
            e = ellipse(f.MajorAxisLength(ff,1)./2,f.MinorAxisLength(ff,1)./2,-deg2rad(f.Orientation(ff,1)),f.Centroid(ff,1),f.Centroid(ff,2),'Clipping','off');
            e.LineWidth = 2.5;
            c_idx = knnsearch(cidx(:),f.Orientation(ff,1));
            e.Color = cmap(c_idx,:);
        end

%% >>>>>>>>>> Behaviour anisotropy map
    xnow = 330;
    ynow = ynow+yplot_buff;
    
    var = 'anisotropy_map';
    v1 = posdata.(var)(:,1); % arena 1 data
    all_amaps_a1 = cat(3,v1{:});
    a1 = mean(all_amaps_a1,3,'omitnan');

    v2 = posdata.(var)(:,2); % hills data (surficial)
    v2 = cellfun(@(x) imresize(x,[48 116],'bilinear'),v2,'UniformOutput',false);
    all_amaps_a2 = cat(3,v2{:});
    a2 = mean(all_amaps_a2,3,'omitnan');

    v3 = posdata.(var)(:,3); % arena 2 data
    all_amaps_a3 = cat(3,v3{:});
    a3 = mean(all_amaps_a3,3,'omitnan');

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
        ah = add_panel_title('j',sprintf('Avg. behavioural anisotropy'),'yoffset',-20,'xoffset',20,'width',400);  
    
        imagesc(a1,'alphadata',~isnan(a1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,cmocean('balance'))
        % colormap(ax1,cmocean('thermal'))        
        ax1.CLim = [-0.8 0.8];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,1,sprintf('N = %d sessions',size(all_amaps_a1,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')

    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
        imagesc(a2,'alphadata',~isnan(a2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)
        ax2.CLim = ax1.CLim;
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(a2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

        % text
        text(0,1,sprintf('N = %d sessions',size(all_amaps_a2,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+210 ynow+80 70 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Anisotropy score'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> Field anisotropy map
    ynow = ynow-yplot_buff;
    
    % get raw data
    var = 'field_anisotropy';
    v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
    d1 = cat(1,v1{:});

    v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data
    d2 = cat(1,v2{:});

    v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 2 data
    d3 = cat(1,v3{:});

    % map data
    % arena 1
    dist_cutoff = 160;
    mapset.binsize = 32; % (mm) firing rate map bin size
    
    dist_cutoff_bins = dist_cutoff ./ mapset.binsize;    
    anisotropy_map_1 = NaN(48,96);
    mapXY = d1(:,3:4);
    for bb = 1:numel(anisotropy_map_1)
        [r,c] = ind2sub(size(anisotropy_map_1),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_1(bb) = mean(d1(box_idx,1),1,'omitnan');
    end
    anisotropy_map_1 = imgaussfilt(anisotropy_map_1,1);

    % hills  
    anisotropy_map_2 = NaN(48,116);
    mapXY = d2(:,3:4);
    for bb = 1:numel(anisotropy_map_2)
        [r,c] = ind2sub(size(anisotropy_map_2),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_2(bb) = mean(d2(box_idx,1),1,'omitnan');
    end
    anisotropy_map_2 = imgaussfilt(anisotropy_map_2,1);

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
        ah = add_panel_title('l',sprintf('Avg. place field anisotropy'),'yoffset',-20,'xoffset',20,'width',400);  
    
        imagesc(anisotropy_map_1,'alphadata',~isnan(anisotropy_map_1)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax1,cmocean('balance'))
        % colormap(ax1,cmocean('thermal'))                
        ax1.CLim = [-0.5 0.5];
        ax1.XTick = [];
        ax1.YTick = [];

        % text
        text(0,1,sprintf('N = %d fields',size(d1,1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.2 ysiz]);
        imagesc(anisotropy_map_2,'alphadata',~isnan(anisotropy_map_2)); hold on;
        axis xy on
        daspect([1 1 1])
        colormap(ax2,ax1.Colormap)
        ax2.CLim = [-0.5 0.5];
        ax2.XTick = [];
        ax2.YTick = [];

        % additional plots
        ps = linspace(0,size(anisotropy_map_2,2),7);
        tops = ps(2:2:end);
        plot([tops; tops],ax2.YLim,'k:')

        % text
        text(0,1,sprintf('N = %d fields',size(d2,1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7,'Color','k')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+210 ynow+80 70 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Anisotropy score'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')






return
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig 3.png']; 
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

return







