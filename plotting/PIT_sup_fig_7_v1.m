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
                    % if pp==2
                        scale_x = 116 ./ size(bmap,2); % we need to scale field centroids to match a common map size
                    % end
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

%% >>>>>>>>>> Field anisotropy vs behaviour (case by case)
    var = 'field_anisotropy';
    v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
    v1 = cell2mat(v1);
    v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % arena 1 data
    v2 = cell2mat(v2);
    v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 1 data
    v3 = cell2mat(v3);

    xnow = 110;
    ynow = 500;
    bin_res = 256;

    % create axis
    ax1 = axes('Units','pixels','Position',[xnow ynow 160 160]);    
        ah = add_panel_title('a',sprintf(''),'yoffset',0,'xoffset',-70,'width',400);  
    
        msiz = 8;
        s1 = scatter(v1(:,1),v1(:,2),msiz,plot_set{1,1},'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none'); hold on;
        s2 = scatter(v2(:,1),v2(:,2),msiz,plot_set{1,2},'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none'); hold on;
    
        daspect([1 1 1])
        grid on
        box off
        axis square
        ax1 = gca;
        ax1.FontSize = 8;
        ax1.XLim = [-1 1];
        ax1.YLim = [-1 1];
        ax1.XTick = -1:0.25:1;
        ax1.YTick = -1:0.25:1;
        ax1.XTickLabel = {};
        ax1.YTickLabel = {};
    
        line([ax1.XLim(1) ax1.XLim(2)],[ax1.YLim(1) ax1.YLim(2)],'Color',[.5 .5 .5]);
    
        [r1,p1] = corr(v1(:,1),v1(:,2),'type','Pearson','rows','pairwise');
        [r2,p2] = corr(v2(:,1),v2(:,2),'type','Pearson','rows','pairwise');
    
        s1f = polyfit(v1(:,1),v1(:,2),1);
        rf1 = refline(s1f(1),s1f(2));
        set(rf1,'Color',plot_set{1,1},'LineWidth',2);
        s2f = polyfit(v2(:,1),v2(:,2),1);
        rf2 = refline(s2f(1),s2f(2));
        set(rf2,'Color',plot_set{1,2},'LineWidth',2);


        text(ax1,0.48,1.02,sprintf('{\\itr} = %.2f, {\\itp} < .001',r1),'Units','normalized','FontSize',7,'Color',plot_set{1,1},'VerticalAlignment','bottom','HorizontalAlignment','right')
        text(ax1,0.52,1.02,sprintf('{\\itr} = %.2f, {\\itp} < .001',r2),'Units','normalized','FontSize',7,'Color',plot_set{1,2},'VerticalAlignment','bottom','HorizontalAlignment','left')
    
        [~,leg] = legendflex([s1 s2],{'Arena 1','Hills'},'anchor',{'sw','sw'},'ncol',1,'box','off','buffer',[-80,-50],'xscale',0.5,'FontSize',8);  

    ax_b = axes('Units','pixels','Position',[ax1.Position(1)-50 ax1.Position(2) 45 ax1.Position(4)]); % behaviour y-axis plot
        xi = linspace(-1,1,bin_res);
        bs = 0.03;
        f1 = ksdensity(v1(:,2),xi,'Bandwidth',bs);
        f2 = ksdensity(v2(:,2),xi,'Bandwidth',bs);
    
        alph = 0.75;
        patch([f1 zeros(size(f1))],[xi fliplr(xi)],plot_set{1,1},'FaceAlpha',alph,'EdgeColor','none'); hold on;
        patch([f2 zeros(size(f2))],[xi fliplr(xi)],plot_set{1,2},'FaceAlpha',alph,'EdgeColor','none'); 

        ylabel(['Behavioural anisotropy'])
        ax_b.YLim = [-1 1];
        ax_b.XTick = [];
        ax_b.YTick = ax1.YTick;        
        ax_b.XColor = 'none';
        box off
        grid on
        ytickformat('%.2f')
        ax_b.FontSize = 8;
    
        % text(ax_b,-0.25,0,sprintf('Biased\nalong x-axis'),'Units','normalized','FontSize',8,'Color','k','HorizontalAlignment','right')
        % text(ax_b,-0.25,1,sprintf('Biased\nalong y-axis'),'Units','normalized','FontSize',8,'Color','k','HorizontalAlignment','right')
    
    ax_f = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-50 ax1.Position(3) 45]); % field, x-axis plot
        f1 = ksdensity(v1(:,1),xi,'Bandwidth',bs);
        f2 = ksdensity(v2(:,1),xi,'Bandwidth',bs);
    
        patch([xi fliplr(xi)],[f1 zeros(size(f1))],plot_set{1,1},'FaceAlpha',alph,'EdgeColor','none'); hold on;
        patch([xi fliplr(xi)],[f2 zeros(size(f2))],plot_set{1,2},'FaceAlpha',alph,'EdgeColor','none'); 

        xlabel(['Place field anisotropy'])
        ax_f.XLim = [-1 1];
        ax_f.YTick = [];
        ax_f.XTick = ax1.XTick;
        ax_f.YColor = 'none';        
        box off
        grid on
        xtickformat('%.2f')
        ax_f.FontSize = 8;

        % text(ax_f,0,-1,sprintf('Elongated\nalong x-axis'),'Units','normalized','FontSize',8,'Color','k','HorizontalAlignment','left')
        % text(ax_f,1,-1,sprintf('Elongated\nalong y-axis'),'Units','normalized','FontSize',8,'Color','k','HorizontalAlignment','right')

%% >>>>>>>>>> Anisotropy scatter density (arena 1)
    xnow = xnow+210;
    ybuff = 70;
    xsiz = 90;
    ysiz = 90;

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow+ybuff xsiz ysiz]);
        ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',0,'width',400);  
    
        % get distribution
        bin_res = 256;
        bw = 0.07;
        bs = linspace(-1,1,bin_res);
        [xx,yy] = ndgrid(bs,bs);
        f = mvksdensity(v1(:,[1 2]),[yy(:) xx(:)],'BandWidth',bw);
        F = reshape(f,size(xx));
    
        % plot data
        imagesc('XData',xx(:),'YData',yy(:),'CData',F);

        % axis settings
        axis xy
        daspect([1 1 1])
        ax1.XLim = [-1 1];
        ax1.YLim = [-1 1];
        ax1.XTick = [];
        
        % additional plots
        line([ax1.XLim(1) ax1.XLim(2)],[ax1.YLim(1) ax1.YLim(2)],'Color',[.5 .5 .5]);

        % text
        text(0.5,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
        text(-0.18,-0.1,sprintf('Behavioural anisotropy'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

    % plot hills
    ax1 = axes('Units','pixels','Position',[xnow ynow-40 xsiz ysiz]);
        % get distribution
        f = mvksdensity(v2(:,[1 2]),[yy(:) xx(:)],'BandWidth',bw);
        F = reshape(f,size(xx));
    
        % plot data
        imagesc('XData',xx(:),'YData',yy(:),'CData',F);

        % axis settings
        axis xy
        daspect([1 1 1])
        xlabel('Field anisotropy')
        ax = gca;
        ax.XLim = [-1 1];
        ax.YLim = [-1 1];

        % additional plots
        line([ax.XLim(1) ax.XLim(2)],[ax.YLim(1) ax.YLim(2)],'Color',[.5 .5 .5]);

        % text
        text(0.5,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S6.png']; 
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







