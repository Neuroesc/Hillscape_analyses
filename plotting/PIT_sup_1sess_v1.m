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
    animals = clumaa.rat;
    clumaa.date = cell(size(clumaa,1),1);
    for ii = 1:size(clumaa,1)
        clumaa.date{ii} = clumaa.session_name{ii}(1:6);
    end
    sessions = clumaa.date;
    isPlace = pidx & clumaa.partn==2;
    
    % Put into a table
    T = table(animals(:), sessions(:), isPlace(:), ...
        'VariableNames', {'Animal','Session','isPlace'});
    
    % Count number of place cells per animal/session
    G = groupsummary(T, {'Animal','Session'}, 'sum', 'isPlace');
    G.Properties.VariableNames{'sum_isPlace'} = 'PlaceCount';
    
    % For each animal, find the session with maximum place cells
    [~,idx] = unique(G.Animal); % indices of first occurrence of each animal
    result = table;
    for k = 1:numel(idx)
        thisAnimal = G.Animal(idx(k));
        rows = strcmp(G.Animal, thisAnimal);
        [~,maxIdx] = max(G.PlaceCount(rows));
        rowIdx = find(rows);
        result = [result; G(rowIdx(maxIdx),:)]; %#ok<AGROW>
    end
    
    % generate index
    pidx_ds = false(size(pidx));
    for rr = 1:size(result,1)
        pidx_ds(ismember(animals,result.Animal(rr)) & ismember(sessions,result.Session(rr)) & pidx) = true;
    end

%% >>>>>>>>>> Fields per m2
    xnow = 50;
    ynow = 700;
    xbuff = 60;

    % create axis
    ax1 = axes('Units','pixels','Position',[xnow ynow 90 125]);
        ah = add_panel_title('A',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        % var = 'planar_fields';
        % v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
        % v2a = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data
        % v2b = clumaa.surficial_fields(pidx & clumaa.partn==2,1); % hills data        
        % v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 2 data

        var = 'planar_fields';
        v1 = clumaa.(var)(pidx_ds & clumaa.partn==1 & clumaa.repetition_score(:,1)<0.28,1); % arena 1 data
        v2a = clumaa.(var)(pidx_ds & clumaa.partn==2 & clumaa.repetition_score(:,1)<0.28,1); % hills data
        v2b = clumaa.surficial_fields(pidx_ds & clumaa.partn==2 & clumaa.repetition_score(:,1)<0.28,1); % hills data        
        v3 = clumaa.(var)(pidx_ds & clumaa.partn==3 & clumaa.repetition_score(:,1)<0.28,1); % arena 2 data

        % convert to fields per m2
        v1 = v1 ./ (3*1.5);
        v2a = v2a ./ (3*1.5);
        v2b = v2b ./ ((3*1.2)*1.5);
        v3 = v3 ./ (3*1.5);

        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Place fields / m^{2}'))    
        ax = gca;
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
        text(5/6,-0.01,sprintf('%s      \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')
        % set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))
        ytickformat('%.1f');

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];        
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2));

%% >>>>>>>>>> field radius
    % create axis
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
        ah = add_panel_title('B',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'planar_field_info';
        v1 = clumaa.(var)(pidx_ds & clumaa.partn==1,1); % arena 1 data
        v2a = clumaa.(var)(pidx_ds & clumaa.partn==2,1); % hills data
        v2b = clumaa.surficial_field_info(pidx_ds & clumaa.partn==2,1); % hills data        
        v3 = clumaa.(var)(pidx_ds & clumaa.partn==3,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds.*10,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,[1 2 2 3]),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Field radius (mm)'))    
        ax = gca;
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        ax.FontSize = 8;
        ax.YScale = 'log';
        ax.YTick = [100 200 400 800];
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s      \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);
            N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];

        % % spatial info stats
        % var = 'planar_spatial_info_shuffles';
        % v1 = clumaa.(var)(pidx_ds & clumaa.partn==1,1); % arena 1 data
        % v2 = clumaa.(var)(pidx_ds & clumaa.partn==2,1); % hills data
        % [z,p,~,obs,~] = stats_perm_test(v1,v2,'iti',1000,'method','cliff','rep',1); % arena 1&hills vs shuffle
        % disp(sprintf('z = %.3f, p = %.1e, Cliff''s delta = %.2f, permutation z-test',z,p,obs))
        % [median(v1(:),"all",'omitmissing') median(v2(:),"all",'omitmissing')];

%% >>>>>>>>>> spatial information content
    xnow = xnow+155;

    % create axis
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+xbuff ax2.Position(2) ax2.Position(3) ax2.Position(4)]);
        ah = add_panel_title('C',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'planar_spatial_info_shuffles';
        v1 = clumaa.(var)(pidx_ds & clumaa.partn==1,2); % arena 1 data
        v2a = clumaa.(var)(pidx_ds & clumaa.partn==2,2); % hills data
        v2b = clumaa.surficial_spatial_info_shuffles(pidx_ds & clumaa.partn==2,2); % arena 1 data        
        [ds,gs] = vectorDATAGROUP([],v1(:),v2a(:),v2b(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dotcolor',plot_set(1,[1 2 2 3]),'dots',1,'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Spatial information (z)'))   
        ax = gca;
        ax.XTick = 1:4;
        ax.XLim = [0.5 3.5]; 
        ax.FontSize = 8;
        % ax.YScale = 'log';
        % ax.YTick = [100 200 400 800];
        ax.XTickLabel = {};
        text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(5/6,-0.01,sprintf('%s      \n(flattened)     ',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Rotation',25,'Units','normalized','VerticalAlignment','top')
        text(0.5,1.1,sprintf('n.s.'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize,'Units','normalized')

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'plot_omnibus',0);
            N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];


%% >>>>>>>>>> repetition score
    % create axis
    ax4 = axes('Units','pixels','Position',[ax3.Position(1)+ax3.Position(3)+xbuff ax3.Position(2) ax3.Position(3) ax3.Position(4)]);
        ah = add_panel_title('D',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'repetition_score';
        v1 = clumaa.(var)(pidx_ds & clumaa.partn==1,1); % arena 1 data
        v2 = clumaa.(var)(pidx_ds & clumaa.partn==2,1); % hills data
        v3 = clumaa.(var)(pidx_ds & clumaa.partn==3,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2(:),v3(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,1:3),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ax = gca;
        ylabel(sprintf('Repetition score'))    
        ax.XTick = 1:3;
        ax.XLim = [0.5 3.5]; 
        ax.FontSize = 8;
        ax.XTickLabel = maze_names(1:3);
        ax.YTick = [-1:0.5:2];
        ytickformat('%.1f');
        % text(1/6,-0.01,sprintf('%s',maze_names{1}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        % text(3/6,-0.01,sprintf('%s',maze_names{2}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')
        % text(5/6,-0.01,sprintf('%s',maze_names{3}),'FontSize',10,'HorizontalAlignment','right','Rotation',25,'Units','normalized','VerticalAlignment','top')

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])

        % stats
        axt = axes('Units','pixels','Position',ax.Position,'Color','none');
            axis off
            axt.XLim = ax.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);
            % N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];

%% >>>>>>>>>> Average field anisotropy
    % collect data
    ucis = unique(clumaa.uci(pidx_ds));
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

            % get field anisotropy
            fanisotropy = (fdata.BoundingBox(:,4) - fdata.BoundingBox(:,3)) ./  (fdata.BoundingBox(:,4) + fdata.BoundingBox(:,3)); % height - width / height + width

            % get elongation/eccentricity
            e = sqrt(1 - (fdata.MinorAxisLength(:,1) ./ fdata.MajorAxisLength(:,1)));

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
            dat = [dat; fdata.MajorAxisLength(:,1) fdata.MinorAxisLength(:,1) fdata.Orientation(:,1) fdata.WeightedCentroid(:,1).*scale_x fdata.WeightedCentroid(:,2) fdata.Area(:,1) ds(:,1) repmat(uu,size(fdata,1),1) ds(:,2) ds(:,3:4) fanisotropy e];
        end
        datn(pp) = { dat };
    end

    % arena 1
    v1 = datn{1}(:,3);
    v2 = datn{2}(:,3);
    v3 = datn{3}(:,3);

    xnow = -5;
    ynow = ynow-180;
    xbuff = 210;

    % map data
    % arena 1
    f = datn{1};
    dist_cutoff = 160;
    mapset.binsize = 32; % (mm) firing rate map bin size
    
    dist_cutoff_bins = dist_cutoff ./ mapset.binsize;    
    anisotropy_map_1 = NaN(48,96);
    mapXY = f(:,4:5);
    for bb = 1:numel(anisotropy_map_1)
        [r,c] = ind2sub(size(anisotropy_map_1),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_1(bb) = mean(f(box_idx,12),1,'omitnan');
    end
    anisotropy_map_1 = imgaussfilt(anisotropy_map_1,1);

    % hills  
    f = datn{2};    
    anisotropy_map_2 = NaN(48,116);
    mapXY = f(:,4:5);
    for bb = 1:numel(anisotropy_map_2)
        [r,c] = ind2sub(size(anisotropy_map_2),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_2(bb) = mean(f(box_idx,12),1,'omitnan');
    end
    anisotropy_map_2 = imgaussfilt(anisotropy_map_2,1);

    % arena 2  
    f = datn{3};    
    anisotropy_map_3 = NaN(48,96);
    mapXY = f(:,4:5);
    for bb = 1:numel(anisotropy_map_3)
        [r,c] = ind2sub(size(anisotropy_map_3),bb);
        box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
        anisotropy_map_3(bb) = mean(f(box_idx,12),1,'omitnan');
    end
    anisotropy_map_3 = imgaussfilt(anisotropy_map_3,1);

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow 220 90]);
        ah = add_panel_title('E',sprintf('Average field anisotropy'),'yoffset',-10,'xoffset',45,'width',400,'fontsize',fs);     
    
        PIT_plot_mazes(ax1,1,1,[],[],anisotropy_map_1,false);
                ve = [0 47];

        view(ve);
        lighting flat
        cmap_now = cmocean('balance');
        colormap(ax1,cmap_now)                
        ax1.CLim = [-0.5 0.5];
        ve = [0 47];
        view(ve);

        % text
        text(0,0,sprintf('%s, %d fields',maze_names{1},size(datn{1},1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
        PIT_plot_mazes(ax2,2,1,[],[],anisotropy_map_2,false)
        view(ve);
        lighting flat
        colormap(ax2,cmap_now)                
        ax2.CLim = [-0.5 0.5];
        view(ve);

        % text
        text(0,0,sprintf('%s, %d fields',maze_names{2},size(datn{2},1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')

    % plot arena 2
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+xbuff ax2.Position(2) ax1.Position(3) ax1.Position(4)]);
        PIT_plot_mazes(ax3,1,1,[],[],anisotropy_map_3,false)
        view(ve);
        lighting flat
        colormap(ax3,cmap_now)                
        ax3.CLim = [-0.5 0.5];
        view(ve);

        % text
        text(0,0,sprintf('%s, %d fields',maze_names{3},size(datn{3},1)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+400 ynow+100 100 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 9;

        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(-0.3,1.5,sprintf('Anisotropy score:'),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')

    % example ellipses
    ax2e = axes('Units','pixels','Position',[axc.Position(1) axc.Position(2)+20 axc.Position(3) axc.Position(3)],'Color','none');
        ax2e.XLim = axc.XLim;
        ax2e.YLim = axc.XLim;
        daspect([1 1 1])
        axis manual off
        ax2e.Clipping = 'off';

        nf = 5;
        xvals = linspace(-0.5,0.5,nf); 
        ye = -0.5;
        alph = 0.5;
        scaler = 0.22;
        c = cmocean('balance',nf);
        for xx = 1:length(xvals)
            emin = 1;
            emax = -(emin + emin*xvals(xx))/(xvals(xx) - 1);
            v = sqrt(emin^2+emax^2);
            emin = emin/v*scaler;
            emax = emax/v*scaler;
            r = rectangle('Position',[xvals(xx)-(emin/2),ye(1)-(emax/2),emin,emax],'Curvature',[1 1],'EdgeColor','k','Clipping','off','FaceColor',[c(xx,:) 0.5]); hold on;
        end
        colormap(ax2e,c);

%% >>>>>>>>>> Behaviour anisotropy map 
    ynow = ynow-140;

    session_subset = unique(clumaa.pos_idx(pidx_ds));

    var = 'anisotropy_map';
    v1 = posdata.(var)(:,1); % arena 1 data
    all_amaps_a1 = cat(3,v1{session_subset});
    a1 = mean(all_amaps_a1,3,'omitnan');

    v2 = posdata.(var)(:,2); % hills data (surficial)
    v2 = cellfun(@(x) imresize(x,[48 116],'bilinear'),v2,'UniformOutput',false);
    all_amaps_a2 = cat(3,v2{session_subset});
    a2 = mean(all_amaps_a2,3,'omitnan');

    v3 = posdata.(var)(:,3); % arena 2 data
    all_amaps_a3 = cat(3,v3{session_subset});
    a3 = mean(all_amaps_a3,3,'omitnan');

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow 220 90]);
        ah = add_panel_title('F',sprintf('Average behavioural anisotropy'),'yoffset',-10,'xoffset',45,'width',400,'fontsize',fs);     
    
        PIT_plot_mazes(ax1,1,1,[],[],a1,false);
        view(ve);
        lighting flat
        colormap(ax1,cmap_now)                
        ax1.CLim = [-0.8 0.8];
        ve = [0 47];
        view(ve);

        % text
        text(0,0,sprintf('%s, %d sessions',maze_names{1},size(all_amaps_a1,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')
        
    % plot hills
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
        PIT_plot_mazes(ax2,2,1,[],[],a2,false)
        view(ve);
        lighting flat
        colormap(ax2,cmap_now)                
        ax2.CLim = [-0.8 0.8];
        view(ve);

        % text
        text(0,0,sprintf('%s, %d sessions',maze_names{2},size(all_amaps_a2,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')

    % plot arena 2
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+xbuff ax2.Position(2) ax1.Position(3) ax1.Position(4)]);
        PIT_plot_mazes(ax3,1,1,[],[],a3,false)
        view(ve);
        lighting flat
        colormap(ax3,cmap_now)                
        ax3.CLim = [-0.8 0.8];
        view(ve);

        % text
        text(0,0,sprintf('%s, %d sessions',maze_names{1},size(all_amaps_a1,3)),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',8,'Color','k')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[xnow+400 ynow+100 100 8]);
        x = linspace(ax1.CLim(1),ax1.CLim(2),100);
        imagesc(x,'XData',ax1.CLim,'YData',[0 1]);
        colormap(axc,ax1.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 9;
        text(0,0.5,sprintf('%.1f ',ax1.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax1.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(-0.3,1.5,sprintf('Anisotropy score:'),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> repetition (autocorrelation plots)
    xnow = 45;
    ynow = ynow-170;
    xbuff = 25;
    ysiz = 100;
    xsiz = 100;
    amap_cmap = cmocean('thermal');
    amap_clim = [-0.2 1];

    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('G',sprintf(''),'yoffset',-20,'xoffset',-5,'width',400,'fontsize',fs);                
    
        % get data
        v1 = clumaa.amap_cent(pidx_ds & clumaa.partn==1); % arena data
        v1 = cell2mat(v1);
        v1idx = clumaa.repetition_score(pidx_ds & clumaa.partn==1,1); % arena data
        [~,sidx] = sort(v1idx,'descend');
        v1s = v1(sidx,:);

        % plot data
        xd = ((1:size(v1s,2))-(size(v1s,2)./2)).*mapset.binsize./1000;
        imagesc(v1s,'XData',xd); hold on;

        % axis settings
        axis ij on
        colormap(ax1,amap_cmap)   
        ax1.CLim = amap_clim;
        ax1.YTick = unique([1 0:50:size(v1s,1) size(v1s,1)]);
        ax1.XTick = [];
        ylabel('Ranked cells')
        ax1.FontSize = 8;

        % additional plots
        spacing_bins = 1;
        field_index = unique([0:spacing_bins:100 -(0:spacing_bins:100)]);        
        plot([field_index; field_index],ax1.YLim,'w:')

        % text
        text(0.5,1,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','bottom')

    ax1b = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-65 xsiz 50]);
        m = mean(v1s(:,2:end-1),1,'omitmissing');
        e = std(v1s(:,2:end-1),[],1,'omitmissing');

        % plot data
        b = boundedline(xd(2:end-1),m,e,'cmap',plot_set{1,1});

        % axis settings
        xlabel('Lag (m)')
        ax1b.XLim = ax1.XLim;
        ax1b.XTick = -3:1:3;
        ytickformat('%.1f')
        ylabel('Avg. r')
        ax1b.YLim = ax1.CLim;
        % line(ax1b.XLim,[0 0],'Color',[.5 .5 .5])
        set(gca,'Layer','top')

    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
        % get data
        v2 = clumaa.amap_cent(pidx_ds & clumaa.partn==2); % hills data
        v2 = cell2mat(v2);
        v2idx = clumaa.repetition_score(pidx_ds & clumaa.partn==2,1); % hills data
        [~,sidx] = sort(v2idx,'descend');
        v2s = v2(sidx,:);

        % plot data
        xd = ((1:size(v2s,2))-(size(v2s,2)./2)).*mapset.binsize./1000;
        imagesc(v2s,'XData',xd); hold on;

        % axis settings
        axis ij on
        colormap(ax2,amap_cmap)   
        ax2.CLim = amap_clim;
        ax2.YTick = unique([1 0:50:size(v2s,1) size(v2s,1)]);
        ax2.XTick = [];
        ax2.FontSize = 8;

        % additional plots
        spacing_bins = 1;
        field_index = unique([0:spacing_bins:100 -(0:spacing_bins:100)]);
        plot([field_index; field_index],ax2.YLim,'w:')

        % text
        text(0.5,1,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','bottom')

    ax2b = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)-65 xsiz 50]);
        m = mean(v2s(:,2:end-1),1,'omitmissing');
        e = std(v2s(:,2:end-1),[],1,'omitmissing');

        % plot data
        b = boundedline(xd(2:end-1),m,e,'cmap',plot_set{1,2});

        % axis settings
        xlabel('Lag (m)')
        ax2b.XLim = ax2.XLim;        
        ax2b.XTick = -3:1:3;
        ytickformat('%.1f')
        ax2b.YLim = ax2.CLim;
        % line(ax2b.XLim,[0 0],'Color',[.5 .5 .5])
        set(gca,'Layer','top')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax2.Position(1)-40 ax2.Position(2)+ax2.Position(4)+18 70 8]);
        x = linspace(ax2.CLim(1),ax2.CLim(2),100);
        imagesc(x,'XData',ax2.CLim,'YData',[0 1]);
        colormap(axc,ax2.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax2.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax2.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Correlation'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> Elongation vs distance from closest boundary
    % collect real data
    if 1
        ucis = unique(clumaa.uci(pidx_ds));
        datn = cell(1,3);
        for pp=1:3    
            dat = table;
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
    
                % get field anisotropy
                fanisotropy = (fdata.BoundingBox(:,4) - fdata.BoundingBox(:,3)) ./  (fdata.BoundingBox(:,4) + fdata.BoundingBox(:,3)); % height - width / height + width

                % get elongation/eccentricity
                e = sqrt(1 - (fdata.MinorAxisLength(:,1) ./ fdata.MajorAxisLength(:,1)));

                % prepare environment polygon
                b = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0]+0.5; 
                a = [90 0 90 0];
                gpoly = [];
                res = 100;
                for jj = 1:size(b,1)-1
                    gpoly = [gpoly; linspace(b(jj,1),b(jj+1,1),res)' linspace(b(jj,2),b(jj+1,2),res)' repmat(a(jj),res,1)];
                end
    
                % get distance to closest wall
                ds = NaN(size(fdata,1),4);
                for jj = 1:size(fdata,1)
                    [i,d] = knnsearch(gpoly(:,[1:2]),[fdata.Centroid(jj,1),fdata.Centroid(jj,2)],'K',1);
                    ds(jj,1:2) = [d gpoly(i,3)];             
                end

                % get field to wall angle
                fun = @(x) abs( rad2deg( atan( (tan(deg2rad(x(:,1)))-tan(deg2rad(x(:,2)))) ./ (1+tan(deg2rad(x(:,2))).*tan(deg2rad(x(:,1)))) ) ) );
                field_wall_angle = fun( [ds(:,2) fdata.Orientation(:,1)] );

                % accumulate
                dat_now = table(fdata.MajorAxisLength(:,1),fdata.MinorAxisLength(:,1),fdata.Orientation(:,1),fdata.WeightedCentroid(:,1).*scale_x,fdata.WeightedCentroid(:,2),fdata.Area(:,1),ds(:,1),repmat(uu,size(fdata,1),1),ds(:,2),fanisotropy,e,field_wall_angle,...
                    'VariableNames',["majoraxis","minoraxis","orientation","centroidx","centroidy","area","distance_to_wall","cell","wall_angle","anisotropy","elongation","field_wall_angle"]);
                dat = [dat; dat_now];
            end
            datn(pp) = { dat };
        end

        % collect model data
        enames = {'3000x1500_arena','3000x1500_hills'};
        scratch_space = 'C:\Users\roddyg\OneDrive - University of Glasgow\Projects in prep\2021 Place cells irregular terrain\associated data\bvc_model';        
        mname = [scratch_space '\3000x1500_arena_maps'];
        load(mname,'env','amap','ename','dmap','-mat');

        zcut = 1.2; % z-score for field detection in z-scored firing rate maps   
        frcut = 1; % place fields must have a peak firing rate at least greater than this value
        arcut = 400; % cm2, place fields must have a total area greater than this value 

        mdatn = cell(1,length(enames));
        model_maps = cell(3,length(enames));
        rscores = cell(1,2); % modelled repetition scores
        for ee = 1:length(enames)
            hname = [scratch_space '\20240221T093430\' enames{ee} '_512b_1028h_hpcs'];   
            load(hname,'hpc_maps','-mat');
            model_maps{1,ee} = hpc_maps;

            out_name = [scratch_space '\20240221T093430\' enames{ee} '_512b_1028h_hpcs_analysis'];  
            if exist(out_name,'file')
                load(out_name,'amaps','amap_cents','rnow','mdat','-mat');    
                mdatn(ee) = { mdat }; % field data
                model_maps{2,ee} = amaps; % autocorrelations
                model_maps{3,ee} = amap_cents; % autocorrelation central rows             
                rscores{ee} = rnow; % repetition scores
            else
                % fields
                nhpcs = size(hpc_maps,3);
                mdat = table;
                for hh = 1:nhpcs
                    map = hpc_maps(:,:,hh);
                    zmap = ( map - mean(map(:),'omitnan') ) / std(map(:),'omitnan'); % zscore ratemap
                    thresh_ratemap = imbinarize(zmap,zcut); % 2 s.d. threshold
        
                    datout = regionprops('table',thresh_ratemap,map,'Area','Centroid','WeightedCentroid','MajorAxisLength','MinorAxisLength','Orientation','MaxIntensity','BoundingBox'); % detect contiguous regions
                    if isempty(datout) % filter out small or low firing regions
                        continue
                    end                     
                    datout(datout.Area(:) < arcut | datout.MaxIntensity(:) < frcut,:) = [];      
        
                    % get field to wall distance
                    mat_cents = round(datout.Centroid(:,[2 1])); % row, column
                    i = sub2ind(size(dmap),mat_cents(:,1),mat_cents(:,2));            
                    datout.field_wall_distance = dmap(i);
        
                    % get field anisotropy
                    datout.anisotropy = (datout.BoundingBox(:,4) - datout.BoundingBox(:,3)) ./  (datout.BoundingBox(:,4) + datout.BoundingBox(:,3)); % height - width / height + width
        
                    % get elongation/eccentricity
                    datout.e = sqrt(1 - (datout.MinorAxisLength(:,1) ./ datout.MajorAxisLength(:,1)));
        
                    % get angle to wall
                    closest_wall_angle = amap(i);     
                    fun = @(x) abs( rad2deg( atan( (tan(deg2rad(x(:,1)))-tan(deg2rad(x(:,2)))) ./ (1+tan(deg2rad(x(:,2))).*tan(deg2rad(x(:,1)))) ) ) );
                    datout.field_wall_angle = fun( [closest_wall_angle(:) datout.Orientation(:)] );
        
                    % accumulate
                    mdat_now = table(datout.MajorAxisLength(:,1),datout.MinorAxisLength(:,1),datout.Orientation(:,1),datout.WeightedCentroid(:,1),datout.WeightedCentroid(:,2),datout.Area(:,1),datout.field_wall_distance,repmat(hh,size(datout,1),1),closest_wall_angle,datout.anisotropy,datout.e,datout.field_wall_angle,...
                        'VariableNames',["majoraxis","minoraxis","orientation","centroidx","centroidy","area","distance_to_wall","cell","wall_angle","anisotropy","elongation","field_wall_angle"]);
                    mdat = [mdat; mdat_now];
                end
                mdatn(ee) = { mdat };
                
                % autocorrelations
                amaps = [];
                amap_cents = [];
                rnow = NaN(nhpcs,1);
                for hh = 1:nhpcs
                    map = hpc_maps(:,:,hh);
                    amap = ndautoCORR(map,map,50); % spatial autocorrelogram
                    if isempty(amaps)
                        amaps = NaN([size(amap) nhpcs]);
                        amap_cents = NaN([nhpcs size(amap,2)]);
                    end
                    amaps(:,:,hh) = amap;
    
                    % repetition score
                    spacing_bins = 100; % x-axis bins per hill
                    cutoff_dist = 12; % cm, area around bin locations to combine                
                    field_index = ceil(sort(unique([size(amap,2)/2:spacing_bins:size(amap,2)-spacing_bins/2,size(amap,2)/2:-spacing_bins:spacing_bins/2])));
                    field_index(3) = []; % remove central field
                    non_field_index = ceil(sort(unique([size(amap,2)/2+spacing_bins/2:spacing_bins:size(amap,2)-spacing_bins/2,size(amap,2)/2-spacing_bins/2:-spacing_bins:spacing_bins/2])));
                                
                    amap_idx = zeros(1,size(amap,2));
                    amap_idx(field_index) = 1;
                    dist_idx = bwdist(amap_idx,'euclidean');
                    dist_log_fields = dist_idx < cutoff_dist;
                    
                    amap_idx = zeros(1,size(amap,2));
                    amap_idx(non_field_index) = 1;
                    dist_idx = bwdist(amap_idx,'euclidean');
                    dist_log_nfields = dist_idx < cutoff_dist;
    
                    central_row = amap( ceil(size(amap,1)/2),: );
                    amap_cents(hh,:) = central_row;
                    rnow(hh,1) = mean(central_row(dist_log_fields),'omitnan') - mean(central_row(dist_log_nfields),'omitnan');                
                end
                model_maps{2,ee} = amaps;
                model_maps{3,ee} = amap_cents;
                rscores{ee} = rnow;
    
                save(out_name,'amaps','amap_cents','rnow','mdat','-mat');    
            end       
        end
    end

    xnow = xnow+280;
    ynow = ynow-40;

    % real data, arena and hills
    ax1 = axes('Units','pixels','Position',[xnow ynow 50 125]);
        ah = add_panel_title('H',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);

        % arena
        d1 = datn{1,1}.distance_to_wall .* 32 .* 0.001; % convert to m
        e1 = datn{1,1}.elongation;

        % hills
        d2 = datn{1,2}.distance_to_wall .* 32 .* 0.001; % convert to m
        e2 = datn{1,2}.elongation;

        % group based on distance
        cut_dist = 0.2;
        idx1 = d1 < cut_dist; % fields close to wall in arena
        idx2 = d2 < cut_dist; % fields close to wall in hills

        % main plot
        xi = [1 2];
        m1 = [mean(e1(idx1),'all','omitmissing') mean(e1(~idx1),'all','omitmissing')]; % arena: walls, center
        s1 = [nansem(e1(idx1)) nansem(e1(~idx1))]; % arena: walls, center
        m2 = [mean(e2(idx2),'all','omitmissing') mean(e2(~idx2),'all','omitmissing')]; % hills: walls, center
        s2 = [nansem(e2(idx2)) nansem(e2(~idx2))]; % hills: walls, center
        
        p1 = errorbar(xi,m1,s1,'Marker',plot_set{2,1},'Color',plot_set{1,1}); hold on
        p2 = errorbar(xi,m2,s2,'Marker',plot_set{2,2},'Color',plot_set{1,2}); hold on

        % axis settings
        ylabel('Elongation')
        ax1.XTick = 1:2;
        ax1.YTick = 0:0.1:1;
        ax1.XLim = [0.5 2.5]; 
        ax1.YLim = [0.5 0.8];
        ax1.XTickLabel = {'Wall','Centre'};
        ax1.FontSize = 8;
        text(0,1.02,sprintf('Real data'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')
        ytickformat('%.1f')
        box off

%         % stats
%         d = [d1(:); d2(:)];
%         g1 = [ones(size(d1)); ones(size(d2)).*2]; % maze group
%         g2 = [idx1(:); idx2(:)]+1;
%         [P,T,STATS,TERMS] = anovan(d,[g1 g2],'display','on','varnames',{'maze','distance'},'model','linear');
% keyboard

    % model data, arena and hills
    ax2 = axes('Units','pixels','Position',[xnow+55 ynow 50 125]);
        % arena
        d1 = mdatn{1,1}.distance_to_wall .* 0.01; % convert to m
        e1 = mdatn{1,1}.elongation;

        % hills
        d2 = mdatn{1,2}.distance_to_wall .* 0.01; % convert to m
        e2 = mdatn{1,2}.elongation;

        % group based on distance
        cut_dist = 0.2;
        idx1 = d1 < cut_dist; % fields close to wall in arena
        idx2 = d2 < cut_dist; % fields close to wall in hills

        % main plot
        xi = [1 2];
        m1 = [mean(e1(idx1),'all','omitmissing') mean(e1(~idx1),'all','omitmissing')]; % arena: walls, center
        s1 = [nansem(e1(idx1)) nansem(e1(~idx1))]; % arena: walls, center
        m2 = [mean(e2(idx2),'all','omitmissing') mean(e2(~idx2),'all','omitmissing')]; % hills: walls, center
        s2 = [nansem(e2(idx2)) nansem(e2(~idx2))]; % hills: walls, center
        
        p1 = errorbar(xi,m1,s1,'Marker',plot_set{2,1},'Color',plot_set{1,1}); hold on
        p2 = errorbar(xi,m2,s2,'Marker',plot_set{2,2},'Color',plot_set{1,2}); hold on

        % axis settings
        ax2.XTick = ax1.XTick;
        ax2.YTick = [];
        ax2.YColor = 'none';
        ax2.XLim = ax1.XLim; 
        ax2.YLim = ax1.YLim;
        ax2.XTickLabel = ax1.XTickLabel;
        ax2.FontSize = ax1.FontSize;
        text(0,1.02,sprintf('BVC model'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')
        box off

%         % stats
%         d = [d1(:); d2(:)];
%         g1 = [ones(size(d1)); ones(size(d2)).*2]; % maze group
%         g2 = [idx1(:); idx2(:)];
%         [P,T,STATS,TERMS] = anovan(d,[g1 g2],'display','on','varnames',{'maze','distance'},'model','full');
% keyboard

%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S14.png']; 
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












