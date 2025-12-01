%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  analysis function written for:
% Grieves, Duvelle and Taube (202X) 
%
% USAGE:
%       [out] = template(in) process with default settings
% 
%       [out] = template(in,optional1) process using optional argument 1
% 
%       [out] = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'param1'          -   (default = X) Scalar value, parameter to do something
% 
%       'param2'          -   (default = X) Scalar value, parameter to do something
% 
% INPUT:
%       in    - input as a vector
%
% OUTPUT:
%       out   - output as a vector
%
% EXAMPLES:
%       % run function using default values
%       out = template(in,varargin)
%
% See also: GIT_audit

% HISTORY:
% version 1.0.0, Release 16/02/23 Code conception
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white
    fsiz = 9;
    fs = [15 10];

    sem = @(x,dim) std(x,0,dim,'omitmissing') ./ sqrt(sum(~isnan(x),dim));

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % collect real data
    if ~exist('model_maps','var') || 0
        ucis = unique(clumaa.uci(pidx));
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

%% >>>>>>>>>> Example cells
    xnow = 20;
    ynow = 780;
    xsiz = [90, 10]; % size, buffer
    ysiz = [70, 10]; % size, buffer
    xvec = xnow : sum(xsiz) : 1000;
    xvec = xvec(1:6);
    % yvec = ynow : -sum(ysiz).*1.45 : -1000;
    yvec = ynow : -sum(ysiz).*3 : -1000;
    
    yvec = yvec(1);
    [x,y] = ndgrid(xvec,yvec);
    ar_siz = [xsiz*0.48 ysiz*0.5];

    % find which cells to plot  
    % idx = [8 18 20 21 25 29];    
    idx = [31 20 21 25 29 40];    

    ucis = unique(clumaa.uci(pidx));
    m1c = clumaa.ratemap_planar(ismember(clumaa.uci,ucis) & clumaa.partn==1);
    m2c = clumaa.ratemap_planar(ismember(clumaa.uci,ucis) & clumaa.partn==2);

    for ii = 1:numel(x)
        %% plot arena (ratemap)
        m1 = model_maps{1,1}(:,:,idx(ii));
        m2 = model_maps{1,2}(:,:,idx(ii));

        ax3 = axes('Units','pixels','Position',[x(ii) y(ii) xsiz(1) ysiz(1)]);
            if ii==1
                ah = add_panel_title('A',sprintf(''),'yoffset',0,'xoffset',20,'width',400,'fontsize',fs);
            end
            imagesc(m1,'alphadata',~isnan(m1)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax3,turbo)   
            ax3.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
            ax3.XTick = [];
            ax3.YTick = [];
            if ii==1 || ii==7
                text(0,1.25,sprintf('BVC model'),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color','k','VerticalAlignment','bottom','rotation',0)                
                text(-0.02,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end
            text(1,1,sprintf('Cell %d, %.1f Hz',idx(ii),ax3.CLim(2)),'Units','normalized','HorizontalAlignment','right','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

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
            % text(0,-0.02,sprintf('%.1f',ax4.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','top')
            if ii==1 || ii==7
                text(-0.02,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

        % find best matching model cell
        rs = NaN(size(m1c,1),1);
        for jj = 1:size(m1c,1)
            a = [imresize(m1c{jj},[100,300],"nearest") imresize(m2c{jj},[100,300],"nearest")];
            b = [imresize(m1,[100,300],"nearest") imresize(m2,[100,300],"nearest")];
            rs(jj,1) = corr(a(:),b(:),'type','Pearson','rows','pairwise');
        end
        [~,mindx] = max(rs,[],'omitmissing');

        ax5 = axes('Units','pixels','Position',[x(ii) ax4.Position(2)-ax3.Position(4)+0 xsiz(1) ysiz(1)]);
            imagesc(m1c{mindx},'alphadata',~isnan(m1c{jj})); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax5,turbo) 
            ax5.XTick = [];
            ax5.YTick = [];
            ax5.CLim = [0 max([max(m1c{mindx}(:),[],'omitnan') max(m2c{mindx}(:),[],'omitnan') 1])];
            if ii==1 || ii==7
                text(0,1,sprintf('Real data (best match)'),'Units','normalized','HorizontalAlignment','left','FontSize',10,'Color','k','VerticalAlignment','bottom','rotation',0)
                text(-0.02,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

        ax6 = axes('Units','pixels','Position',[x(ii) ax5.Position(2)-ax3.Position(4)+20 xsiz(1) ysiz(1)]);
            imagesc(m2c{mindx},'alphadata',~isnan(m2c{jj})); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax6,turbo) 
            ax6.XTick = [];
            ax6.YTick = [];
            ax6.CLim = ax5.CLim;
            if ii==1 || ii==7
                text(-0.02,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

            % additional plots
            ps = linspace(0,size(m2c{mindx},2),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax6.YLim,'w:')

            if ii==1
                % colorbar
                axc = axes('Units','pixels','Position',[ax3.Position(1)+250 ax3.Position(2)+85 100 8]);
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

    end
% return

%% >>>>>>>>>> Repetition score
    xnow = 330;
    ynow = ynow-325;  

    % real data
    ax1r = axes('Units','pixels','Position',[xnow ynow 50 125]);
        ah = add_panel_title('C',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);               
    
        var = 'repetition_score';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 real data
        v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % hills real data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,1:3),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Repetition score'))    
        ax1r.XTick = 1:2;
        ax1r.XLim = [0.5 2.5]; 
        ax1r.FontSize = 8;
        ax1r.XTickLabel = maze_names(1:2);
        ax1r.YTick = [-1:0.5:2];
        ytickformat('%.1f');
        text(0,1.08,sprintf('Real data'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')

        % additional plots
        line(ax1r.XLim,[0 0],'Color',[.5 .5 .5])

        % stats
        axt = axes('Units','pixels','Position',ax1r.Position,'Color','none');
            axis off
            axt.XLim = ax1r.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);
            % N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];

    % model data
    ax2 = axes('Units','pixels','Position',[xnow+55 ynow 50 125]);
        v1 = rscores{1}(:); % arena 1 model data
        v2 = rscores{2}(:); % hills model data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,1:3),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ylabel(sprintf('Repetition score'))    
        ax2.XTick = ax1r.XTick;
        ax2.XLim = ax1r.XLim; 
        ax2.FontSize = 8;
        ax2.XTickLabel = ax1r.XTickLabel;
        ax2.YTick = ax1r.YTick;
        ax2.YColor = 'none';
        ax2.YLim = ax1r.YLim;
        ytickformat('%.1f');
        text(0,1.08,sprintf('BVC model'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')

        % additional plots
        line(ax2.XLim,[0 0],'Color',[.5 .5 .5])

        % stats
        axt = axes('Units','pixels','Position',ax2.Position,'Color','none');
            axis off
            axt.XLim = ax2.XLim;
            axt.YLim = [0 1];
            [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);
            % N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];
% keyboard

%% >>>>>>>>>> Remapping between arena and hills
    var = 'between_session_stability';
    v1 = clumaa.(var)(pidx & clumaa.partn==1,2); % arena 1 vs arena 2
    v2 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 vs hills
    v3 = NaN(size(model_maps{1,1},3),1);
    frate_cutoff = 1;
    for kk = 1:size(model_maps{1,1},3)
        m1 = model_maps{1,1}(:,:,kk);
        m2 = model_maps{1,2}(:,:,kk);            
        if max(m1(:),[],'omitnan')>=frate_cutoff || max(m2(:),[],'omitnan')>=frate_cutoff       
            v3(kk,1) = corr(m1(:),m2(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills (troughs only)
        end
    end

    % real data
    xnow = xnow+160;
    ynow = ynow;  
    ax1 = axes('Units','pixels','Position',[xnow ynow 50 125]);
        ah = add_panel_title('D',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);
            
        % vectorise
        [ds,gs] = vectorDATAGROUP([],v1(:),v2(:)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#012030'),hex2rgb('#45C4B0')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);
    
        % axis settings
        ylabel(sprintf('Correlation (\\it{r})'))
        ax1.XTick = 1:2;
        ax1.YTick = -1:0.5:1;  
        ax1.YLim = [-0.5 1];
        ax1.XLim = [0.5 2.5]; 
        ax1.XTickLabel = {};
        ax1.FontSize = 8;
        ytickformat('%.1f');
        ax1.XTickLabel = {sprintf('%s\\times%s',maze_names{4},maze_names{6}),sprintf('%s\\times%s',maze_names{4},maze_names{5})};
        text(0,1.08,sprintf('Real data'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')

        % text(0.5,-0.18,sprintf('Arena 1'),'FontSize',9,'HorizontalAl','center','VerticalAl','top','Units','normalized')
        % text(0.5,-0.09,sprintf('\\times'),'FontSize',9,'HorizontalAl','center','VerticalAl','top','Units','normalized')
        % text(0.25,0,sprintf('Arena 2'),'FontSize',8,'HorizontalAl','center','VerticalAl','top','Units','normalized')
        % text(0.75,0,sprintf('Hills'),'FontSize',8,'HorizontalAl','center','VerticalAl','top','Units','normalized')

        % % stats
        % axt = axes('Units','pixels','Position',ax1.Position,'Color','none');
        %     axis off
        %     axt.XLim = ax1.XLim;
        %     axt.YLim = [0 1];        
        %     [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.1,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);

    % model data
    ax2r = axes('Units','pixels','Position',[xnow+55 ynow 50 125]);
        % vectorise
        [ds,gs] = vectorDATAGROUP([],ones(size(v3(:))),v3(:)); % linearise data
    
        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',{hex2rgb('#012030'),hex2rgb('#45C4B0')},'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
        ax2r.XTick = 1:2;
        ax2r.YTick = [];  
        ax2r.YColor = 'none';
        ax2r.YLim = ax1.YLim;
        ax2r.XLim = ax1.XLim; 
        ax2r.XTickLabel = ax1.XTickLabel;
        ax2r.FontSize = 8;
        ax2r.XTickLabel = {sprintf('%s\\times%s',maze_names{4},maze_names{6}),sprintf('%s\\times%s',maze_names{4},maze_names{5})};
        text(0,1.08,sprintf('BVC model'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')

        % text(0.5,-0.18,sprintf('Arena 1'),'FontSize',9,'HorizontalAl','center','VerticalAl','top','Units','normalized')
        % text(0.5,-0.09,sprintf('\\times'),'FontSize',9,'HorizontalAl','center','VerticalAl','top','Units','normalized')
        % text(0.25,0,sprintf('Arena 2'),'FontSize',8,'HorizontalAl','center','VerticalAl','top','Units','normalized')
        % text(0.75,0,sprintf('Hills'),'FontSize',8,'HorizontalAl','center','VerticalAl','top','Units','normalized')

        % % stats
        % axt = axes('Units','pixels','Position',ax2.Position,'Color','none');
        %     axis off
        %     axt.XLim = ax2.XLim;
        %     axt.YLim = [0 1];        
        %     [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4,'bracket_y_base',axt.YLim(2)*1.1,'plot_omnibus',1,'omnibus_text_y_gap_coeff',2);
% keyboard
%% >>>>>>>>>> repetition (autocorrelation plots)
    xnow = 45;
    ynow = ynow-130;
    xbuff = 25;
    ysiz = 240;
    xsiz = 100;
    amap_cmap = cmocean('thermal');
    amap_clim = [-0.2 1];

    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('B',sprintf(''),'yoffset',-10,'xoffset',-5,'width',400,'fontsize',fs);                
    
        % get data
        v1 = clumaa.amap_cent(pidx & clumaa.partn==2); % arena data
        v1 = cell2mat(v1);
        v1idx = clumaa.repetition_score(pidx & clumaa.partn==2,1); % arena data
        [~,sidx] = sort(v1idx,'descend');
        v1s = v1(sidx,:);

        % plot data
        xd = ((1:size(v1s,2))-(size(v1s,2)./2)).*32./1000;
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
        text(0.5,1,sprintf('Real data'),'FontSize',8,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','bottom')

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
        v2 = model_maps{3,2};
        v2idx = rscores{2}(:);

        nindx = ~isnan(rscores{2}(:));
        v2 = v2(nindx,:);
        v2idx = v2idx(nindx);

        [~,sidx] = sort(v2idx,'descend');
        v2s = v2(sidx,:);

        % plot data
        xd = ((1:size(v2s,2))-(size(v2s,2)./2)).*10./1000;
        imagesc(v2s,'XData',xd); hold on;

        % axis settings
        axis ij on
        colormap(ax2,amap_cmap)   
        ax2.CLim = amap_clim;
        ax2.YTick = unique([1 0:200:size(v2s,1) size(v2s,1)]);
        ax2.XTick = [];
        ax2.FontSize = 8;

        % additional plots
        spacing_bins = 1;
        field_index = unique([0:spacing_bins:100 -(0:spacing_bins:100)]);
        plot([field_index; field_index],ax2.YLim,'w:')

        % text
        text(0.5,1,sprintf('BVC model'),'FontSize',8,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','bottom')

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
    axc = axes('Units','pixels','Position',[ax2.Position(1)-40 ax2.Position(2)+ax2.Position(4)+22 70 8]);
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

%% >>>>>>>>>> Field area vs distance from closest boundary
    xnow = ax1r.Position(1);
    ynow = ax1r.Position(2)-200;

    % real data, arena and hills
    ax1 = axes('Units','pixels','Position',[xnow ynow 50 125]);
        ah = add_panel_title('E',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);

        % arena
        d1 = datn{1,1}.distance_to_wall .* 32 .* 0.001; % convert to mm then m
        a1 = datn{1,1}.area .* 0.0001; % convert to m2

        % hills
        d2 = datn{1,2}.distance_to_wall .* 32 .* 0.001; % convert to m
        a2 = datn{1,2}.area .* 0.0001; % convert to m2

        % group based on distance
        cut_dist = 0.2;
        idx1 = d1 < cut_dist; % fields close to wall in arena
        idx2 = d2 < cut_dist; % fields close to wall in hills

        % main plot
        xi = [1 2];
        m1 = [mean(a1(idx1),'all','omitmissing') mean(a1(~idx1),'all','omitmissing')]; % arena: walls, center
        s1 = [sem(a1(idx1),'all') sem(a1(~idx1),'all')]; % arena: walls, center
        m2 = [mean(a2(idx2),'all','omitmissing') mean(a2(~idx2),'all','omitmissing')]; % hills: walls, center
        s2 = [sem(a2(idx2),'all') sem(a2(~idx2),'all')]; % hills: walls, center
        
        e1 = errorbar(xi,m1,s1,'Marker',plot_set{2,1},'Color',plot_set{1,1}); hold on
        e2 = errorbar(xi,m2,s2,'Marker',plot_set{2,2},'Color',plot_set{1,2}); hold on

        % axis settings
        ylabel('Area (m^{2})')
        ax1.XTick = 1:2;
        ax1.YTick = 0:0.1:1;
        ax1.XLim = [0.5 2.5]; 
        ax1.YLim = [0 0.35];
        ax1.XTickLabel = {'Wall','Centre'};
        ax1.FontSize = 8;
        text(0,1.02,sprintf('Real data'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')
        ytickformat('%.1f')
        box off

        text(0.1,1,sprintf(maze_names{1}),'HorizontalAlignment','left','VerticalAlignment','top','FontSize',9,'Color',plot_set{1,1},'Units','normalized')
        text(0.1,0.9,sprintf(maze_names{2}),'HorizontalAlignment','left','VerticalAlignment','top','FontSize',9,'Color',plot_set{1,2},'Units','normalized')

        % % stats
        % d = [d1(:); d2(:)];
        % g1 = [ones(size(d1)); ones(size(d2)).*2]; % maze group
        % g2 = [idx1(:); idx2(:)];
        % [P,T,STATS,TERMS] = anovan(d,[g1 g2],'display','on','varnames',{'maze','distance'},'model','full');

    % model data, arena and hills
    ax2 = axes('Units','pixels','Position',[xnow+55 ynow 50 125]);
        % arena
        d1 = mdatn{1,1}.distance_to_wall .* 0.01; % convert to m
        a1 = mdatn{1,1}.area .* 0.0001; % convert to m2

        % hills
        d2 = mdatn{1,2}.distance_to_wall .* 0.01; % convert to m
        a2 = mdatn{1,2}.area .* 0.0001; % convert to m2

        % group based on distance
        cut_dist = 0.2;
        idx1 = d1 < cut_dist; % fields close to wall in arena
        idx2 = d2 < cut_dist; % fields close to wall in hills

        % main plot
        xi = [1 2];
        m1 = [mean(a1(idx1),'all','omitmissing') mean(a1(~idx1),'all','omitmissing')]; % arena: walls, center
        s1 = [sem(a1(idx1),'all') sem(a1(~idx1),'all')]; % arena: walls, center
        m2 = [mean(a2(idx2),'all','omitmissing') mean(a2(~idx2),'all','omitmissing')]; % hills: walls, center
        s2 = [sem(a2(idx2),'all') sem(a2(~idx2),'all')]; % hills: walls, center
        
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
%% >>>>>>>>>> Elongation vs distance from closest boundary
    xnow = xnow+160;
    ynow = ynow;

    % real data, arena and hills
    ax1 = axes('Units','pixels','Position',[xnow ynow 50 125]);
        ah = add_panel_title('F',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);

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
        s1 = [sem(e1(idx1),'all') sem(e1(~idx1),'all')]; % arena: walls, center
        m2 = [mean(e2(idx2),'all','omitmissing') mean(e2(~idx2),'all','omitmissing')]; % hills: walls, center
        s2 = [sem(e2(idx2),'all') sem(e2(~idx2),'all')]; % hills: walls, center
        
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
        s1 = [sem(e1(idx1),'all') sem(e1(~idx1),'all')]; % arena: walls, center
        m2 = [mean(e2(idx2),'all','omitmissing') mean(e2(~idx2),'all','omitmissing')]; % hills: walls, center
        s2 = [sem(e2(idx2),'all') sem(e2(~idx2),'all')]; % hills: walls, center
        
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

%% >>>>>>>>>> Field-to-wall angle vs distance from closest boundary
    xnow = xnow-430;
    ynow = ynow-200;

    % real data, arena and hills
    ax1 = axes('Units','pixels','Position',[xnow ynow 70 125]);
        ah = add_panel_title('G',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);

        % arena
        d1 = datn{1,1}.distance_to_wall .* 32 .* 0.001; % convert to m
        e1 = datn{1,1}.field_wall_angle;

        % hills
        d2 = datn{1,2}.distance_to_wall .* 32 .* 0.001; % convert to m
        e2 = datn{1,2}.field_wall_angle;

        % group based on distance
        cut_dist = 0.2;
        idx1 = d1 < cut_dist; % fields close to wall in arena
        idx2 = d2 < cut_dist; % fields close to wall in hills

        % main plot
        xi = [1 2];
        m1 = [mean(e1(idx1),'all','omitmissing') mean(e1(~idx1),'all','omitmissing')]; % arena: walls, center
        s1 = [sem(e1(idx1),'all') sem(e1(~idx1),'all')]; % arena: walls, center
        m2 = [mean(e2(idx2),'all','omitmissing') mean(e2(~idx2),'all','omitmissing')]; % hills: walls, center
        s2 = [sem(e2(idx2),'all') sem(e2(~idx2),'all')]; % hills: walls, center
        
        p1 = errorbar(xi,m1,s1,'Marker',plot_set{2,1},'Color',plot_set{1,1}); hold on
        p2 = errorbar(xi,m2,s2,'Marker',plot_set{2,2},'Color',plot_set{1,2}); hold on

        % axis settings
        ylabel(sprintf('Field-to-wall angle (%c)',176))
        ax1.XTick = 1:2;
        ax1.YTick = 0:30:90;
        ax1.XLim = [0.5 2.5]; 
        ax1.YLim = [0 90];
        ax1.XTickLabel = {'Wall','Centre'};
        ax1.FontSize = 8;
        text(0,1.02,sprintf('Real data'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')
        ytickformat('%.f')
        box off

        text(0.1,1,sprintf(maze_names{1}),'HorizontalAlignment','left','VerticalAlignment','top','FontSize',9,'Color',plot_set{1,1},'Units','normalized')
        text(0.1,0.9,sprintf(maze_names{2}),'HorizontalAlignment','left','VerticalAlignment','top','FontSize',9,'Color',plot_set{1,2},'Units','normalized')

        % % stats
        % d = [d1(:); d2(:)];
        % g1 = [ones(size(d1)); ones(size(d2)).*2]; % maze group
        % g2 = [idx1(:); idx2(:)];
        % [P,T,STATS,TERMS] = anovan(d,[g1 g2],'display','on','varnames',{'maze','distance'},'model','full');
        % keyboard

    % model data, arena and hills
    ax2 = axes('Units','pixels','Position',[xnow+80 ynow 70 125]);
        % arena
        d1 = mdatn{1,1}.distance_to_wall .* 0.01; % convert to m
        e1 = mdatn{1,1}.field_wall_angle;

        % hills
        d2 = mdatn{1,2}.distance_to_wall .* 0.01; % convert to m
        e2 = mdatn{1,2}.field_wall_angle;

        % group based on distance
        cut_dist = 0.2;
        idx1 = d1 < cut_dist; % fields close to wall in arena
        idx2 = d2 < cut_dist; % fields close to wall in hills

        % main plot
        xi = [1 2];
        m1 = [mean(e1(idx1),'all','omitmissing') mean(e1(~idx1),'all','omitmissing')]; % arena: walls, center
        s1 = [sem(e1(idx1),'all') sem(e1(~idx1),'all')]; % arena: walls, center
        m2 = [mean(e2(idx2),'all','omitmissing') mean(e2(~idx2),'all','omitmissing')]; % hills: walls, center
        s2 = [sem(e2(idx2),'all') sem(e2(~idx2),'all')]; % hills: walls, center
        
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

        % % stats
        % d = [d1(:); d2(:)];
        % g1 = [ones(size(d1)); ones(size(d2)).*2]; % maze group
        % g2 = [idx1(:); idx2(:)];
        % [P,T,STATS,TERMS] = anovan(d,[g1 g2],'display','on','varnames',{'maze','distance'},'model','full');
        % keyboard

% %% >>>>>>>>>> Anisotropy map
%     xnow = xnow+205;
%     ynow = ynow;
% 
%     % Field density vs distance to boundary
%     ax1 = axes('Units','pixels','Position',[xnow ynow 140 125]);
% 
%         % map data
%         % arena 1
%         f = table2array(mdatn{2});
%         dist_cutoff = 160;
%         mapset.binsize = 10; % (mm) firing rate map bin size
% 
%         dist_cutoff_bins = dist_cutoff ./ mapset.binsize;    
%         anisotropy_map_1 = NaN(150,300);
%         mapXY = f(:,4:5);
%         for bb = 1:numel(anisotropy_map_1)
%             [r,c] = ind2sub(size(anisotropy_map_1),bb);
%             box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins);
%             anisotropy_map_1(bb) = mean(f(box_idx,10),1,'omitnan');
%         end
%         anisotropy_map_1 = imgaussfilt(anisotropy_map_1,1);
% 
%         figure
%         cmap_now = cmocean('balance');
% 
%         imagesc(anisotropy_map_1,'alphadata',~isnan(anisotropy_map_1)); hold on;
%         axis xy on
%         daspect([1 1 1])
%         ax1 = gca;
%         colormap(ax1,cmap_now)                
%         ax1.CLim = [-0.5 0.5];
%         ax1.XTick = [];
%         ax1.YTick = [];

%% >>>>>>>>>> Median distance from walls
    xnow = xnow+205;
    ynow = ynow;

    % Field density vs distance to boundary
    ax1 = axes('Units','pixels','Position',[xnow ynow 140 125]);
        ah = add_panel_title('H',sprintf(''),'yoffset',-10,'xoffset',-15,'width',400,'fontsize',fs);

        % % real data
        % % maze edges
        % maze_x = [0 0 3000 3000 0];
        % maze_y = [0 1500 1500 0 0];
        % 
        % % inner half area
        % maze_xh = maze_x.*(sqrt(2)/2);
        % maze_yh = maze_y.*(sqrt(2)/2);
        % maze_xh = maze_xh + (3000 - max(maze_xh))/2;
        % maze_yh = maze_yh + (1500 - max(maze_yh))/2;
        % 
        % % arena
        % arena_x = datn{1,1}.centroidx.*32-32;
        % arena_y = datn{1,1}.centroidy.*32-32;
        % hidx = inpolygon(arena_x,arena_y,maze_xh,maze_yh);
        % r_arena = [sum(hidx)./numel(hidx).*100 sum(~hidx)./numel(hidx).*100];
        % 
        % % hills
        % hills_x = datn{1,2}.centroidx.*32-32;
        % hills_y = datn{1,2}.centroidy.*32-32;
        % maze_x2 = maze_x .* 1.24;        
        % maze_xh2 = maze_xh .* 1.24;        
        % hidx = inpolygon(hills_x,hills_y,maze_xh2,maze_yh);
        % r_hills = [sum(hidx)./numel(hidx).*100 sum(~hidx)./numel(hidx).*100];
        % 
        % % model data
        % % arena
        % arena_x = mdatn{1,1}.centroidx.*10;
        % arena_y = mdatn{1,1}.centroidy.*10;
        % hidx = inpolygon(arena_x,arena_y,maze_xh,maze_yh);
        % m_arena = [sum(hidx)./numel(hidx).*100 sum(~hidx)./numel(hidx).*100];
        % 
        % % hills
        % hills_x = mdatn{1,2}.centroidx.*10;
        % hills_y = mdatn{1,2}.centroidy.*10;
        % hidx = inpolygon(hills_x,hills_y,maze_xh,maze_yh);
        % m_hills = [sum(hidx)./numel(hidx).*100 sum(~hidx)./numel(hidx).*100];

        % real data
        % arena
        d1 = datn{1,1}.distance_to_wall .* 32 ./ 1e3; % convert to m
        % hills
        d2 = datn{1,2}.distance_to_wall .* 32 ./ 1e3; % convert to m

        % group based on distance
        md1r = median(d1(:),'all');
        md2r = median(d2(:),'all');

        % main plot
        xi = [1 2];        
        p1 = plot(xi,[md1r md2r],'Marker','.','MarkerSize',25,'Color',hex2rgb('#0081A7')); hold on

        % model data
        % arena
        d1 = mdatn{1,1}.distance_to_wall .* 0.01; % convert to m
        % hills
        d2 = mdatn{1,2}.distance_to_wall .* 0.01; % convert to m

        % group based on distance
        md1m = median(d1(:),'all');
        md2m = median(d2(:),'all');

        % main plot
        p2 = plot(xi,[md1m md2m],'Marker','.','MarkerSize',25,'Color',hex2rgb('#F07167')); hold on

        % axis settings
        ylabel('Median field distance (m)')
        ax1.XTick = 1:2;
        ax1.YTick = 0:0.05:1;
        ax1.XLim = [0.5 2.5]; 
        ax1.YLim = [0.24 0.36];
        ax1.XTickLabel = {sprintf('%s',maze_names{1}),sprintf('%s',maze_names{2})};
        ax1.FontSize = 8;
        ytickformat('%.2f')
        box off

        % chance line
        iti = 100;
        nf = numel(d2);
        mr = NaN(iti,1);
        dmap = bwdist(padarray(zeros(150,300),[1 1],1,'both'),'euclidean');
        for rr = 1:iti
            xr = round( (299.*rand(nf,1))+2 );        
            yr = round( (149.*rand(nf,1))+2 );
            ds = dmap(sub2ind(size(dmap),yr(:),xr(:)));
            mr(rr,1) = median(ds(:),'all');
        end
        med_dist = mean(mr(:))./100;
        s1 = prctile(mr(:),99)./100;
        s2 = prctile(mr(:),1)./100;        
        l1 = line(ax1.XLim,[med_dist med_dist],'Color',[.5 .5 .5],'LineWidth',1.5);
        l2 = patch(ax1.XLim([1 1 2 2 1]),[s2 s1 s1 s2 s2],'k','FaceAlpha',0.2,'EdgeColor','none');
        uistack(l1,'bottom');
        uistack(l2,'bottom');

        text(0.1,1.0,sprintf('Real data'),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9,'Color',p1.Color,'Units','normalized')
        text(0.1,0.9,sprintf('BVC model'),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9,'Color',p2.Color,'Units','normalized')
        text(0.1,0.8,sprintf('Uniform'),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9,'Color','k','Units','normalized')

%% >>>>>>>>>> Firing rate vs distance from closest boundary
    xnow = xnow+180;
    ynow = ynow;

    % real data
    ax1 = axes('Units','pixels','Position',[xnow ynow 70 125]);
        ah = add_panel_title('I',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);

        clumaa.date = cell(size(clumaa,1),1);
        for kk = 1:size(clumaa,1)
            x = strsplit(clumaa.session_name{kk},'_');
            clumaa.date{kk} = x{1}(1:end-1);
        end
        % arena 1
        r = unique(clumaa.ratn);
        dmax = cell(size(r,1),2);
        pidx_max = false(size(pidx));        
        for rr = 1:length(r)
            d = unique(clumaa.date(ismember(clumaa.ratn,r(rr))));
            n = NaN(size(d,1),1);
            for dd = 1:length(d)
                n(dd,1) = sum( pidx & clumaa.partn==1 & ismember(clumaa.ratn,r(rr)) & ismember(clumaa.date,d(dd)) );
            end
            [~,maxidx] = max(n);
            dmax(rr,:) = {r(rr) d(maxidx)};
            pidx_max(pidx & ismember(clumaa.ratn,r(rr)) & ismember(clumaa.date,d(maxidx))) = true;
        end

        m1 = clumaa.ratemap_planar(pidx_max & clumaa.partn==1);        
        % m1 = clumaa.ratemap_planar(pidx & clumaa.partn==1);        
        m1 = cat(3,m1{:});

        dmap = zeros(size(m1(:,:,1)));
        dmap = padarray(dmap,[1 1],1,'both');
        dmap = bwdist(dmap);
        dmap = dmap(2:end-1,2:end-1);
        dmap = dmap .* 32 ./ 1000; % in m
              
        fun_mean = @(x) mean(x(:),"all",'omitmissing');
        fun_sem = @(x) std(x(:),[],"all",'omitmissing') ./ sqrt(sum(~isnan(x(:)),'all'));       

        % centre vs walls
        edg2 = [0,cut_dist,max(dmap(:))];   
        [~,~,bidx] = histcounts(dmap(:),edg2);  
        dat_a = NaN(length(edg2)-1,size(m1,3));            
        for jj = 1:size(m1,3)
            map_now = reshape(m1(:,:,jj),[],1);
            dat_a(:,jj) = accumarray(bidx(bidx>0),map_now(bidx>0),[],fun_mean);
        end

        % ridge
        m2 = clumaa.ratemap_planar(pidx_max & clumaa.partn==2);        
        % m2 = clumaa.ratemap_planar(pidx & clumaa.partn==2);
        m2 = cat(3,m2{:});

        dmap = zeros(size(m2(:,:,1)));
        dmap = padarray(dmap,[1 1],1,'both');
        dmap = bwdist(dmap);
        dmap = dmap(2:end-1,2:end-1);
        dmap = dmap .* 32 ./ 1000; % in m

        % centre vs walls
        [~,~,bidx] = histcounts(dmap(:),edg2);   
        dat_r = NaN(length(edg2)-1,size(m1,3));            
        for jj = 1:size(m2,3)
            map_now = reshape(m2(:,:,jj),[],1);
            dat_r(:,jj) = accumarray(bidx(bidx>0),map_now(bidx>0),[],fun_mean);
        end
        
        % plot data
        m = mean(dat_a,2,'omitmissing');
        e = sem(dat_a,2);
        errorbar(1:2,m,e,'Color',plot_set{1,1},'Marker',plot_set{2,1}); hold on;
        m = mean(dat_r,2,'omitmissing');
        e = sem(dat_r,2);        
        errorbar(1:2,m,e,'Color',plot_set{1,2},'Marker',plot_set{2,2});

        % axis settings
        ylabel(sprintf('Mean firing rate (Hz)'))
        ax1.YLim = [0 1.2];
        ax1.XTick = 1:2;
        ax1.YTick = ax1.YTick;        
        ax1.XLim = [0.5 2.5]; 
        ax1.XTickLabel = {'Wall','Centre'};
        ax1.FontSize = 8;
        ytickformat('%.1f');
        grid off
        box off
        text(0,1.02,sprintf('Real data'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')

        % % stats
        % d1 = dat_a(:);
        % d2 = dat_r(:);
        % d = [d1(:); d2(:)];
        % g1 = [ones(size(d1)); ones(size(d2)).*2]; % maze group
        % g2 = [reshape(ones(size(dat_a)).*[1;2],[],1); reshape(ones(size(dat_r)).*[1;2],[],1)];
        % [P,T,STATS,TERMS] = anovan(d,[g1 g2],'display','on','varnames',{'maze','distance'},'model','full');
        % keyboard

        text(0.1,0.2,sprintf(maze_names{1}),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9,'Color',plot_set{1,1},'Units','normalized')
        text(0.1,0.1,sprintf(maze_names{2}),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9,'Color',plot_set{1,2},'Units','normalized')

    % model data
    ax2 = axes('Units','pixels','Position',[xnow+80 ynow 70 125]);
        % arena 1
        m1 = model_maps{1,1};

        dmap = zeros(size(m1(:,:,1)));
        dmap = padarray(dmap,[1 1],1,'both');
        dmap = bwdist(dmap);
        dmap = dmap(2:end-1,2:end-1);
        dmap = dmap .* 10 ./ 1000; % in m
              
        fun_mean = @(x) mean(x(:),"all",'omitmissing');
        fun_sem = @(x) std(x(:),[],"all",'omitmissing') ./ sqrt(sum(~isnan(x(:)),'all'));       

        % centre vs walls
        edg2 = [0,cut_dist,max(dmap(:))];   
        [~,~,bidx] = histcounts(dmap(:),edg2);  
        dat_a = NaN(length(edg2)-1,size(m1,3));            
        for jj = 1:size(m1,3)
            map_now = reshape(m1(:,:,jj),[],1);
            dat_a(:,jj) = accumarray(bidx(bidx>0),map_now(bidx>0),[],fun_mean);
        end

        % ridge
        m1 = model_maps{1,2};

        dmap = zeros(size(m2(:,:,1)));
        dmap = padarray(dmap,[1 1],1,'both');
        dmap = bwdist(dmap);
        dmap = dmap(2:end-1,2:end-1);
        dmap = dmap .* 10 ./ 1000; % in m

        % centre vs walls
        [~,~,bidx] = histcounts(dmap(:),edg2);   
        dat_r = NaN(length(edg2)-1,size(m1,3));            
        for jj = 1:size(m2,3)
            map_now = reshape(m2(:,:,jj),[],1);
            dat_r(:,jj) = accumarray(bidx(bidx>0),map_now(bidx>0),[],fun_mean);
        end
        
        % plot data
        m = mean(dat_a,2,'omitmissing');
        e = sem(dat_a,2);
        errorbar(1:2,m,e,'Color',plot_set{1,1},'Marker',plot_set{2,1}); hold on;
        m = mean(dat_r,2,'omitmissing');
        e = sem(dat_r,2);        
        errorbar(1:2,m,e,'Color',plot_set{1,2},'Marker',plot_set{2,2},'LineStyle','-');

        % axis settings
        ylabel(sprintf('Mean firing rate (Hz)'))
        ax2.YLim = ax1.YLim;
        ax2.YColor = 'none';
        ax2.XTick = ax1.XTick;
        ax2.YTick = ax1.YTick;        
        ax2.XLim = ax1.XLim; 
        ax2.XTickLabel = ax1.XTickLabel;
        ax2.FontSize = 8;
        ytickformat('%.1f');
        grid off
        box off
        text(0,1.02,sprintf('BVC model'),'FontSize',8,'HorizontalAl','left','VerticalAl','bottom','Units','normalized')

        % % stats
        % d1 = dat_a(:);
        % d2 = dat_r(:);
        % d = [d1(:); d2(:)];
        % g1 = [ones(size(d1)); ones(size(d2)).*2]; % maze group
        % g2 = [reshape(ones(size(dat_a)).*[1;2],[],1); reshape(ones(size(dat_r)).*[1;2],[],1)];
        % [P,T,STATS,TERMS] = anovan(d,[g1 g2],'display','on','varnames',{'maze','distance'},'model','full');
        % keyboard

        % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig 5.png']; 
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
















