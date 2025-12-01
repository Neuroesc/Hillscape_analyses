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
    %% Map settings - must match hill_analysis_v2
    mapset.ppm          = 1000;
    mapset.method       = 'histogram';
    mapset.binsize      = 32; % (mm) firing rate map bin size
    mapset.ssigma       = 64; % (mm) firing rate map smoothing sigma
    mapset.padding      = 10; % (mm) how much to pad the edges of the firing rate map
    mapset.mindwell     = 0.01; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    mapset.mindist      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.smethod      = 1; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing    
    mapset.frcut        = 0.2; % cutoff % for regions to be considered a place field
    mapset.arcut        = 36; % cutoff area for regions to be considered a place field
    mapset.drive_height_mm = 20;
    mapset.srate = 50;

    % settings for anisotropy calculation
    dist_cutoff = 120;
    dist_cutoff_bins = dist_cutoff ./ mapset.binsize;
    min_run_speed = 5; % cm/s



    % % Create figure
    % fig_now = figure('Units','pixels','Position',[100 100 1000 900],'visible','on');
    % set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    % set(gcf,'color','w'); % makes the background colour white

    if ~any(ismember(clumaa.Properties.VariableNames,'field_anisotropy')) % if the column(s) do not exist yet
        clumaa.field_anisotropy = cell(size(clumaa,1),1); % preallocate
    end 

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Field anisotropy
if 0
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
                clumaa.field_anisotropy(idx) = { single([fanisotropy banisotropy repmat(clumaa.repetition_score(idx,1),size(fanisotropy,1),1)]) };

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

%% >>>>>>>>>> Field anisotropy vs behaviour (overall)
    % clumaa.planar_field_aspect_ratio = clumaa.planar_field_info(:,4) ./ clumaa.planar_field_info(:,5); % aspect ratio = field height / width
    clumaa.planar_field_aspect_ratio = (clumaa.planar_field_info(:,4) - clumaa.planar_field_info(:,5)) ./ (clumaa.planar_field_info(:,4) + clumaa.planar_field_info(:,5)); % aspect ratio = field height / width
    
    % pp = 1;
    % cidx1 = pidx & clumaa.partn==pp & clumaa.f_rate_hz>0.1;
    % x_field_anisotropy_1 = clumaa.planar_field_aspect_ratio(cidx1);
    % [~,bidx1] = ismember(clumaa.pos_idx(cidx1),posdata.pos_idx);
    % y_behav_anisotropy_1 = posdata.anisotropy_score(bidx1,pp);
    % 
    % pp = 2;
    % cidx1 = pidx & clumaa.partn==pp & clumaa.f_rate_hz>0.1;
    % x_field_anisotropy_2 = clumaa.planar_field_aspect_ratio(cidx1);
    % [~,bidx1] = ismember(clumaa.pos_idx(cidx1),posdata.pos_idx);
    % y_behav_anisotropy_2 = posdata.anisotropy_score(bidx1,pp);
    % 
    % v1 = [x_field_anisotropy_1(:) y_behav_anisotropy_1(:)];
    % v2 = [x_field_anisotropy_2(:) y_behav_anisotropy_2(:)];
    % v3 = [v1;v2];
    % 
    % figure('Units','pixels','Position',[100 100 1000 900],'visible','on');
    % ax = axes('Units','pixels','Position',[250 350 400 400]);
    % 
    %     s1 = scatter(v1(:,1),v1(:,2),50,plot_set{1,1},'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none'); hold on;
    %     s2 = scatter(v2(:,1),v2(:,2),50,plot_set{1,2},'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none'); hold on;
    % 
    %     daspect([1 1 1])
    %     grid on
    %     box off
    %     % ax1.XScale = 'log';
    %     % ax1.YScale = 'log';            
    %     axis square
    %     ax1 = gca;
    %     ax1.FontSize = 10;
    %     % ax1.XTick = [1 2 5 10 20];
    %     % ax1.YTick = ax1.XTick;
    %     ax.XLim = [-0.8 1];
    %     ax.YLim = [-0.3 0.1];
    %     ax.XTick = -1:0.2:1;
    %     ax.YTick = -1:0.05:1;
    %     xlabel(['Field anisotropy'])
    %     ylabel(['Behaviour anisotropy'])
    % 
    %     line([ax1.XLim(1) ax1.XLim(2)],[ax1.YLim(1) ax1.YLim(2)],'Color',[.5 .5 .5]);
    % 
    %     [r1,p1] = corr(v3(:,1),v3(:,2),'type','Pearson','rows','pairwise');
    % 
    %     s1f = polyfit(v3(:,1),v3(:,2),1);
    %     rf1 = refline(s1f(1),s1f(2));
    %     set(rf1,'Color','k','LineWidth',2);
    % 
    %     text(ax1,1.1,0.9,sprintf('{\\itr} = %.2f, {\\itr} = %.3f',r1,p1),'Units','normalized','FontSize',10,'Color','k')
    %     text(ax1,0,-0.1,sprintf('Elongated\nalong x-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','left')
    %     text(ax1,1,-0.1,sprintf('Elongated\nalong y-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','right')
    %     text(ax1,-0.1,0,sprintf('Biased\nalong x-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','right')
    %     text(ax1,-0.1,1,sprintf('Biased\nalong y-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','right')
    % 
    %     [~,leg] = legendflex([s1 s2],{'Arena 1','Hills'},'anchor',{'e','e'},'ncol',1,'box','off','buffer',[130,0],'xscale',0.5);  
    %     ytickformat('%.2f')
    %     xtickformat('%.1f')

%% >>>>>>>>>> Field anisotropy vs behaviour (case by case)
    var = 'field_anisotropy';
    v1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.f_rate_hz>0.1,1); % arena 1 data
    v1 = cell2mat(v1);
    v2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.1,1); % arena 1 data
    v2 = cell2mat(v2);
    v3 = clumaa.(var)(pidx & clumaa.partn==3 & clumaa.f_rate_hz>0.1,1); % arena 1 data
    v3 = cell2mat(v3);
    
    
    figure('Units','pixels','Position',[100 100 1000 900],'visible','on');
    ax = axes('Units','pixels','Position',[250 350 400 400]);
    
        s1 = scatter(v1(:,1),v1(:,2),50,plot_set{1,1},'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none'); hold on;
        s2 = scatter(v2(:,1),v2(:,2),50,plot_set{1,2},'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none'); hold on;
    
        daspect([1 1 1])
        grid on
        box off
        axis square
        ax1 = gca;
        ax1.FontSize = 10;
        ax.XLim = [-1 1];
        ax.YLim = [-1 1];
        ax.XTick = -1:0.2:1;
        ax.YTick = -1:0.2:1;
        ax.XTickLabel = {};
        ax.YTickLabel = {};
    
        line([ax1.XLim(1) ax1.XLim(2)],[ax1.YLim(1) ax1.YLim(2)],'Color',[.5 .5 .5]);
    
        [r1,p1] = corr(v1(:,1),v1(:,2),'type','Pearson','rows','pairwise');
        [r2,p2] = corr(v2(:,1),v2(:,2),'type','Pearson','rows','pairwise');
    
        s1f = polyfit(v1(:,1),v1(:,2),1);
        rf1 = refline(s1f(1),s1f(2));
        set(rf1,'Color',plot_set{1,1},'LineWidth',2);
        s2f = polyfit(v2(:,1),v2(:,2),1);
        rf2 = refline(s2f(1),s2f(2));
        set(rf2,'Color',plot_set{1,2},'LineWidth',2);
    % p1
    % p2
        text(ax1,0.48,1.02,sprintf('{\\itr} = %.2f, {\\itp} < .001',r1),'Units','normalized','FontSize',12,'Color',plot_set{1,1},'VerticalAlignment','bottom','HorizontalAlignment','right')
        text(ax1,0.52,1.02,sprintf('{\\itr} = %.2f, {\\itp} < .001',r2),'Units','normalized','FontSize',12,'Color',plot_set{1,2},'VerticalAlignment','bottom','HorizontalAlignment','left')
    
        [~,leg] = legendflex([s1 s2],{'Arena 1','Hills'},'anchor',{'se','se'},'ncol',1,'box','off','buffer',[0,0],'xscale',0.5,'FontSize',12);  

    ax_b = axes('Units','pixels','Position',[ax.Position(1)-60 ax.Position(2) 55 ax.Position(4)]); % behaviour y-axis plot
        xi = linspace(-1,1,bin_res);
        bs = 0.03;
        f1 = ksdensity(v1(:,2),xi,'Bandwidth',bs);
        f2 = ksdensity(v2(:,2),xi,'Bandwidth',bs);
    
        % plot(f1,xi,'Color',plot_set{1,1}); hold on;
        % plot(f2,xi,'Color',plot_set{1,2});
        alph = 0.75;
        patch([f1 zeros(size(f1))],[xi fliplr(xi)],plot_set{1,1},'FaceAlpha',alph,'EdgeColor','none'); hold on;
        patch([f2 zeros(size(f2))],[xi fliplr(xi)],plot_set{1,2},'FaceAlpha',alph,'EdgeColor','none'); 

        ylabel(['Behaviour anisotropy'])
        ax_b.YLim = [-1 1];
        ax_b.XTick = [];
        ax_b.YTick = ax.YTick;        
        ax_b.XColor = 'none';
        box off
        grid on
        ytickformat('%.1f')
        ax_b.FontSize = 12;
    
        text(ax,-0.25,0,sprintf('Biased\nalong x-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','right')
        text(ax,-0.25,1,sprintf('Biased\nalong y-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','right')
    
    ax_f = axes('Units','pixels','Position',[ax.Position(1) ax.Position(2)-60 ax.Position(3) 55]); % field, x-axis plot
        f1 = ksdensity(v1(:,1),xi,'Bandwidth',bs);
        f2 = ksdensity(v2(:,1),xi,'Bandwidth',bs);
    
        % plot(xi,f1,'Color',plot_set{1,1}); hold on;
        % plot(xi,f2,'Color',plot_set{1,2});
        patch([xi fliplr(xi)],[f1 zeros(size(f1))],plot_set{1,1},'FaceAlpha',alph,'EdgeColor','none'); hold on;
        patch([xi fliplr(xi)],[f2 zeros(size(f2))],plot_set{1,2},'FaceAlpha',alph,'EdgeColor','none'); 

        xlabel(['Field anisotropy'])
        ax_f.XLim = [-1 1];
        ax_f.YTick = [];
        ax_f.XTick = ax.XTick;
        ax_f.YColor = 'none';        
        box off
        grid on
        ytickformat('%.1f')
        ax_f.FontSize = 12;

        text(ax_f,0,-1,sprintf('Elongated\nalong x-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','left')
        text(ax_f,1,-1,sprintf('Elongated\nalong y-axis'),'Units','normalized','FontSize',10,'Color','k','HorizontalAlignment','right')

% figure
% subplot(1,2,1)
%     bin_res = 256;
%     bs = linspace(-1,1,bin_res);
%     [xx,yy] = ndgrid(bs,bs);
%     f = mvksdensity(v1,[yy(:) xx(:)],'BandWidth',0.07);
%     F = reshape(f,size(xx));
% 
%     imagesc('XData',xx(:),'YData',yy(:),'CData',F);
%     axis xy
%     daspect([1 1 1])
%     xlabel('Field anisotropy')
%     ylabel('Behaviour anisotropy')
%     ax = gca;
%     ax.XLim = [-1 1];
%     ax.YLim = [-1 1];
%     line([ax.XLim(1) ax.XLim(2)],[ax.YLim(1) ax.YLim(2)],'Color',[.5 .5 .5]);
% 
% subplot(1,2,2)
%     bin_res = 256;
%     bs = linspace(-1,1,bin_res);
%     [xx,yy] = ndgrid(bs,bs);
%     f = mvksdensity(v2,[yy(:) xx(:)],'BandWidth',0.07);
%     F = reshape(f,size(xx));
% 
%     imagesc('XData',xx(:),'YData',yy(:),'CData',F);
%     axis xy
%     daspect([1 1 1])
%     xlabel('Field anisotropy')
%     ylabel('Behaviour anisotropy')
%     ax = gca;
%     ax.XLim = [-1 1];
%     ax.YLim = [-1 1];
%     line([ax.XLim(1) ax.XLim(2)],[ax.YLim(1) ax.YLim(2)],'Color',[.5 .5 .5]);

%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\PIT_field_anisotropy.png']; 
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


 


















