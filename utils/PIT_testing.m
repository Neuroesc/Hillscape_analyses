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




































return
close all;


clumaa = PIT_correlations(config,pidx,clumaa,posdata,0);







% h1 = clumaa.hd_3d_ratemap(pidx & clumaa.partn==1); % arena 1 data
% h1 = cat(3,h1{:});
% h2 = clumaa.hd_3d_ratemap(pidx & clumaa.partn==2); % hills data
% h2 = cat(3,h2{:});
% 
% 
% 
% figure
% subplot(1,2,1)
% m1 = mean(h1,3,'omitmissing');
% imagesc(m1)
% daspect([1 1 1])
% 
% subplot(1,2,2)
% m2 = mean(h2,3,'omitmissing');
% imagesc(m2)
% daspect([1 1 1])



figure
        var = 'pitch_stability';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1
        v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % hills
        % v4 = squeeze(shuffs_within(:,1,1)); % arena shuffles
        % v5 = squeeze(shuffs_within(:,1,2)); % hills shuffles

        xi = -1:0.01:1;
        bw = 0.05;
        f1 = ksdensity(v1(:),xi(:),"Bandwidth",bw);
        f2 = ksdensity(v2(:),xi(:),"Bandwidth",bw);
        % f4 = ksdensity(v4(:),xi(:),"Bandwidth",bw);
        % f5 = ksdensity(v5(:),xi(:),"Bandwidth",bw);

        % main plot   
        alph = 0.7;          
        a1 = area(xi,f1,'FaceColor',plot_set{1,1},'EdgeColor','k','FaceAlpha',alph); hold on;
        a2 = area(xi,f2,'FaceColor',plot_set{1,2},'EdgeColor','k','FaceAlpha',alph);       
        % a4 = area(xi,f4,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--'); hold on;       
        % a5 = area(xi,f5,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle',':'); hold on;       

        % axis settings
        ax = gca;
        ax.XLim = [-0.5 1]; 
        ax.FontSize = 10;
        set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))
        box off
        xlabel(sprintf('Correlation ({\\itr})'))    
        ylabel(sprintf('PDF'))    

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 






















keyboard
return


        % real data
        ucis = unique(clumaa.uci(pidx));
        dat = NaN(95,211,length(ucis));
        dat_control = NaN(95,191,length(ucis));
        for ii = 1:length(ucis)
            % normal correlations
            m1 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1}; % arena 1
            m2 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==2}; % hills
            m2s = clumaa.ratemap_surficial{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==2}; % hills, surficial map                    
            m3 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==3}; % arena 2
            amap = ndautoCORR(m1,m2s);
            dat(:,:,ii) = imresize(amap,[95 211],'nearest');
            amap = ndautoCORR(m1,m3);
            dat_control(:,:,ii) = imresize(amap,[95 191],'nearest');
        end

        % shuffles
        iti = length(ucis);
        dats = NaN(95,211,length(ucis));
        for ii = 1:iti
            uix = randperm(numel(ucis),2); % 2 random place cells
            m1 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{uix(1)}) & clumaa.partn==1}; % arena 1, cell 1
            m2 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==2}; % hills, cell 2
            m2s = clumaa.ratemap_surficial{ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==2}; % hills, cell 2    
            m3 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{uix(1)}) & clumaa.partn==3}; % arena 1, cell 1 
            amap = ndautoCORR(m1,m2s);
            dats(:,:,ii) = imresize(amap,[95 211],'nearest');
        end

        % mean maps
        mean_map = mean(dat,3,'omitnan');
        mean_map_control = mean(dat_control,3,'omitnan');        
        mean_shuffle = mean(dats,3,'omitnan');


        figure
        subplot(2,2,1)
        imagesc(mean_map)
        daspect([1 1 1])
        colorbar;
        ax = gca;
        ax.CLim = [0 0.1];

        subplot(2,2,2)
        imagesc(mean_shuffle)
        daspect([1 1 1])
        colorbar;
        ax = gca;
        ax.CLim = [0 0.1];

        subplot(2,2,3)
        m = mean_map-mean_shuffle;
        imagesc(m); hold on;
        daspect([1 1 1])
        colorbar;
        ax = gca;
        ax.CLim = [0 0.1];

        mx = size(m,2)/2;
        line([mx mx],ax.YLim,'Color','w')
        line([mx mx]+20,ax.YLim,'Color','w','LineStyle',':')
        line([mx mx]-20,ax.YLim,'Color','w','LineStyle',':')
        my = size(m,1)/2;
        line(ax.XLim,[my my],'Color','w')

        subplot(2,2,4)
        imagesc(mean_map_control)
        daspect([1 1 1])
        colorbar;
        ax = gca;
        ax.CLim = [0 0.8];


keyboard

























return
    % Create figure
    % fig_now = figure('Units','pixels','Position',[100 100 900 800],'visible','on');
    % set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    % set(gcf,'color','w'); % makes the background colour white


        var = 'ratemap_planar_half';
        v1_h1 = clumaa.(var)(pidx & clumaa.partn==1,1); 
        v1_h2 = clumaa.(var)(pidx & clumaa.partn==1,2); 
        v2_h1 = clumaa.(var)(pidx & clumaa.partn==2,1); 
        v2_h2 = clumaa.(var)(pidx & clumaa.partn==2,2); 
        v3_h1 = clumaa.(var)(pidx & clumaa.partn==3,1); 
        v3_h2 = clumaa.(var)(pidx & clumaa.partn==3,2); 

        v1h1 = cat(3,v1_h1{:});
        v1h2 = cat(3,v1_h2{:});
        % v1_diff = (v1h1 - v1h2) ./ (v1h1 + v1h2);
        % v1_mean_diff = mean(v1_diff,3,'omitnan');
        v1_mean_diff = NaN(size(v1h1(:,:,1)));
        for jj = 1:numel(v1h1(:,:,1))
            [r,c] = ind2sub(size(v1h1(:,:,1)),jj);
            v1_mean_diff(jj) = corr(reshape(v1h1(r,c,:),[],1),reshape(v1h2(r,c,:),[],1),'rows','pairwise','type','Pearson');
        end

        v2h1 = cat(3,v2_h1{:});
        v2h2 = cat(3,v2_h2{:});
        % v2_diff = (v2h1 - v2h2) ./ (v2h1 + v2h2);
        % v2_mean_diff = mean(v2_diff,3,'omitnan');
        v2_mean_diff = NaN(size(v2h1(:,:,1)));
        for jj = 1:numel(v2h1(:,:,1))
            [r,c] = ind2sub(size(v2h1(:,:,1)),jj);
            v2_mean_diff(jj) = corr(reshape(v2h1(r,c,:),[],1),reshape(v2h2(r,c,:),[],1),'rows','pairwise','type','Pearson');
        end

        v3h1 = cat(3,v3_h1{:});
        v3h2 = cat(3,v3_h2{:});
        % v2_diff = (v2h1 - v2h2) ./ (v2h1 + v2h2);
        % v2_mean_diff = mean(v2_diff,3,'omitnan');
        v3_mean_diff = NaN(size(v3h1(:,:,1)));
        for jj = 1:numel(v3h1(:,:,1))
            [r,c] = ind2sub(size(v3h1(:,:,1)),jj);
            v3_mean_diff(jj) = corr(reshape(v3h1(r,c,:),[],1),reshape(v3h2(r,c,:),[],1),'rows','pairwise','type','Pearson');
        end

figure
subplot(1,3,1)
imagesc(v1_mean_diff); hold on;
daspect([1 1 1])
colormap(turbo)
colorbar
ax = gca;
c = [0.7 1];
ax.CLim = c;

subplot(1,3,2)
imagesc(v2_mean_diff); hold on;
daspect([1 1 1])
colormap(turbo)
colorbar

ax = gca;
ps = linspace(0,size(v2_mean_diff,2),7);
tops = ps(2:2:end);
plot([tops; tops],ax.YLim,'w:')
ax.CLim = c;

subplot(1,3,3)
imagesc(v3_mean_diff); hold on;
daspect([1 1 1])
colormap(turbo)
colorbar
ax = gca;
c = [0.7 1];
ax.CLim = c;

keyboard
        var = 'planar_field_aspect_ratio';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); 
        var = 'repetition_score';
        v2 = clumaa.(var)(pidx & clumaa.partn==2,1); 


ucis = clumaa.uci(pidx);
d = NaN(length(ucis),2);
for uu = 1:length(ucis)
    m1 = clumaa.planar_field_aspect_ratio(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==1,1); 
    m2 = clumaa.planar_field_aspect_ratio(ismember(clumaa.uci,ucis{uu}) & clumaa.partn==2,1); 
    d(uu,:) = [m1 m2];
end


figure
scatter(d(:,1),d(:,2),50,'k','filled')
























%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> 
    xnow = 60;
    ynow = 500;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 400 200]);
        PIT_plot_mazes(ax,1); % plot arena

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow-250 400 200]);
        PIT_plot_mazes(ax,2); % plot hills

%% >>>>>>>>>> Save the overall figure
    if 0
        fname = [config.fig_dir '\PIT_testing.png']; 
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

close all
figure
pid = 2;
p = posdata.pos{pid}; % [session, pot, pox, poy, poh, rx, ry, gx, gy, bx, by, poz, pox_planar, poy_planar, poz_planar, pox_surficial, poy_surficial, poz_surficial, poz_curve]

sidx = p.session==2;
c = cline(p.pox(sidx),p.poy(sidx),-p.poz(sidx),-p.poz(sidx));
set(c,'LineWidth',2)
daspect([1 1 1])
view(3)
colormap('turbo')
axis off

ax = gca;
ax.CLim = [100 500];
camlight(ax)
% lighting(ax,'phong')
ax.Projection = 'perspective';

    if 1
        fname = [config.fig_dir '\PIT_testing_traj.png']; 
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















