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

%% >>>>>>>>>> Field to wall angle, expected maps
    % collect data
    ucis = unique(clumaa.uci(pidx));
    datn = cell(1,3);
    for pp=1:3    
        idx = ismember(clumaa.uci,ucis{1}) & clumaa.partn==pp;
        if pp==2
            rmap = clumaa.ratemap_surficial{idx};       
        else
            rmap = clumaa.ratemap_planar{idx};                       
        end

        b = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0]+0.5; 
        a = [90 0 90 0];
        gpoly = [];
        res = 1000;
        for jj = 1:size(b,1)-1
            gpoly = [gpoly; linspace(b(jj,1),b(jj+1,1),res)' linspace(b(jj,2),b(jj+1,2),res)' repmat(a(jj),res,1)];
        end
        gpoly = unique(round(gpoly),'rows');

        ps = linspace(0,size(rmap,2),7);
        tops = ps(2:2:end);
        bottoms = ps(3:2:end-2);
        b2 = [0 0;0 size(rmap,1); size(rmap,2) size(rmap,1);size(rmap,2) 0;0 0;tops(1) 0;tops(1) size(rmap,1);tops(2) size(rmap,1);tops(2) 0;tops(3) 0;tops(3) size(rmap,1);tops(3) 0;bottoms(1) 0;bottoms(1) size(rmap,1);bottoms(2) size(rmap,1);bottoms(2) 0]+0.5; 
        a2n = [90 0 90 0 0 90 0 90 0 90 90 0 90 0 90];
        gpoly2 = [];
        res = 1000;
        for jj = 1:size(b2,1)-1
            gpoly2 = [gpoly2; linspace(b2(jj,1),b2(jj+1,1),res)' linspace(b2(jj,2),b2(jj+1,2),res)' repmat(a2n(jj),res,1)];
        end
        gpoly2 = unique(round(gpoly2),'rows');

        dmap = NaN([size(rmap),3,2]);

        angle_sigma = 10; % bigger = distal walls have more of an effect
        dist_sigma = 16; % bigger = walls have more distal effects on anisotropy
        for jj = 1:numel(dmap(:,:,1))
            %% arena
            % current coordinate
            [y,x] = ind2sub(size(dmap),jj);

            % distance to all walls
            ds = pdist2([x,y],gpoly(:,[1:2]),'euclidean');

            % wall angles
            ws = gpoly(:,3);
            weights = normpdf(ds,0,angle_sigma); % Gaussian weight wall angles
            mean_angle = sum(weights(:) .* ws(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');
            mean_distance = sum(weights(:) .* ds(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');

            dmap(y,x,1,1) = (mean_angle./45)-1;
            dmap(y,x,2,1) = normpdf(mean_distance,0,dist_sigma);
            dmap(y,x,3,1) = dmap(y,x,1,1) .* dmap(y,x,2,1);

            %% hills
            % current coordinate
            [y,x] = ind2sub(size(dmap),jj);

            % distance to all walls
            ds = pdist2([x,y],gpoly2(:,[1:2]),'euclidean');

            % wall angles
            ws = gpoly2(:,3);
            weights = normpdf(ds,0,angle_sigma); % Gaussian weight wall angles
            mean_angle = sum(weights(:) .* ws(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');
            mean_distance = sum(weights(:) .* ds(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');

            dmap(y,x,1,2) = (mean_angle./45)-1;
            dmap(y,x,2,2) = normpdf(mean_distance,0,dist_sigma);
            dmap(y,x,3,2) = dmap(y,x,1,2) .* dmap(y,x,2,2);
        end
        datn(pp) = { dmap };
    end

%% >>>>>>>>>> Rat in arena, weighted walls
    xnow = 50;
    ynow = 600;
    xsiz = 125;
    ysiz = 70;
    xbuff = 20;
    ybuff = 100;
    cmap1 = viridis(256);
    cmap = cmocean('balance');

    % plot arena
    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);    
        % ah = add_panel_title('g',sprintf('Anisotropy is not explained\nby geometry alone'),'yoffset',0,'xoffset',20,'width',400); 

        p1 = plot(gpoly(:,1),gpoly(:,2),'k.','LineStyle','none','MarkerSize',10); hold on;

        % axis settings
        axis off xy
        daspect([1 1 1])

        % rat scale inset
        x = 50;
        y = 15;        
        h = 400;
        scale = (ax1.Position(4)/1500)*h;
        ir = imread('C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated media\rat_silh_v2.png');
        w = scale*(size(ir,2)/size(ir,1));
        h = scale;

        xp = ax1.Position(3)./ax1.XLim(2);
        yp = ax1.Position(4)./ax1.YLim(2);
        axp = axes('Units','pixels','Position',[ax1.Position(1)+(x*xp)-(w/2) ax1.Position(2)+(y*yp)-(scale/2) w h]); 
            image(ir);
            axis off ij
            daspect([1 1 1])

    % Gaussian weights
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz ysiz]);    
        % ah = add_panel_title('g',sprintf('Anisotropy is not explained\nby geometry alone'),'yoffset',0,'xoffset',20,'width',400); 

        % distance to all walls
        ds = pdist2([x,y],gpoly(:,[1:2]),'euclidean');

        % wall angles
        ws = gpoly(:,3);
        weights = normpdf(ds,0,angle_sigma); % Gaussian weight wall angles
        mean_angle = sum(weights(:) .* ws(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');
        mean_distance = sum(weights(:) .* ds(:),'all','omitmissing') / sum(weights(:),'all','omitmissing');

        % plot geometry/weights
        s1 = scatter(gpoly(:,1),gpoly(:,2),30,weights,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5); hold on;

        % axis settings
        axis off xy
        daspect([1 1 1])
        colormap(ax2,cmap1)

        % rat scale inset        
        xp = ax2.Position(3)./ax1.XLim(2);
        yp = ax2.Position(4)./ax1.YLim(2);
        axp = axes('Units','pixels','Position',[ax2.Position(1)+(x*xp)-(w/2) ax2.Position(2)+(y*yp)-(scale/2) w h]); 
            image(ir);
            axis off ij
            daspect([1 1 1])

    % walla angles
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+xbuff ynow xsiz ysiz]);    
        % ah = add_panel_title('g',sprintf('Anisotropy is not explained\nby geometry alone'),'yoffset',0,'xoffset',20,'width',400); 

        % plot geometry/weights
        s1 = scatter(gpoly(:,1),gpoly(:,2),30,gpoly(:,3),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5); hold on;

        % axis settings
        axis off xy
        daspect([1 1 1])
        colormap(ax3,cmap1)

        % rat scale inset        
        xp = ax3.Position(3)./ax3.XLim(2);
        yp = ax3.Position(4)./ax3.YLim(2);
        axp = axes('Units','pixels','Position',[ax3.Position(1)+(x*xp)-(w/2) ax3.Position(2)+(y*yp)-(scale/2) w h]); 
            image(ir);
            axis off ij
            daspect([1 1 1])

    % walla angles
    ax3 = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-ybuff xsiz ysiz]);    
        % ah = add_panel_title('g',sprintf('Anisotropy is not explained\nby geometry alone'),'yoffset',0,'xoffset',20,'width',400); 

        % plot data
        mnow = dmap(:,:,1,1);
        imagesc(mnow,'alphadata',~isnan(mnow)); hold on;

        % axis settings
        axis xy on
        daspect([1 1 1])
        colormap(ax3,cmap1)
        ax1.XTick = [];
        ax1.YTick = [];

    % walla angles
    ax4 = axes('Units','pixels','Position',[ax3.Position(1)+ax3.Position(3)+xbuff ax3.Position(2) xsiz ysiz]);    
        % ah = add_panel_title('g',sprintf('Anisotropy is not explained\nby geometry alone'),'yoffset',0,'xoffset',20,'width',400); 

        % plot data
        mnow = dmap(:,:,2,1);
        imagesc(mnow,'alphadata',~isnan(mnow)); hold on;

        % axis settings
        axis xy on
        daspect([1 1 1])
        colormap(ax4,cmap1)
        ax1.XTick = [];
        ax1.YTick = [];

    % walla angles
    ax5 = axes('Units','pixels','Position',[ax4.Position(1)+ax4.Position(3)+xbuff ax4.Position(2) xsiz ysiz]);    
        % ah = add_panel_title('g',sprintf('Anisotropy is not explained\nby geometry alone'),'yoffset',0,'xoffset',20,'width',400); 

        % plot data
        mnow = dmap(:,:,3,1);
        imagesc(mnow,'alphadata',~isnan(mnow)); hold on;

        % axis settings
        axis xy on
        daspect([1 1 1])
        colormap(ax5,cmap)
        ax1.XTick = [];
        ax1.YTick = [];





























