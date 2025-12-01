%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  analysis function written for:
% Grieves, Duvelle and Taube (2024) 

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
    fs = [15 10];

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    clumaa.rat = cell(size(clumaa,1),1);
    for ii = 1:size(clumaa,1)
        uci = clumaa.uci{ii};
        info = strsplit(uci,'_');
        clumaa.rat{ii} = info{1};
    end

%% >>>>>>>>>> Cells per rat
    xnow = 50;
    ynow = 650;

    % cells per rat
    ax = axes('Units','pixels','Position',[xnow ynow 200 120]);
        ah = add_panel_title('A',sprintf('Place cells recorded per rat'),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);                

        ris = clumaa.rat(pidx & clumaa.partn==1);
        % rs = unique(ris);
        rs = {'RG6','RG11','RG19','RG26','RG40'}; % so order matches other figures
        c_num = NaN(length(rs),1);
        for rr = 1:length(rs) % for each rat
            c_num(rr) = sum(ismember(ris,rs{rr}));
        end
% keyboard
        % plot data
        b = bar(1:length(rs),c_num,0.5,'k');

        % axis settings
        ax.XLim = [0.5 length(c_num)+0.5];
        box off
        ylabel('Place cells')
        xlabel('Rat')

        % text
        for rr = 1:length(rs) % for each rat
            text(rr,c_num(rr)+3,sprintf('%.1f%%',c_num(rr)./sum(c_num(:)).*100),'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','bottom');
        end

%% >>>>>>>>>> Sessions per rat
    xnow = xnow+270;

    % sessions per rat
    ax = axes('Units','pixels','Position',[xnow ynow 200 120]);
        ah = add_panel_title('B',sprintf('Sessions recorded per rat'),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);                
    
        % ris = clumaa.rat(pidx & clumaa.partn==1);
        % rs = unique(ris);
        s_num = NaN(length(rs),1);
        for rr = 1:length(rs) % for each rat
            s_num(rr) = length( unique( clumaa.session_name( ismember(clumaa.rat,rs{rr}) & clumaa.partn==1 ) ) );
        end

        % plot data
        b = bar(1:length(rs),s_num,0.5,'k');

        % axis settings
        ax.XLim = [0.5 length(s_num)+0.5];
        box off
        ylabel('Sessions')
        xlabel('Rat')

        % text
        for rr = 1:length(rs) % for each rat
            text(rr,s_num(rr)+0.25,sprintf('%.1f%%',s_num(rr)./sum(s_num(:)).*100),'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','bottom');
        end

%% >>>>>>>>>> Cell activity per maze
    xnow = 50;
    ynow = ynow-270;

    % cells per maze
    ax = axes('Units','pixels','Position',[xnow ynow 180 180]);
        ah = add_panel_title('C',sprintf('Active place cells per maze'),'yoffset',0,'xoffset',0,'width',400,'fontsize',fs);                

        p_ucis = unique(clumaa.uci(pidx),'stable');
        f = NaN(numel(p_ucis),2);
        active_cutoff_hz = 0.5;
        for uu = 1:numel(p_ucis)
            f(uu,1) = clumaa.f_rate_hz( ismember(clumaa.uci,p_ucis{uu}) & clumaa.partn==1 ) > active_cutoff_hz;
            f(uu,2) = clumaa.f_rate_hz( ismember(clumaa.uci,p_ucis{uu}) & clumaa.partn==2 ) > active_cutoff_hz;
        end

        v = [sum( f(:,1) & ~f(:,2) ), sum( f(:,1) & f(:,2) ), sum( ~f(:,1) & f(:,2) )];

        % work out square properties
        vp = [v(1)+v(2) v(2) v(2)+v(3)];
        alph = 0.5;
        s1 = sqrt(vp(1));
        s2 = sqrt(vp(2));
        s3 = sqrt(vp(3));
        
        % plot data
        p1 = patch([-s1 -s1 0 0],[0 s1 s1 0],plot_set{1,1},'EdgeColor','none','FaceAlpha',alph); hold on; % bottom right corner is 0,0
        p3 = patch([0 0 s3 s3]-s2,[-s3 0 0 -s3]+s2,plot_set{1,2},'EdgeColor','none','FaceAlpha',alph);
        
        % axis settings
        daspect([1 1 1])
        axis tight off
        
        % text
        fsiz = 9;
        text(-s1,s1,sprintf('  %s only (N = %d)',maze_names{1},v(1)),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fsiz,'Color',p1.FaceColor);
        text(s3-s2,-s3+s2,sprintf('%s only (N = %d)  ',maze_names{2},v(3)),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',fsiz,'Color',p3.FaceColor);
        text(mean([s3-s2,-s1]),mean([s1,-s3+s2]),sprintf('Both\n(N = %d)',v(2)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fsiz);
% keyboard
%% >>>>>>>>>> Fields per maze
    xnow = xnow+260;

    % cells per maze
    ax = axes('Units','pixels','Position',[xnow ynow 180 180]);
        ah = add_panel_title('D',sprintf('Place cells with at least 1 place field'),'yoffset',0,'xoffset',10,'width',400,'fontsize',fs);                

        p_ucis = unique(clumaa.uci(pidx),'stable');
        f = NaN(numel(p_ucis),2);
        for uu = 1:numel(p_ucis)
            f(uu,1) = clumaa.planar_fields( ismember(clumaa.uci,p_ucis{uu}) & clumaa.partn==1 ) > 1;
            f(uu,2) = clumaa.planar_fields( ismember(clumaa.uci,p_ucis{uu}) & clumaa.partn==2 ) > 1;
        end

        v = [sum( f(:,1) & ~f(:,2) ), sum( f(:,1) & f(:,2) ), sum( ~f(:,1) & f(:,2) )];

        % work out square properties
        vp = [v(1)+v(2) v(2) v(2)+v(3)];
        alph = 0.5;
        s1 = sqrt(vp(1));
        s2 = sqrt(vp(2));
        s3 = sqrt(vp(3));
        
        % plot data
        p1 = patch([-s1 -s1 0 0],[0 s1 s1 0],plot_set{1,1},'EdgeColor','none','FaceAlpha',alph); hold on; % bottom right corner is 0,0
        p3 = patch([0 0 s3 s3]-s2,[-s3 0 0 -s3]+s2,plot_set{1,2},'EdgeColor','none','FaceAlpha',alph);
        
        % axis settings
        daspect([1 1 1])
        axis tight off
        
        % text
        text(-s1,s1,sprintf('  %s only (N = %d)',maze_names{1},v(1)),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fsiz,'Color',p1.FaceColor);
        text(s3-s2,-s3+s2,sprintf('%s only (N = %d)  ',maze_names{2},v(3)),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',fsiz,'Color',p3.FaceColor);
        text(mean([s3-s2,-s1]),mean([s1,-s3+s2]),sprintf('Both\n(N = %d)',v(2)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fsiz);

%% >>>>>>>>>> Histology
    xnow = 50;
    ynow = ynow-390;

    % cells per maze
    ax = axes('Units','pixels','Position',[xnow ynow 190 400]);
        ah = add_panel_title('E',sprintf('Histology'),'yoffset',-80,'xoffset',0,'width',400,'fontsize',fs);                
    
        iname = [config.data_out_dir 'Histology\R40.png'];
        [img1,m] = imread(iname,'png');
        img1 = img1(:,:,1:3);
        imshow(img1(:,:,1:3)); hold on;

        text(0,1,sprintf('Rat 4'),'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','bottom');

    xvec = [xnow+220 xnow+320 xnow+220 xnow+320];
    yvec = [ynow+165 ynow+165 ynow+35 ynow+35];
    nmes = {'R11','R19','R26','R6'};
    rn = [1 2 3 5];
    for ii = 1:4
        ax = axes('Units','pixels','Position',[xvec(ii) yvec(ii) 88 200]);
            iname = [config.data_out_dir 'Histology\' nmes{ii} '.png'];
            [img1,m] = imread(iname,'png');
            img1 = img1(:,:,1:3);
    
            imshow(img1(:,:,1:3)); hold on;

            text(0,1,sprintf('Rat %d',rn(ii)),'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','bottom');            
    end
% keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S1.png']; 
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



























