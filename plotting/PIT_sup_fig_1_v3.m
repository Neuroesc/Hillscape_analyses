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
    fs = [15 10];
    f_r = 15;
    f_c = 6;
    tiledlayout(f_r,f_c,"TileSpacing","compact",'Padding','compact')


ucis_to_exclude = {'RG19_221013_t1_c10',...
    'RG19_221026_t1_c4',...
    'RG19_221013_t3_c13',...
    'RG26_230119_t3_c1',...
    'RG26_230119_t4_c5',...
    'RG19_221014_t2_c9',...
    'RG19_221014_t3_c10',...
    'RG26_230119_t1_c9',...
    'RG19_221021_t2_c1',...
    'RG26_230130_t2_c1',...
    'RG19_221012_t4_c11',...
    'RG19_221014_t3_c1',...
    'RG19_221013_t1_c5',...
    'RG19_221013_t3_c20',...
    'RG19_221012_t1_c3',...
    'RG19_221012_t2_c11',...
    'RG11_220519_t2_c3',...
    'RG19_221013_t1_c14',...
    'RG19_221013_t2_c11',...
    'RG19_221012_t2_c11',...
    'RG11_220519_t2_c4',...
    'RG19_221012_t1_c2'...
    'RG11_220520_t2_c3',...
    'RG11_220520_t2_c4',...
    'RG11_220520_t2_c8',...
    'RG19_221012_t1_c10'...
    'RG19_221012_t2_c10',...
    'RG26_230119_t1_c9',...
    'RG26_230114_t2_c2',...
    'RG19_221012_t3_c16',...
    'RG26_230126_t4_c1',...
    'RG19_221013_t4_c5',...
    'RG19_221014_t3_c16',...
    'RG26_230120_t1_c6',...
    'RG19_221014_t3_c12',...
    'RG19_221012_t3_c3',...
    'RG26_230130_t1_c7',...
    'RG11_220517_t2_c4',...
    'RG19_221012_t2_c13',...
    'RG26_230120_t3_c1',...
    'RG19_221014_t2_c1',...
    'RG19_221012_t4_c13',...
    'RG19_221012_t1_c7',...
    'RG26_230114_t2_c7'};

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example cells
    % find which cells to plot  
    ucis = unique(clumaa.uci(pidx));
    main_uci = ucis{236};
    idx = [1 ...
        90 103 112 118 123 127 131 133 143 145 146 149 152 163 165 171 172 ...
        7 36 37 43 46 54 56 57 59 67 76 77 81 82 84 ...
        1 12 21 22 23 31 30 32 34 37 40 42 50 51 52 53 55 64 71 72 87 89 91 93 94 96 100 ...
        102 110 113 116 117 119 120 122 125 129 130 138 141 142 154 163 161 166 167 168 169 170 173 180 197 ... 
        201 202 206 209 213 214 230 234 235 238 252 257 260 267 268 275 277 283 ...
        322 328 330 333 340 341 353 354 357 367 369 379 377 390 393 397 ...
        403 406 408 412 423 429 431 434 441 443 448 457 458 459 460 462 468 480 476 485 488 491 497 504];
    % idx = unique(idx);
    ucis = ucis(idx);

    ix = ismember(ucis,ucis_to_exclude);
    idx(ix) = [];
    length(idx) % should be 111
    % return
    ucis = unique(clumaa.uci(pidx));
    % ucis = ucis(unique(idx));
    ucis = ucis((idx));

% keyboard
    count = 0;
    fig_count = 1;
    disp(sprintf('Fig %d...',fig_count))    
    max_plots_now = f_r .* f_c;
    for ii = 1:numel(ucis)
    % for ii = 1:20
        if fig_count==1 && count==0
            count = count+18;

            disp(sprintf('\t...%s',main_uci))
            m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,main_uci) & clumaa.partn==1)};        
            m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,main_uci) & clumaa.partn==2)};
            m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,main_uci) & clumaa.partn==3)};
    
            rn = clumaa.rat{(ismember(clumaa.uci,main_uci) & clumaa.partn==1)}; 
            rs = {'RG6','RG11','RG19','RG26','RG40'}; % so order matches other figures
            [a,b] = ismember(rs,rn);
            rnum = find(a);

            nexttile([3,2])
                % plot map
                imagesc(m1,'alphadata',~isnan(m1)); hold on;
    
                % axis settings
                axis xy on
                daspect([1 1 1])
                ax1 = gca;
                colormap(gca,'turbo')   
                ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') max(m3(:),[],'omitnan') 1])];
                ax1.XTick = [];
                ax1.YTick = [];
                t = title(ax1,'A','FontWeight','normal','FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom','Color','k');
                ax1.TitleHorizontalAlignment='left';
                t.Position(2) = t.Position(2)*1.1;

                % text
                text(0,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',9,'Color','w','VerticalAlignment','top','rotation',0,'FontWeight','normal')
                text(0,0,sprintf('Max frate:\n%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',9,'Color','w','VerticalAlignment','bottom','FontWeight','normal')
                text(1,1,sprintf('Cell 236'),'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color','k','VerticalAlignment','bottom','FontWeight','normal')
                text(0,1,sprintf('Rat %d',rnum),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','bottom','FontWeight','normal')

            nexttile([3,2])
                % plot map
                imagesc(m2,'alphadata',~isnan(m2)); hold on;
    
                % axis settings
                axis xy on
                daspect([1 1 1])
                ax2 = gca;
                colormap(ax2,ax1.Colormap)   
                ax2.CLim = ax1.CLim;
                ax2.XTick = [];
                ax2.YTick = [];
    
                % additional plots
                ps = linspace(0,size(m1,2),7);
                tops = ps(2:2:end);
                plot([tops; tops],ax1.YLim,'w:','LineWidth',1)
                text(0,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',9,'Color','w','VerticalAlignment','top','rotation',0)
                text(0.825,0.5,sprintf('Ridge tops'),'Units','normalized','HorizontalAlignment','center','FontSize',12,'Color','w','VerticalAlignment','bottom','rotation',90)

            nexttile([3,2])
                % plot map
                imagesc(m3,'alphadata',~isnan(m3)); hold on;
    
                % axis settings
                ax3 = gca;
                axis xy on
                daspect([1 1 1])
                colormap(ax3,ax1.Colormap)   
                ax3.CLim = ax1.CLim;
                ax3.XTick = [];
                ax3.YTick = [];
                text(0,1,sprintf('%s',maze_names{3}),'Units','normalized','HorizontalAlignment','left','FontSize',9,'Color','w','VerticalAlignment','top','rotation',0)

        else
            uci = ucis{ii};
            disp(sprintf('\t...%s',uci))
            m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};        
            m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
            m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};
    
            rn = clumaa.rat{(ismember(clumaa.uci,uci) & clumaa.partn==1)}; 
            rs = {'RG6','RG11','RG19','RG26','RG40'}; % so order matches other figures
            [a,b] = ismember(rs,rn);
            rnum = find(a);
    
            % arena 1 ratemap
            nexttile           
                % plot map
                imagesc(m1,'alphadata',~isnan(m1)); hold on;
    
                % axis settings
                axis xy on
                daspect([1 1 1])
                ax1 = gca;
                colormap(gca,'turbo')   
                ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') max(m3(:),[],'omitnan') 1])];
                ax1.XTick = [];
                ax1.YTick = [];
    
                % text
                % text(0,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','k','VerticalAlignment','bottom','rotation',0)
                text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','w','VerticalAlignment','bottom')
                text(1,1,sprintf('Cell %d',idx(ii)),'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color','k','VerticalAlignment','bottom')
                text(0,1,sprintf('Rat %d',rnum),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','bottom')
                if count==18 && ii==2
                    t = title(ax1,'B','FontWeight','normal','FontSize',20,'HorizontalAlignment','left','VerticalAlignment','bottom','Color','k');
                    ax1.TitleHorizontalAlignment='left';
                    t.Position(2) = t.Position(2)*0.95;
                end
% keyboard
            % hills ratemap
            nexttile
                % plot map
                imagesc(m2,'alphadata',~isnan(m2)); hold on;
    
                % axis settings
                axis xy on
                daspect([1 1 1])
                ax2 = gca;
                colormap(ax2,ax1.Colormap)   
                ax2.CLim = ax1.CLim;
                ax2.XTick = [];
                ax2.YTick = [];
    
                % additional plots
                ps = linspace(0,size(m1,2),7);
                tops = ps(2:2:end);
                plot([tops; tops],ax1.YLim,'w:','LineWidth',1)
    
                % text
                % text(0,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','k','VerticalAlignment','bottom','rotation',0)
                % text(0,0,sprintf('%s',uci),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','k','VerticalAlignment','bottom','rotation',90,'Interpreter','none')
    
            % arena 2 ratemap
            % ax3 = axes('Units','pixels','Position',[x(jj)+ax1.Position(3)-(ax1.Position(3).*0.48) y(jj)+ax1.Position(4)+ysiz(2) ax1.Position(3).*0.48 ax1.Position(4).*0.48]);
            nexttile
                % plot map
                imagesc(m3,'alphadata',~isnan(m3)); hold on;
    
                % axis settings
                ax3 = gca;
                axis xy on
                daspect([1 1 1])
                colormap(ax3,ax1.Colormap)   
                ax3.CLim = ax1.CLim;
                ax3.XTick = [];
                ax3.YTick = [];

            count = count+3;

            % keyboard
        end

        % count how many plots we have done
        if count==max_plots_now || ii == numel(ucis) % max we can fit on a page
            disp(sprintf('...saving'))      
            if fig_count==1
                % horizontal colormap
                axc = axes('Units','pixels','Position',[470 830 95 12]);
                    xn = flipud(linspace(ax1.CLim(1),ax1.CLim(2),100));
                    imagesc(xn,'XData',ax1.CLim,'YData',[0 1]);
                    colormap(axc,ax1.Colormap);
                    axis on
                    axc.XTick = [];
                    axc.YTick = [];
                    axc.FontSize = 7;
                    text(1,0.5,sprintf(' Max'),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8,'Units','normalized')
                    text(0,0.5,sprintf('0 Hz '),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',8,'Units','normalized')
            end
            % keyboard
            if 1
                fname = [config.fig_dir '\' sprintf('Fig S2_%d.png',fig_count)]; 
                exportgraphics(gcf,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',res);  
                clf(gcf);   
                tiledlayout(f_r,f_c,"TileSpacing","compact",'Padding','compact')
                count = 0;
                fig_count = fig_count+1;
            end
            disp(sprintf('Fig %d...',fig_count))    
        end
    end
% keyboard
% return
    close(gcf);
















