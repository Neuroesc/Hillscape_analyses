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

    if 0
        enames = {'3000x1500_arena','3000x1500_hills'};
        scratch_space = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2021 Place cells irregular terrain\associated data\bvc_model'; 
        model_maps = cell(3,length(enames));
        for ee = 1:length(enames)
            hname = [scratch_space '\20240221T093430\' enames{ee} '_512b_1028h_hpcs'];   
            load(hname,'hpc_maps','-mat');
            model_maps{1,ee} = hpc_maps;
        end
    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example cells
    xnow = 20;
    ynow = 750;
    xsiz = [105, 15]; % size, buffer
    ysiz = [60, 0]; % size, buffer
    xvec = xnow : sum(xsiz) : 10000;
    xvec = xvec(1:5);
    yvec = ynow : -sum(ysiz).*2.15 : -10000;
    yvec = yvec(1:8);
    [x,y] = ndgrid(xvec,yvec);
    max_plots_now = length(xvec).*length(yvec);

    % find which cells to plot  
    main_uci = 7;
    idx = [2 8 11 22 32 36 38 41 42 49 52 55 59 68 75 76 81 84 87 91 95 97 98 99 101];

    count = 0;
    fig_count = 1;
    disp(sprintf('Fig %d...',fig_count))    
    jj = 1;
    for ii = 1:numel(idx)
        if ii==1 % first cell
            disp(sprintf('\t...%d',main_uci))
            m1 = model_maps{1,1}(:,:,main_uci);        
            m2 = model_maps{1,2}(:,:,main_uci);
            % m3 = model_maps{1,1}(:,:,main_uci);

            max_plots_now = 25;

            % hills ratemap
            coef = 2.5;
            ysiz2 = ysiz*coef;
            xsiz2 = xsiz*coef;
            xnow2 = x(2)+170;
            ynow2 = y(1,2)+80;
            ax1 = axes('Units','pixels','Position',[xnow2 ynow2 xsiz2(1) ysiz2(1)]);
                % plot map
                imagesc(m2,'alphadata',~isnan(m2)); hold on;
    
                % axis settings
                axis xy on
                daspect([1 1 1])
                colormap(ax1,turbo)   
                ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') max(m3(:),[],'omitnan') 1])];
                ax1.XTick = [];
                ax1.YTick = [];
    
                % additional plots
                ps = linspace(0,size(m1,2),7);
                tops = ps(2:2:end);
                plot([tops; tops],ax1.YLim,'w:','LineWidth',1.5)
    
                % text
                text(0,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color','k','VerticalAlignment','bottom','rotation',0)
                % text(-0.01,0,sprintf('%s',uci),'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color','k','VerticalAlignment','bottom','rotation',90,'Interpreter','none')
                text(0,0,sprintf('Max:\n%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color','w','VerticalAlignment','bottom')
                text(0.825,0.5,sprintf('Ridge tops'),'Units','normalized','HorizontalAlignment','center','FontSize',12,'Color','w','VerticalAlignment','bottom','rotation',90)
                text(1,-0.02,sprintf('Cell %d',main_uci),'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color','k','VerticalAlignment','top')

                % vertical colormap
                axc = axes('Units','pixels','Position',[xnow2+ax1.Position(3)+20 ynow2+25 12 ax1.Position(4)*0.7]);
                    xn = flipud(linspace(ax1.CLim(1),ax1.CLim(2),100)');
                    imagesc(xn,'YData',ax1.CLim,'XData',[0 1]);
                    colormap(axc,ax1.Colormap);
                    axis on
                    axc.XTick = [];
                    axc.YTick = [];
                    axc.FontSize = 7;
                    text(0.5,1,sprintf('Max'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Units','normalized')
                    text(0.5,0,sprintf('0 Hz'),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'Units','normalized')

            % arena 1 ratemap
            ax2 = axes('Units','pixels','Position',[xnow2-ax1.Position(3)-5 ynow2 ax1.Position(3) ax1.Position(4)]);
                ah = add_panel_title('A',sprintf(''),'yoffset',-20,'xoffset',-10,'width',400,'fontsize',fs); 
            
                % plot map
                imagesc(m1,'alphadata',~isnan(m1)); hold on;
    
                % axis settings
                axis xy on
                daspect([1 1 1])
                colormap(ax2,ax1.Colormap)   
                ax2.CLim = ax1.CLim;
                ax2.XTick = [];
                ax2.YTick = [];
    
                % text
                text(-0.02,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color','k','VerticalAlignment','bottom','rotation',0)
    
            % % arena 2 ratemap
            % ax3 = axes('Units','pixels','Position',[xnow2+ax1.Position(3)-(ax1.Position(3).*0.48) ynow2+ax1.Position(4)+ysiz2(2)+5 ax1.Position(3).*0.48 ax1.Position(4).*0.48]);
            %     % plot map
            %     imagesc(m3,'alphadata',~isnan(m3)); hold on;
            % 
            %     % axis settings
            %     axis xy on
            %     daspect([1 1 1])
            %     colormap(ax3,ax1.Colormap)   
            %     ax3.CLim = ax1.CLim;
            %     ax3.XTick = [];
            %     ax3.YTick = [];
            % 
            %     % text
            %     text(0,1,sprintf('%s',maze_names{3}),'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color','k','VerticalAlignment','bottom','rotation',0)

            % shift next plots down
            x = x(:,3:end);
            y = y(:,3:end);

        end 

        disp(sprintf('\t...%d',idx(ii)))
        m1 = model_maps{1,1}(:,:,idx(ii));        
        m2 = model_maps{1,2}(:,:,idx(ii));
        % m3 = model_maps{1,1}(:,:,idx(ii));

        % hills ratemap
        yoffset = 50;
        ax1 = axes('Units','pixels','Position',[x(jj) y(jj)+yoffset xsiz(1) ysiz(1)]);
            % plot map
            imagesc(m2,'alphadata',~isnan(m2)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax1,turbo)   
            ax1.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') max(m3(:),[],'omitnan') 1])];
            ax1.XTick = [];
            ax1.YTick = [];

            % additional plots
            ps = linspace(0,size(m1,2),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax1.YLim,'w:','LineWidth',1)

            % text
            % text(0,1,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','k','VerticalAlignment','bottom','rotation',0)
            % text(0,0,sprintf('%s',uci),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','k','VerticalAlignment','bottom','rotation',90,'Interpreter','none')
            text(0,0,sprintf('%.1f',ax1.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','w','VerticalAlignment','bottom')
            text(1,-0.02,sprintf('Cell %d',idx(ii)),'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color','k','VerticalAlignment','top')

            if x(jj)==x(1)
                text(-0.02,0.5,maze_names{2},'Units','normalized','HorizontalAlignment','center','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

        % arena 1 ratemap
        ax2 = axes('Units','pixels','Position',[x(jj) y(jj)+ax1.Position(4)+ysiz(2)+yoffset ax1.Position(3) ax1.Position(4)]);
            if count==0
                ah = add_panel_title('B',sprintf(''),'yoffset',-15,'xoffset',20,'width',400,'fontsize',fs); 
            end
        
            % plot map
            imagesc(m1,'alphadata',~isnan(m1)); hold on;

            % axis settings
            axis xy on
            daspect([1 1 1])
            colormap(ax2,ax1.Colormap)   
            ax2.CLim = ax1.CLim;
            ax2.XTick = [];
            ax2.YTick = [];

            if x(jj)==x(1)            
                text(-0.02,0.5,maze_names{1},'Units','normalized','HorizontalAlignment','center','FontSize',7,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

            % text
            % text(0,1,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','k','VerticalAlignment','bottom','rotation',0)

        % % arena 2 ratemap
        % ax3 = axes('Units','pixels','Position',[x(jj)+ax1.Position(3)-(ax1.Position(3).*0.48) y(jj)+ax1.Position(4)+ysiz(2) ax1.Position(3).*0.48 ax1.Position(4).*0.48]);
        %     % plot map
        %     imagesc(m3,'alphadata',~isnan(m3)); hold on;
        % 
        %     % axis settings
        %     axis xy on
        %     daspect([1 1 1])
        %     colormap(ax3,ax1.Colormap)   
        %     ax3.CLim = ax1.CLim;
        %     ax3.XTick = [];
        %     ax3.YTick = [];

            % text
            % text(0,1,sprintf('%s',maze_names{3}),'Units','normalized','HorizontalAlignment','left','FontSize',6,'Color','k','VerticalAlignment','bottom','rotation',0)
% keyboard
        % count how many plots we have done
        jj = jj+1;
        count = count+1;
        if count==max_plots_now % max we can fit on a page
            disp(sprintf('...saving'))     
            if 1
                fname = [config.fig_dir '\' sprintf('Fig S10_%d.png',fig_count)]; 
                exportgraphics(gcf,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',res);  
                clf(gcf);   
                count = 0;
                fig_count = fig_count+1;
                jj = 1;
                max_plots_now = length(xvec).*length(yvec);   
                [x,y] = ndgrid(xvec,yvec);
            end
            disp(sprintf('Fig %d...',fig_count))    
        end
    end
% return
    close(gcf);

 


















