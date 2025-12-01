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
    fs = [15 10];

    mapset.binsize = 32; % (mm) firing rate map bin siz

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example cells
    xnow = 20;
    ynow = 780;
    xsiz = [90, 10]; % size, buffer
    ysiz = [70, 10]; % size, buffer
    xvec = xnow : sum(xsiz) : 1000;
    xvec = xvec(1:6);
    yvec = ynow : -sum(ysiz).*1.9 : -1000;
    yvec = yvec(1);
    [x,y] = ndgrid(xvec,yvec);
    ar_siz = [xsiz*0.48 ysiz*0.5];

    % find which cells to plot  
    idx = find(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5 & clumaa.planar_spatial_info_shuffles(:,2)>30); % place cells in the hills
    [~,sidx] = sort(clumaa.repetition_score(idx,1),'descend');
    idx = idx(sidx);
    % idx = idx([22 29 2 17 16 35]);
    idx = idx([21 29 2 17 16 36]);

    for ii = 1:numel(x)
        uci = clumaa.uci{idx(ii)};

        %% plot arena (spikes and positions)
        ax1 = axes('Units','pixels','Position',[x(ii) y(ii) xsiz(1) ysiz(1)]);
            if ii==1
                ah = add_panel_title('A',sprintf(''),'yoffset',0,'xoffset',20,'width',400,'fontsize',fs);                
                text(0,0.97,maze_names{2},'Units','normalized','FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','left');
            end    

            % get data
            part_to_plot = 1;
            idxn = find( ismember(clumaa.uci,uci) & clumaa.partn==part_to_plot );
            pos_idx = clumaa.pos_idx(idxn);
            pos = posdata.pos{posdata.pos_idx==pos_idx};
            mf = posdata.maze_frame{posdata.pos_idx==pos_idx}; 

            session_times = posdata.session_times{posdata.pos_idx==pos_idx};
            pidxn = pos.pot > session_times(part_to_plot,1) & pos.pot < session_times(part_to_plot,2); % index for position data in this part

            spt = clumaa.spike_times_s{idxn};
            sidx = spt > session_times(part_to_plot,1) & spt < session_times(part_to_plot,2); % index for spike data in this part, first half
            spike_index = clumaa.spike_index{idxn};

            % plot data
            msiz = 4;
            lwidth = 0.5;
            pos_plot = plot(pos.pox_planar(pidxn),pos.poy_planar(pidxn),'Color',[.5 .5 .5],'LineWidth',lwidth); hold on;
            spk_plot = plot(pos.pox_planar(spike_index(sidx)),pos.poy_planar(spike_index(sidx)),'Marker','.','Color','r','MarkerSize',msiz,'LineStyle','none');
            plot(mf(:,1),mf(:,2),'Color','w');

            % axis settings
            daspect([1 1 1])
            axis off xy
            ax1.XLim = [min(mf(:,1)) max(mf(:,1))];
            ax1.YLim = [min(mf(:,2)) max(mf(:,2))];

            % text
            text(0,1.1,sprintf('Cell %d',ii),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom')
            if ii==1
                text(-0.02,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)

                % legend
                [~,leg] = legendflex([pos_plot spk_plot],{'Positions','Spikes'},'anchor',{'e','e'},'ncol',2,'box','off','buffer',[90,55],'xscale',.5,'fontsize',9); 
                leg(3).LineWidth = 2;
                leg(6).MarkerSize = 18;
                leg(6).Color = spk_plot.Color;
                leg(3).Color = pos_plot.Color;

                % colorbar
                axc = axes('Units','pixels','Position',[ax1.Position(1)+250 ax1.Position(2)+85 100 8]);
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

        %% plot hills (spikes and positions)
        ax2 = axes('Units','pixels','Position',[x(ii) ax1.Position(2)-ax1.Position(4)+20 xsiz(1) ysiz(1)]);
            % get data
            part_to_plot = 2;
            idxn = find( ismember(clumaa.uci,uci) & clumaa.partn==part_to_plot );
            pos_idx = clumaa.pos_idx(idxn);
            pos = posdata.pos{posdata.pos_idx==pos_idx};
            mf = posdata.maze_frame{posdata.pos_idx==pos_idx}; 

            session_times = posdata.session_times{posdata.pos_idx==pos_idx};
            pidxn = pos.pot > session_times(part_to_plot,1) & pos.pot < session_times(part_to_plot,2); % index for position data in this part

            spt = clumaa.spike_times_s{idxn};
            sidx = spt > session_times(part_to_plot,1) & spt < session_times(part_to_plot,2); % index for spike data in this part, first half
            spike_index = clumaa.spike_index{idxn};
            
            % plot data
            msiz = 4;
            lwidth = 0.5;
            pos_plot = plot(pos.pox_planar(pidxn),pos.poy_planar(pidxn),'Color',[.5 .5 .5],'LineWidth',lwidth); hold on;
            spk_plot = plot(pos.pox_planar(spike_index(sidx)),pos.poy_planar(spike_index(sidx)),'Marker','.','Color','r','MarkerSize',msiz,'LineStyle','none');
            plot(mf(:,1),mf(:,2),'Color','k');

            % axis settings
            daspect([1 1 1])
            axis off xy
            ax2.XLim = [min(mf(:,1)) max(mf(:,1))];
            ax2.YLim = [min(mf(:,2)) max(mf(:,2))];

            % additional plots
            ps = linspace(min(mf(:,1)),max(mf(:,1)),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax2.YLim,'k:')            

            % text
            if ii==1
                text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

        %% plot arena (ratemap)
        m2 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
        m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==1)};
        m3 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==3)};

        ax3 = axes('Units','pixels','Position',[x(ii) ax2.Position(2)-ax2.Position(4)+20 xsiz(1) ysiz(1)]);
            imagesc(m1,'alphadata',~isnan(m1)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax3,turbo)   
            ax3.CLim = [0 max([max(m1(:),[],'omitnan') max(m2(:),[],'omitnan') 1])];
            ax3.XTick = [];
            ax3.YTick = [];
            if ii==1
                text(0,0.5,sprintf('%s',maze_names{1}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

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
            text(0,-0.02,sprintf('%.1f',ax4.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','k','VerticalAlignment','top')
            if ii==1
                text(0,0.5,sprintf('%s',maze_names{2}),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end
    end
% keyboard
% %% >>>>>>>>>> autocorrelation schematic
%     xnow = 20;
%     ynow = ynow-250;
%     xbuff = 15;
%     xsiz = 90;
%     ysiz = 70;
%     amap_cmap = cmocean('thermal');
%     amap_clim = [-0.2 1];
% 
%     idx = find(pidx & clumaa.partn==2 & clumaa.f_rate_hz>0.5 & clumaa.planar_spatial_info_shuffles(:,2)>30); % place cells in the hills
%     [~,sidx] = sort(clumaa.repetition_score(idx,1),'descend');
%     idx = idx(sidx);
%     uci = clumaa.uci{idx(1)};
%     idxn = find( ismember(clumaa.uci,uci) & clumaa.partn==2 );
% 
%     % ratemap
%     ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
%         ah = add_panel_title('b',sprintf(''),'yoffset',-15,'xoffset',20,'width',400,'fontsize',[18 13]);                
% 
%         % get ratemap
%         m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
% 
%         % plot data
%         imagesc(m1,'alphadata',~isnan(m1)); hold on;
% 
%         % axis settings
%         axis xy on
%         daspect([1 1 1])
%         colormap(ax1,turbo)   
%         ax1.CLim = [0 max([max(m1(:),[],'omitnan') 1])];
%         ax1.XTick = [];
%         ax1.YTick = [];
%         text(0,1,sprintf('Ratemap'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
% 
%         % additional plots
%         ps = linspace(0,size(m1,2),7);
%         tops = ps(2:2:end);
%         plot([tops; tops],ax1.YLim,'w:')

    % % ratemap
    % ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow xsiz.*1.3 ysiz]);
    %     % get ratemap
    %     m1 = clumaa.ratemap_planar{(ismember(clumaa.uci,uci) & clumaa.partn==2)};
    % 
    %     % plot data
    %     imagesc(m1,'XData',[1 size(m1,2)],'YData',[1 size(m1,1)],'alphadata',ones(size(m1)).*0.5); hold on;
    %     imagesc(m1,'XData',[1 size(m1,2)]+32,'YData',[1 size(m1,1)],'alphadata',ones(size(m1)).*0.5); hold on;
    % 
    %     % axis settings
    %     axis xy on tight
    %     daspect([1 1 1])
    %     colormap(ax2,turbo)   
    % 
    %     ax2.CLim = [0 max([max(m1(:),[],'omitnan') 1])];
    %     ax2.XTick = [];
    %     ax2.YTick = [];
    %     text(0,1,sprintf('Autocorrelation procedure'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
    %     r = corr(reshape(m1(:,32:end),[],1),reshape(m1(:,1:end-31),[],1),'rows','pairwise','type','Pearson');
    %     text(0,0,sprintf('Lag = 1m, r = %.1f',r),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','top','rotation',0)
    % 
    %     % % additional plots
    %     arrow([1 23],[33 23],'Length',8,'TipAngle',30)

    % % autocorr central component
    % ax3 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ynow+30 xsiz 15]);
    %     a1 = clumaa.planar_amap{(ismember(clumaa.uci,uci) & clumaa.partn==2)};    
    %     rindx = ceil(size(a1,1)./2);
    %     a3 = repmat(a1(rindx,:),size(a1,1),1);
    % 
    %     % plot data
    %     xd = ([1 size(a3,2)]-size(a3,2)/2).*32./1000;
    %     imagesc(a3,'XData',xd); hold on;
    % 
    %     % axis settings
    %     axis xy on
    %     box off
    %     % daspect([1 1 1])
    %     colormap(ax3,amap_cmap)   
    %     ax3.CLim = amap_clim;
    %     ax3.YTick = [];
    %     ax3.YColor = 'none';
    %     ax3.XTick = -10:1:10;
    %     xlabel('Lag (m)')
    %     ax3.YTick = [];
    %     ax3.XTickLabelRotation = 0;
    %     text(0,1,sprintf('Autocorrelation'),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
    % 
    % % vertical colormap
    % axc = axes('Units','pixels','Position',[ax3.Position(1)+ax3.Position(3)+5 ynow+10 8 50]);
    %     x = linspace(ax3.CLim(1),ax3.CLim(2),100)';
    %     imagesc(x,'YData',ax3.CLim,'XData',[0 1]);
    %     colormap(axc,ax3.Colormap);
    %     axis on
    %     axc.XTick = [];
    %     axc.YTick = [];
    %     axc.FontSize = 7;
    %     text(0,1,sprintf('%.1f',ax3.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')
    %     text(0,0,sprintf('%.1f',ax3.CLim(1)),'HorizontalAlignment','left','VerticalAlignment','top','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> repetition score
    xnow = 50;
    ynow = ynow-320;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow 110 125]);
        ah = add_panel_title('B',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'repetition_score';
        v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data
        v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % hills data
        v3 = clumaa.(var)(pidx & clumaa.partn==3,1); % arena 2 data
        [ds,gs] = vectorDATAGROUP([],v1(:),v2(:),v3(:)); % linearise data

        % main plot
        meanplot(ds,gs,'linecolor',{'k'},'meancolor',{'k'},'meansize',8,'meanlinewidth',1,'dots',1,'dotcolor',plot_set(1,1:3),'dotsize',20,'dotalpha',0.5,'dotsigma',0.1);

        % axis settings
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

%% >>>>>>>>>> Repetition vs behavioural anisotropy
    %% count repeating cells per rat
    clumaa.rat = cell(size(clumaa,1),1);
    clumaa.ratn = NaN(size(clumaa,1),1);    
    for ii = 1:size(clumaa,1)
        uci = clumaa.uci{ii};
        info = strsplit(uci,'_');
        clumaa.rat{ii} = info{1};
        clumaa.ratn(ii) = str2double( info{1}(3:end) );
    end
    [tmp,~,r] = unique(clumaa.ratn(pidx & clumaa.partn==2,1),'rows');
    % tmp
    x = clumaa.repetition_score(pidx & clumaa.partn==2,1)>prctile(v1,99);

    v = accumarray(r,x,[5,1],@sum,NaN);
    vc = accumarray(r,ones(size(x)),[5,1],@sum,NaN);
    vprop = v ./ vc .* 100;

    %% behavioural anisotropy per rat    
    i = knnsearch(clumaa.pos_idx,posdata.pos_idx,'K',1);
    r = clumaa.rat(i);
    posdata.rat = r;
    posdata.ratn = NaN(size(posdata,1),1);
    for ii = 1:size(posdata,1)
        posdata.ratn(ii) = str2double( posdata.rat{ii}(3:end) );
    end
    [tmp,~,r] = unique(posdata.ratn,'rows');
    % tmp
    x = double( posdata.anisotropy_score(:,2) );

    m = accumarray(r,x,[5,1],@mean,NaN);
    s = accumarray(r,x,[5,1],@nansem,NaN);

    % create axis
    xnow = xnow+170;
    ynow = ynow;
    cols = winter(5);

    ax = axes('Units','pixels','Position',[xnow ynow 125 125]);
        ax.XLim = [-2 40];
        ax.YLim = [-0.1 0];
        ylabel(sprintf('XY-axis behaviour bias'))
        ax.XColor = 'none';

    ax = axes('Units','pixels','Position',[xnow+10 ynow 110 125]);
        ah = add_panel_title('C',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        scatter(vprop(:),m,30,cols,'filled','Marker','o'); hold on;
        for jj = 1:length(tmp)
            ee = errorbar(vprop(jj),m(jj),s(jj),'vertical','Color',cols(jj,:));
            text(vprop(jj),m(jj)+s(jj)+0.003,num2str(jj),'HorizontalAl','center','VerticalAl','bottom','Color',cols(jj,:))
        end
        ax.XLim = [0 40];
        ax.YLim = [-0.1 0];

        xlabel(sprintf('Repetition score >%.2f (%%)',prctile(v1,99)))
        ax.YColor = 'none';

        refline();

        % correlate % cells in panel c with behavioural bias in panel d
        [r,p] = corr(vprop(:),m(:),'rows','pairwise','type','Spearman');
        text(1,1,sprintf('{\\itr} = %.2f, {\\itp} = %.2f',r,p),'FontSize',7,'HorizontalAl','right','VerticalAl','bottom','Units','normalized')

        % % plot data
        % cols = winter(5);
        % for jj = 1:5
        %     s = stem(jj,vprop(jj),'filled','Color',cols(jj,:),'Marker','o'); hold on;
        % end
        % 
        % % axis settings
        % ax.XLim = [0.5 5.5];
        % ax.XTick = 1:5;
        % box off
        % ylabel(sprintf('Repetition score >%.2f (%%)',prctile(v1,99)))
        % xlabel('Rat')

% %% >>>>>>>>>> behavioural anisotropy per rat
%     xnow = xnow+120;
%     ynow = ynow;
% 
%     i = knnsearch(clumaa.pos_idx,posdata.pos_idx,'K',1);
%     r = clumaa.rat(i);
%     posdata.rat = r;
%     posdata.ratn = NaN(size(posdata,1),1);
%     for ii = 1:size(posdata,1)
%         posdata.ratn(ii) = str2double( posdata.rat{ii}(3:end) );
%     end
%     [~,~,r] = unique(posdata.ratn,'rows');
%     x = double( posdata.anisotropy_score(:,2) );
% 
%     m = accumarray(r,x,[5,1],@mean,NaN);
%     s = accumarray(r,x,[5,1],@nansem,NaN);
% 
%     % create axis
%     ax = axes('Units','pixels','Position',[xnow ynow 70 125]);
%         ah = add_panel_title('d',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',[18 13]);                
% 
%         % plot data
%         for jj = 1:5
%             e = errorbar(jj,m(jj),s(jj),'Color',cols(jj,:),'Marker','.','LineStyle','none','MarkerSize',15); hold on;
%         end
% 
%         % axis settings
%         ax.XLim = [0.5 5.5];
%         ax.XTick = 1:5;
%         box off
%         ylabel(sprintf('XY-axis behaviour bias'))
%         xlabel('Rat')
%         ax.YLim = [-0.1 0];
% 
%         % text
%         text(1.05,1.05,sprintf('No bias'),'FontSize',7,'HorizontalAlignment','right','Units','normalized','VerticalAlignment','top')
%         text(1.05,0,sprintf('Bias for X'),'FontSize',7,'HorizontalAlignment','right','Units','normalized','VerticalAlignment','bottom')
% 
%         % stats
%         [r,p] = corr(vprop,m,"rows",'pairwise','type','Pearson');
% 
%         % stats
%         axt = axes('Units','pixels','Position',ax.Position,'Color','none');
%             axis off
%             axt.XLim = ax.XLim;
%             axt.YLim = [0 1];
%             [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.4);
%             % N = [sum(~isnan(v1(:))) sum(~isnan(v2a(:))) sum(~isnan(v2b(:))) sum(~isnan(v3(:)))];
% 
%     % correlate % cells in panel c with behavioural bias in panel d
%     % [r,p] = corr(vprop(:),m(:),'rows','pairwise','type','Spearman')
% % keyboard

%% >>>>>>>>>> repetition (autocorrelation plots)
    xnow = 45;
    ynow = ynow-290;
    xbuff = 25;
    ysiz = 200;
    xsiz = 100;
    amap_cmap = cmocean('thermal');
    amap_clim = [-0.2 1];

    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]);
        ah = add_panel_title('D',sprintf(''),'yoffset',-20,'xoffset',-5,'width',400,'fontsize',fs);                
    
        % get data
        v1 = clumaa.amap_cent(pidx & clumaa.partn==1); % arena data
        v1 = cell2mat(v1);
        v1idx = clumaa.repetition_score(pidx & clumaa.partn==1,1); % arena data
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
        v2 = clumaa.amap_cent(pidx & clumaa.partn==2); % hills data
        v2 = cell2mat(v2);
        v2idx = clumaa.repetition_score(pidx & clumaa.partn==2,1); % hills data
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

%% >>>>>>>>>> Example cells
    % find which cells to plot  
    v2 = clumaa.ratemap_planar(pidx & clumaa.partn==2); % hills data
    v2idx = clumaa.repetition_score(pidx & clumaa.partn==2,1); % hills data
    [~,sidx] = sort(v2idx,'descend');
    v2s = v2(sidx,:);

    vals = fliplr([1 61 159 266 361 456]);

    xnow = xnow+250;
    ynow = ynow-30;
    xsiz = 70;
    ysiz = 60;
    yvecn = ynow:ysiz-20:ax1.Position(2)+ax1.Position(4);

    for jj = 1:length(yvecn)
        ax = axes('Units','pixels','Position',[xnow yvecn(jj) xsiz ysiz]);
            m1 = v2s{vals(jj)};
        
            imagesc(m1,'alphadata',~isnan(m1)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax,turbo)   
            ax.CLim = [0 max([max(m1(:),[],'omitnan') 1])];
            ax.XTick = [];
            ax.YTick = [];
            if jj==length(yvecn)
                text(0,1,sprintf('Example cells'),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')
            end
    end

    axl = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3) ax2.Position(2) 100 ax2.Position(4)],'Clipping','off');
        axl.Color = 'none';
        axl.YLim = ax2.YLim;
        axl.XLim = [0 1];
        axis manual
        axis ij off

        plot_x = [0 0.25];
        vs = [vals' [472 375 290 190 100 1]'];
        line(plot_x,vs,'Color','k')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[axl.Position(1)+30 axl.Position(2)-48 50 8]);
        x = linspace(ax.CLim(1),ax.CLim(2),100);
        imagesc(x,'XData',ax.CLim,'YData',[0 1]);
        colormap(axc,ax.Colormap);
        axis on
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.f ',0),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' Max'),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1.8,sprintf('Firing rate (Hz)'),'FontSize',axc.FontSize,'HorizontalAl','center','Units','normalized')

%% >>>>>>>>>> repetition (field plots)
    % repetition scores
    var = 'repetition_score';
    % rs_cutoff = -inf;
    v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 data    
    rs_cutoff = prctile(v1,99);
    r1 = clumaa.(var)(pidx & clumaa.partn==1 & clumaa.repetition_score(:,1)>rs_cutoff,1); % arena 1 data
    r2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.repetition_score(:,1)>rs_cutoff,1); % hills data

    % ratemaps, averaged, zscored
    var = 'ratemap_planar';
    m2 = clumaa.(var)(pidx & clumaa.partn==2 & clumaa.repetition_score(:,1)>rs_cutoff,1); % arena 1 data
    m2 = cellfun(@(x) imresize(x,[48 96],'nearest'),m2,'UniformOutput',0);
    m2 = cat(3,m2{:});   
    m2 = squeeze(sum(m2,1,'omitnan'));
    z2 = zscore(rot90(m2),[],2);

    u = ismember(clumaa.uci,clumaa.uci(pidx & clumaa.partn==2 & clumaa.repetition_score(:,1)>rs_cutoff)); % hills data
    m1 = clumaa.(var)(u & clumaa.partn==1); % arena 1 data
    m1 = cellfun(@(x) imresize(x,[48 96],'nearest'),m1,'UniformOutput',0);
    m1 = cat(3,m1{:});   
    m1 = squeeze(sum(m1,1,'omitnan'));
    z1 = zscore(rot90(m1),[],2);

    % ratemap peak vals
    [~,midx1] = max(z1,[],2,'omitnan');
    [~,idx1] = sort(midx1,'ascend');
    z1 = z1(idx1,:);

    [~,midx2] = max(z2,[],2,'omitnan');
    [~,idx2] = sort(midx2,'ascend');
    z2 = z2(idx2,:);

    % plotting settings
    xnow = xnow+120;
    ynow = yvec(1)-385;
    xbuff = 10;
    ysiz = 200;
    xsiz = 95;

    ax1 = axes('Units','pixels','Position',[xnow ynow xsiz ysiz]); % x-axis plot 
        ah = add_panel_title('E',sprintf(''),'yoffset',-20,'xoffset',-5,'width',400,'fontsize',fs);                
    
        imagesc(z1,'XData',0:size(z1,2).*mapset.binsize./1000); hold on;
        axis ij on
        colormap(ax1,turbo)   
        ax1.CLim = [-1 4];
        ylabel('Ordered cells')
        ax1.YTick = unique([1 0:50:size(z1,1) size(z1,1)]);
        ax1.FontSize = 8;
        xlabel('Position (m)')
    
        % additional plots
        maze_sections = 0:0.5:3;
        hill_x = maze_sections(2:2:end);
        plot([hill_x; hill_x],ax1.YLim,'w:','LineWidth',1)
    
        % text
        text(0.5,1,sprintf('%s',maze_names{1}),'FontSize',8,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','bottom')

    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+xbuff ax1.Position(2) ax1.Position(3) ax1.Position(4)]); % x-axis plot    
        imagesc(z2,'XData',0:size(z2,2).*mapset.binsize./1000); hold on;
        axis ij on
        colormap(ax2,turbo)   
        ax2.CLim = [-1 4];
        xlabel('Position (m)')
        ax2.YTick = [];
        ax2.FontSize = 8;
    
        % additional plots
        maze_sections = 0:0.5:3;
        hill_x = maze_sections(2:2:end);
        plot([hill_x; hill_x],ax2.YLim,'w:','LineWidth',1)

        % text
        text(0.5,1,sprintf('%s',maze_names{2}),'FontSize',8,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','bottom')

    % horizontal colormap
    axc = axes('Units','pixels','Position',[ax2.Position(1)-40 ax2.Position(2)-50 70 8]);
        x = linspace(ax2.CLim(1),ax2.CLim(2),100);
        imagesc(x,'XData',ax2.CLim,'YData',[0 1]);
        colormap(axc,ax2.Colormap);
        axis off
        axc.XTick = [];
        axc.YTick = [];
        axc.FontSize = 7;
        text(0,0.5,sprintf('%.1f ',ax2.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(1,0.5,sprintf(' %.1f',ax2.CLim(2)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
        text(0.5,1,sprintf('Norm. rate (z)'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')

%% >>>>>>>>>> Field distribution along hills
    ynow = ynow-ysiz+50;
    xnow = xnow;

    % create axis
    ax = axes('Units','pixels','Position',[xnow ynow xsiz*2+xbuff 75]);
        ah = add_panel_title('F',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                

        % get data
        v1 = (((midx1-1) .* mapset.binsize) + mapset.binsize/2) ./ 1000; % convert bins to m
        v2 = (((midx2-1) .* mapset.binsize) + mapset.binsize/2) ./ 1000; % convert bins to m
% keyboard
        m = max( [max(v1) max(v2)] );
        v1 = v1 ./ max(v1) .* 3;
        v2 = v2 ./ max(v2) .* 3;

        xi = linspace(0,3,32);
        f1 = histcounts(v1,xi,'Normalization','probability') .* 100;
        f2 = histcounts(v2,xi,'Normalization','probability') .* 100;
        f3 = ones(size(f1)) ./ numel(f1) .* 100; % uniform

        % plot data
        xd = movmean(xi,2,'Endpoints','discard');
        plot(xd,f1,'Color',plot_set{1,1},'LineWidth',1.5,'Marker','.'); hold on
        plot(xd,f2,'Color',plot_set{1,2},'LineWidth',1.5,'Marker','.'); hold on

        ncos = 100;
        xi1 = linspace(-pi,5*pi,100);                
        xi2 = linspace(xi(1),xi(end),100);                        
        a3 = area(xi2,(cos(xi1)+1).*(ax.YLim(2)./2),'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        uistack(a3,'bottom')

        % axis settings
        ax.XLim = [0 3];
        xlabel('Position (m)')
        ylabel('Place cells (%)')
        ytickformat('%.f')
        box off
        ax.FontSize = 8;

        [h,p,k] = kstest2(v1(:),v2(:));
        text(0,1.05,sprintf('D = %.1f, {\\itp} = %.2f',k,p),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)
% keyboard

%% >>>>>>>>>> Field distribution along hills
    ynow = ynow-135;
    xnow = xnow;

    % create axis
    ax = axes('Units','pixels','Position',[ax.Position(1) ynow ax.Position(3) ax.Position(4)]);
        ah = add_panel_title('G',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                

        v1p = wrapTo2Pi(v1.*(2*pi)); % convert to phase and wrap to circle
        v2p = wrapTo2Pi(v2.*(2*pi)); % convert to phase and wrap to circle

        % get density of data
        xi = linspace(0,2*pi,32);
        f1 = histcounts(v1p,xi,'Normalization','probability') .* 100;
        f2 = histcounts(v2p,xi,'Normalization','probability') .* 100;
        f3 = ones(size(f1)) ./ numel(f1) .* 100; % uniform

        % main plot
        xd = movmean(xi,2,'EndPoints','discard')./pi.*0.5-0.5;
        a1 = plot(xd,f1,'Color',plot_set{1,1},'LineWidth',1.5,'Marker','.'); hold on
        a2 = plot(xd,f2,'Color',plot_set{1,2},'LineWidth',1.5,'Marker','.'); hold on

        ncos = 100;
        xi1 = linspace(-pi,pi,100);                
        xi2 = linspace(xd(1),xd(end),100);                        
        a3 = area(xi2,(cos(xi1)+1).*(ax.YLim(2)./2),'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        uistack(a3,'bottom')

        % axis settings
        ax.XLim = [-0.5 0.5]; 
        ax.FontSize = 8;
        box off
        xlabel(sprintf('Position relative to hill peak (m)'))    
        ylabel('Place cells (%)')
        ytickformat('%.f');

        % additional plots
        [~,leg] = legendflex([a1 a2 a3],{sprintf('%s',maze_names{1}),sprintf('%s',maze_names{2}),'Terrain'},'anchor',{'ne','ne'},'ncol',1,'box','off','buffer',[5,45],'xscale',.5,'fontsize',8); 
        leg(8).Children.FaceAlpha = a3.FaceAlpha;
        leg(8).Children.FaceColor = a3.FaceColor;
        leg(8).Children.EdgeColor = a3.EdgeColor;
        leg(5).Color = a1.Color;
        leg(7).Color = a2.Color;
        leg(4).Color = a1.Color;
        leg(6).Color = a2.Color;

        % stats
        % rv1 = circ_r(v1p);
        % rv2 = circ_r(v2p);
        % [ckp,ckf] = circ_ktest(v1p(:),v2p(:));
        [h,p,k] = kstest2(v1p(:),v2p(:));
        % numel(v1p)
        text(0,1.05,sprintf('D = %.1f, {\\itp} = %.2f',k,p),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',0)

        % keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig 4.png']; 
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


 


















