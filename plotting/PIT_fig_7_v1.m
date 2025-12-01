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
%% >>>>>>>>>> Check that correlations and shuffles are complete
    clumaa = PIT_correlations(config,pidx,clumaa,posdata,0);
    fs = [15 10];

    if 1
%% Parse inputs
        % Create figure
        fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
        set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
        set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example cells (low azimuth stability)
        xnow = 30;
        ynow = 720;
        xsiz = 107;
        ysiz = 70;
        xbuff = xsiz+10;
        
        % find which cells to plot  
        [~,sidx] = sort(clumaa.azimuth_stability(pidx & clumaa.partn==2),'ascend');
        idx = sidx([11 16 26 27 30]);
        % ucis = {'RG19_221021_t2_c1','RG11_220520_t2_c5','RG19_221013_t3_c7'};
        ucis_c = unique(clumaa.uci(pidx & clumaa.partn==2),'stable');
        ucis = ucis_c(idx);
    
        xvec = [xnow xnow+xbuff xnow+(xbuff.*2) xnow+(xbuff.*3) xnow+(xbuff.*4)];
        for uu = 1:length(ucis)
            uci = ucis{uu};
            m2a = clumaa.ratemap_azimuth_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),1}; % hills half 1, left to right 
            m2b = clumaa.ratemap_azimuth_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),2}; % hills half 2, right to left
        
            % plot hills half 1
            ax2a = axes('Units','pixels','Position',[xvec(uu) ynow+60 xsiz ysiz]);
                imagesc(m2a,'alphadata',~isnan(m2a)); hold on;
                axis xy on
                daspect([1 1 1])
                colormap(ax2a,turbo)   
                ax2a.CLim = [0 max([max(m2a(:),[],'omitnan') max(m2b(:),[],'omitnan') 0.1])];
                ax2a.XTick = [];
                ax2a.YTick = [];
        
                % additional plots
                ps = linspace(ax2a.XLim(1),ax2a.XLim(2),7);
                tops = ps(2:2:end);
                valleys = ps(1:2:end);
                plot([tops; tops],ax2a.YLim,'w:')
        
                % text
                if uu==1
                    ah = add_panel_title('A',sprintf(''),'yoffset',-5,'xoffset',10,'width',400,'fontsize',fs);                
                    text(-0.02,0.5,sprintf('Right -> left'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
                end
                text(0,0,sprintf('%.1f',ax2a.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
                text(0,1.1,sprintf('Cell %d ({\\itr} = %.2f)',uu,clumaa.azimuth_stability(ismember(clumaa.uci,uci) & clumaa.partn==2)),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom')
                
            % plot hills half 2
            ax2b = axes('Units','pixels','Position',[xvec(uu) ynow xsiz ysiz]);
                imagesc(m2b,'alphadata',~isnan(m2b)); hold on;
                axis xy on
                daspect([1 1 1])
                colormap(ax2b,turbo)   
                ax2b.CLim = ax2a.CLim;
                ax2b.XTick = [];
                ax2b.YTick = [];
        
                % additional plots
                ps = linspace(ax2b.XLim(1),ax2b.XLim(2),7);
                tops = ps(2:2:end);
                plot([tops; tops],ax2b.YLim,'w:')
        
                % text
                if uu==1
                    text(-0.02,0.5,sprintf('Left -> right'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
                end
    
            if uu==1
                % horizontal colormap
                axc = axes('Units','pixels','Position',[xnow+100 ynow-20 70 8]);
                    x = linspace(ax2a.CLim(1),ax2a.CLim(2),100);
                    imagesc(x,'XData',ax2a.CLim,'YData',[0 1]);
                    colormap(ax2a.Colormap);
                    axis on
                    axc.XTick = [];
                    axc.YTick = [];
                    axc.FontSize = 7;
                    text(0,0.5,sprintf('%.1f ',ax2a.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
                    text(1,0.5,sprintf(' Max'),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
                    text(0.5,1,sprintf('Firing rate (Hz)'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')                
            end
        end
% return
%% >>>>>>>>>> Example cells (high repetition score)
        xnow = 30;
        ynow = ynow-170;
        xsiz = 107;
        ysiz = 70;
        xbuff = xsiz+10;
        
        % find which cells to plot  
        [~,sidx] = sort(clumaa.repetition_score(pidx & clumaa.partn==2),'descend');
        idx = sidx([1 2 3 4 5]);
        % ucis = {'RG19_221021_t2_c1','RG11_220520_t2_c5','RG19_221013_t3_c7'};
        ucis_c = unique(clumaa.uci(pidx & clumaa.partn==2),'stable');
        ucis = ucis_c(idx);
    
        xvec = [xnow xnow+xbuff xnow+(xbuff.*2) xnow+(xbuff.*3) xnow+(xbuff.*4)];
        for uu = 1:length(ucis)
            uci = ucis{uu};
            m2a = clumaa.ratemap_azimuth_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),1}; % hills half 1     
            m2b = clumaa.ratemap_azimuth_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),2}; % hills half 2     
        
            % plot hills half 1
            ax2a = axes('Units','pixels','Position',[xvec(uu) ynow+60 xsiz ysiz]);
                imagesc(m2a,'alphadata',~isnan(m2a)); hold on;
                axis xy on
                daspect([1 1 1])
                colormap(ax2a,turbo)   
                ax2a.CLim = [0 max([max(m2a(:),[],'omitnan') max(m2b(:),[],'omitnan') 0.1])];
                ax2a.XTick = [];
                ax2a.YTick = [];
        
                % additional plots
                ps = linspace(ax2a.XLim(1),ax2a.XLim(2),7);
                tops = ps(2:2:end);
                valleys = ps(1:2:end);
                plot([tops; tops],ax2a.YLim,'w:')
        
                % text
                if uu==1
                    ah = add_panel_title('B',sprintf(''),'yoffset',-5,'xoffset',10,'width',400,'fontsize',fs);                
                    text(-0.02,0.5,sprintf('Right -> left'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
                end
                text(0,0,sprintf('%.1f',ax2a.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
                text(0,1.1,sprintf('Cell %d ({\\itr} = %.2f)',uu,clumaa.azimuth_stability(ismember(clumaa.uci,uci) & clumaa.partn==2)),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom')
                
            % plot hills half 2
            ax2b = axes('Units','pixels','Position',[xvec(uu) ynow xsiz ysiz]);
                imagesc(m2b,'alphadata',~isnan(m2b)); hold on;
                axis xy on
                daspect([1 1 1])
                colormap(ax2b,turbo)   
                ax2b.CLim = ax2a.CLim;
                ax2b.XTick = [];
                ax2b.YTick = [];
        
                % additional plots
                ps = linspace(ax2b.XLim(1),ax2b.XLim(2),7);
                tops = ps(2:2:end);
                plot([tops; tops],ax2b.YLim,'w:')
        
                % text
                if uu==1
                    text(-0.02,0.5,sprintf('Left -> right'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
                end
    
        end

%% >>>>>>>>>> Azimuth correlation results
        xnow = xnow+20;
        ynow = ynow-140;
        fname = [config.data_out_dir 'PIT_correlations_azimuth_shuffles.mat'];
        load(fname,'shuffs_azimuth');
    
        % create axis
        ax2 = axes('Units','pixels','Position',[xnow ynow 220 110]);  
            ah = add_panel_title('C',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
        
            pp = 2;
            var = 'within_session_stability';
            v2a = clumaa.(var)(pidx & clumaa.partn==pp,1); % within session stability,
            var = 'azimuth_stability';
            v2b = clumaa.(var)(pidx & clumaa.partn==pp,1); % pitch stability,
            v2s = shuffs_azimuth(:,2); % pitch stability shuffle,
            r = clumaa.repetition_score(pidx & clumaa.partn==1,1); % repetition score, arena 1
            cutoff = prctile(r,99);
            v3 = v2a(clumaa.repetition_score(pidx & clumaa.partn==pp,1)>cutoff);

            xi = -1:0.01:1;
            bw = 0.05;
            f1 = ksdensity(v2a(:),xi(:),"Bandwidth",bw);
            f2 = ksdensity(v2b(:),xi(:),"Bandwidth",bw);
            f3 = ksdensity(v2s(:),xi(:),"Bandwidth",bw);
            f4 = ksdensity(v3(:),xi(:),"Bandwidth",bw);

            % main plot    
            alph = 0.7;          
            a1 = area(xi,f1,'FaceColor',plot_set{1,2},'EdgeColor','k','FaceAlpha',alph); hold on;
            a2 = area(xi,f2,'FaceColor','m','EdgeColor','k','FaceAlpha',alph);       
            a3 = area(xi,f3,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--'); hold on;       
            a4 = area(xi,f4,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','-'); hold on;       

            % axis settings
            ax2.XLim = [-0.5 1]; 
            ax2.FontSize = 8;
            ytickformat('%.1f')
            box off
            xlabel(sprintf('Correlation ({\\itr})'))    
            ylabel(sprintf('PDF'))  

            [~,leg] = legendflex([a1 a2 a3 a4],...
                {sprintf('Session half 1 vs 2'),...
                sprintf('Leftward vs rightward'),...
                sprintf('Leftward vs rightward shuffle'),sprintf('Leftward vs rightward (terrain cells only)')},...
                'anchor',{'sw','sw'},'ncol',1,'box','off','buffer',[-10,-100],'xscale',.5,'fontsize',9); 
            
            leg(5).Children.FaceAlpha = a1.FaceAlpha;
            leg(6).Children.FaceAlpha = a2.FaceAlpha;
            leg(7).Children.FaceAlpha = a3.FaceAlpha;
            leg(8).Children.FaceAlpha = a4.FaceAlpha;  

            leg(5).Children.FaceColor = a1.FaceColor;
            leg(6).Children.FaceColor = a2.FaceColor;
            leg(7).Children.FaceColor = a3.FaceColor;
            leg(8).Children.FaceColor = a4.FaceColor;


%% >>>>>>>>>> Between session correlation results (Cohen's d values)
        % create axis
        ax = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+50 ynow 50 100]);
            % ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',-20,'width',400);
            e1 = computeEffectSize(v2a,v2b); % within vs pitch
            g1 = e1.cliff_delta;
            e3 = computeEffectSize(v2a,v3); % within vs pitch (terrain only)
            g3 = e3.cliff_delta;            
            e2 = computeEffectSize(v2b,v2s); % pitch vs shuffle
            g2 = e2.cliff_delta;
    
            stem(1:3,[g1 g3 g2],'filled','ko'); hold on;
    
            % axis settings
            ax.YLim = [-0.1 1];
            ax.XLim = [0.5 3.5];
            ylabel(sprintf('Cliff''s Delta'))   
            ax.XTick = 1:3;
            ax.XTickLabel = {};
            ytickformat('%.1f');
            ax.FontSize = 8;
            box off
            ax.XColor = 'none';
    
            % additional plots
            line(ax.XLim,[0 0],'Color',[.5 .5 .5]);
    
            % stats
            vcoef = 1;
            fsiz = 12;
            dat = {v2a v2b;v2a v3;v2b v2s};
            for ii = 1:3
                [z,p,~,~,~] = permutationZTest(dat{ii,1},dat{ii,2},'method','cliff');
                if p<.001
                    text(ii,ax.YLim(2).*vcoef,'***','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
                elseif p<0.01
                    text(ii,ax.YLim(2).*vcoef,'**','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
                elseif p<.05
                    text(ii,ax.YLim(2).*vcoef,'*','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
                end
            end

            axi = axes('Units','pixels','Position',[ax.Position(1) ax.Position(2)-ax.Position(3)+10 ax.Position(3) ax.Position(3)]);
                axi.Color = 'none';
                axi.XLim = ax.XLim;
                axi.YLim = ax.XLim;
                rsiz = 0.6;
                rectangle('Position',[1-rsiz/2 1.75 rsiz rsiz],'FaceColor',[a1.FaceColor a1.FaceAlpha],'EdgeColor',a1.EdgeColor,'LineStyle',a1.LineStyle);
                rectangle('Position',[1-rsiz/2 0.6 rsiz rsiz],'FaceColor',[a2.FaceColor a2.FaceAlpha],'EdgeColor',a2.EdgeColor,'LineStyle',a2.LineStyle);
                ty = 1.1;
                text(1,ty,'vs','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize)

                rectangle('Position',[2-rsiz/2 1.75 rsiz rsiz],'FaceColor',[a1.FaceColor a1.FaceAlpha],'EdgeColor',a1.EdgeColor,'LineStyle',a1.LineStyle);
                rectangle('Position',[2-rsiz/2 0.6 rsiz rsiz],'FaceColor',[a4.FaceColor a4.FaceAlpha],'EdgeColor',a4.EdgeColor,'LineStyle',a4.LineStyle);
                text(2,ty,'vs','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize)

                rectangle('Position',[3-rsiz/2 1.75 rsiz rsiz],'FaceColor',[a2.FaceColor a2.FaceAlpha],'EdgeColor',a2.EdgeColor,'LineStyle',a2.LineStyle);
                rectangle('Position',[3-rsiz/2 0.6 rsiz rsiz],'FaceColor',[a3.FaceColor a3.FaceAlpha],'EdgeColor',a3.EdgeColor,'LineStyle',a3.LineStyle);
                text(3,ty,'vs','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize)

                axis off

%% >>>>>>>>>> Relationship azimuth stability repetition score
        xnow = xnow+380;
        ynow = ynow;
    
        % create axis
        ax2 = axes('Units','pixels','Position',[xnow ynow 160 110]);  
            ah = add_panel_title('D',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
        
            var = 'azimuth_stability';
            v1 = clumaa.(var)(pidx & clumaa.partn==2,1); % pitch stability, hills
            var = 'repetition_score';
            v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % repetition score, hills
    
            % main plot
            v = [v2(:) v1(:)];
            d = NaN(length(v(:,1)),1);
            v_norm = v ./ max(v,[],1,'omitmissing');
            sigma = 0.04;
            for j = 1:size(v_norm,1) % loop over every point
                D = sqrt(sum(bsxfun(@minus,v_norm(j,:),v_norm).^2,2));
                d(j) = sum( normpdf(D,0,sigma),'all','omitmissing');
            end            
            scatter(v2,v1,30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5); hold on;
    
            % axis settings
            xlabel(sprintf('Repetition score'))            
            ylabel(sprintf('Directional stability (r)'))    
            % ax.XTick = [0.001 0.01 0.1 1 10 100];
            % ax.YTick = [0.001 0.01 0.1 1 10 100];        
            ax2.XLim = [-0.5 1.5]; 
            ax2.YLim = [0 1]; 
            ax2.CLim(1) = 0;
            ax2.FontSize = 8;
            % ax.YScale = 'log';
            % ax.XScale = 'log';
            ytickformat('%.1f')
            xtickformat('%.1f')        
            rf = refline;
            set(rf,'Color','k');
            grid on
    
            % % additional plots
            % line([ax.XLim(1) 100],[ax.XLim(1) 100],'Color',[.5 .5 .5],'LineStyle','--')
            % 
            % stats
            [r,p] = corr(v2,v1,'rows','pairwise','type','Spearman');
            text(1,0,sprintf('r = %.2f',r),'FontSize',8,'HorizontalAlignment','right','Units','normalized','VerticalAlignment','bottom','Color','k')
            % [r,p] = corr(fs(:,1),fs(:,3),'rows','pairwise','type','Pearson');    
            % text(0,1.02,sprintf('Arena 1 vs Arena 2: r = %.1f',r),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','#012030')

            % horizontal colormap
            axc = axes('Units','pixels','Position',[xnow ynow-60 100 10]);
                x = linspace(ax2.CLim(1),ax2.CLim(2),100);
                imagesc(x,'XData',ax2.CLim,'YData',[0 1]);
                colormap(ax2.Colormap);
                axis on
                axc.XTick = [];
                axc.YTick = [];
                axc.FontSize = 7;
                text(0,0.5,sprintf('%.1f ',0),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
                text(1,0.5,sprintf(' Max'),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
                text(0.5,1,sprintf('Density'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized') 

% % stats comparing directionality in terrain cells vs non terrain cells
% var = 'repetition_score';
% r = clumaa.(var)(pidx & clumaa.partn==1,1); % repetition score, hills
% cutoff = prctile(r,99);
% [ds,gs] = vectorDATAGROUP([],v1(v2<cutoff),v1(v2>cutoff));
% [p,a,s] = anova1(ds,gs,'off')
                keyboard

%% >>>>>>>>>> Save the overall figure
        if 1
            fname = [config.fig_dir '\Fig S8.png']; 
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
    end


%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 210.*3 297.*3],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Example cells (low azimuth stability)
    xnow = 40;
    ynow = 720;
    xsiz = 107;
    ysiz = 70;
    xbuff = xsiz+10;
    
    % find which cells to plot  
    [~,sidx] = sort(clumaa.pitch_stability(pidx & clumaa.partn==2),'ascend');
    idx = sidx([7 12 15 17 22]);
    ucis_c = unique(clumaa.uci(pidx & clumaa.partn==2),'stable');
    ucis = ucis_c(idx);

    xvec = [xnow xnow+xbuff xnow+(xbuff.*2) xnow+(xbuff.*3) xnow+(xbuff.*4)];
    for uu = 1:length(ucis)
        uci = ucis{uu};
        m2a = clumaa.ratemap_pitch_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),1}; % hills half 1, pitch down   
        m2b = clumaa.ratemap_pitch_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),2}; % hills half 2, pitch up    
        mpit = clumaa.hd_3d_ratemap{ismember(clumaa.uci,uci) & clumaa.partn==2};
    
        % plot hills half 1
        ax2a = axes('Units','pixels','Position',[xvec(uu) ynow+60 xsiz ysiz]);
            imagesc(m2a,'alphadata',~isnan(m2a)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax2a,turbo)   
            ax2a.CLim = [0 max([max(m2a(:),[],'omitnan') max(m2b(:),[],'omitnan') 0.1])];
            ax2a.XTick = [];
            ax2a.YTick = [];
    
            % additional plots
            ps = linspace(ax2a.XLim(1),ax2a.XLim(2),7);
            tops = ps(2:2:end);
            valleys = ps(1:2:end);
            plot([tops; tops],ax2a.YLim,'w:')
    
            % text
            if uu==1
                ah = add_panel_title('A',sprintf(''),'yoffset',-10,'xoffset',-5,'width',400,'fontsize',fs);                
                text(-0.02,0.5,sprintf('Nose down'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end
            text(0,0,sprintf('%.1f',ax2a.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
            text(0,1.1,sprintf('Cell %d ({\\itr} = %.2f)',uu,clumaa.pitch_stability(ismember(clumaa.uci,uci) & clumaa.partn==2)),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom')
            
        % plot hills half 2
        ax2b = axes('Units','pixels','Position',[xvec(uu) ynow xsiz ysiz]);
            imagesc(m2b,'alphadata',~isnan(m2b)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax2b,turbo)   
            ax2b.CLim = ax2a.CLim;
            ax2b.XTick = [];
            ax2b.YTick = [];
    
            % additional plots
            ps = linspace(ax2b.XLim(1),ax2b.XLim(2),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax2b.YLim,'w:')
    
            % text
            if uu==1
                text(-0.02,0.5,sprintf('Nose up'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

        % plot pitch tuning
        ax3b = axes('Units','pixels','Position',[xvec(uu) ynow-70 xsiz ysiz]);
            % plot data
            imagesc([-180 180],[-90 90],mpit,'alphadata',~isnan(mpit))
    
            % axis settings
            ax3b.XTick = -180:90:180;
            ax3b.YTick = -90:45:90;            
            if uu==1
                xlabel(sprintf('Azimuth (%c)',176))
                ylabel(sprintf('Tilt (%c)',176))                
            else
                xlabel(sprintf('Azimuth (%c)',176))
                ax3b.YTickLabel = {};                
            end
            axis xy tight
            view(0,90);
            ax3b.FontSize = 8;
            ax3b.CLim = [0,max(mpit(:))];
            colormap(gca,turbo)
            % ax3b.XTickLabelRotation = 0;

        if uu==1
            % horizontal colormap
            axc = axes('Units','pixels','Position',[xnow+150 ynow+145 70 8]);
                x = linspace(ax2a.CLim(1),ax2a.CLim(2),100);
                imagesc(x,'XData',ax2a.CLim,'YData',[0 1]);
                colormap(ax2a.Colormap);
                axis on
                axc.XTick = [];
                axc.YTick = [];
                axc.FontSize = 7;
                text(0,0.5,sprintf('%.1f ',ax2a.CLim(1)),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
                text(1,0.5,sprintf(' Max'),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
                text(0.5,1,sprintf('Firing rate (Hz)'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized')                
        end
    end

%% >>>>>>>>>> Example cells (high repetition score)
    xnow = 40;
    ynow = ynow-270;
    xsiz = 107;
    ysiz = 70;
    xbuff = xsiz+10;
    
    % find which cells to plot  
    [~,sidx] = sort(clumaa.repetition_score(pidx & clumaa.partn==2),'descend');
    idx = sidx([1 2 3 4 5]);
    ucis_c = unique(clumaa.uci(pidx & clumaa.partn==2),'stable');
    ucis = ucis_c(idx);

    xvec = [xnow xnow+xbuff xnow+(xbuff.*2) xnow+(xbuff.*3) xnow+(xbuff.*4)];
    for uu = 1:length(ucis)
        uci = ucis{uu};
        m2a = clumaa.ratemap_pitch_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),1}; % hills half 1, pitch down 
        m2b = clumaa.ratemap_pitch_half{(ismember(clumaa.uci,uci) & clumaa.partn==2),2}; % hills half 2, pitch up  
        mpit = clumaa.hd_3d_ratemap{ismember(clumaa.uci,uci) & clumaa.partn==2};

        % plot hills half 1
        ax2a = axes('Units','pixels','Position',[xvec(uu) ynow+60 xsiz ysiz]);
            imagesc(m2a,'alphadata',~isnan(m2a)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax2a,turbo)   
            ax2a.CLim = [0 max([max(m2a(:),[],'omitnan') max(m2b(:),[],'omitnan') 0.1])];
            ax2a.XTick = [];
            ax2a.YTick = [];
    
            % additional plots
            ps = linspace(ax2a.XLim(1),ax2a.XLim(2),7);
            tops = ps(2:2:end);
            valleys = ps(1:2:end);
            plot([tops; tops],ax2a.YLim,'w:')
    
            % text
            if uu==1
                ah = add_panel_title('B',sprintf(''),'yoffset',-10,'xoffset',-5,'width',400,'fontsize',fs);                
                text(-0.02,0.5,sprintf('Nose down'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end
            text(0,0,sprintf('%.1f',ax2a.CLim(2)),'Units','normalized','HorizontalAlignment','left','FontSize',7,'Color','w','VerticalAlignment','bottom')
            text(0,1.1,sprintf('Cell %d ({\\itr} = %.2f)',uu,clumaa.pitch_stability(ismember(clumaa.uci,uci) & clumaa.partn==2)),'Units','normalized','HorizontalAlignment','left','FontSize',8,'Color','k','VerticalAlignment','bottom')
            
        % plot hills half 2
        ax2b = axes('Units','pixels','Position',[xvec(uu) ynow xsiz ysiz]);
            imagesc(m2b,'alphadata',~isnan(m2b)); hold on;
            axis xy on
            daspect([1 1 1])
            colormap(ax2b,turbo)   
            ax2b.CLim = ax2a.CLim;
            ax2b.XTick = [];
            ax2b.YTick = [];
    
            % additional plots
            ps = linspace(ax2b.XLim(1),ax2b.XLim(2),7);
            tops = ps(2:2:end);
            plot([tops; tops],ax2b.YLim,'w:')
    
            % text
            if uu==1
                text(-0.02,0.5,sprintf('Nose up'),'Units','normalized','HorizontalAlignment','center','FontSize',8,'Color','k','VerticalAlignment','bottom','rotation',90)
            end

        % plot pitch tuning
        ax3b = axes('Units','pixels','Position',[xvec(uu) ynow-70 xsiz ysiz]);
            % plot data
            imagesc([-180 180],[-90 90],mpit,'alphadata',~isnan(mpit))
    
            % axis settings
            ax3b.XTick = -180:90:180;
            ax3b.YTick = -90:45:90;            
            if uu==1
                xlabel(sprintf('Azimuth (%c)',176))
                ylabel(sprintf('Tilt (%c)',176))                
            else
                xlabel(sprintf('Azimuth (%c)',176))
                ax3b.YTickLabel = {};                
            end
            axis xy tight
            view(0,90);
            ax3b.FontSize = 8;
            ax3b.CLim = [0,max(mpit(:))];
            colormap(gca,turbo)
            % ax3b.XTickLabelRotation = 0;
    end

%% >>>>>>>>>> Azimuth correlation results
    xnow = xnow+20;
    ynow = ynow-280;
    fname = [config.data_out_dir 'PIT_correlations_pitch_shuffles.mat'];
    load(fname,'shuffs_pitch');

    % create axis
    ax2 = axes('Units','pixels','Position',[xnow ynow 220 130]);  
        ah = add_panel_title('C',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        pp = 2;
        var = 'within_session_stability';
        v2a = clumaa.(var)(pidx & clumaa.partn==pp,1); % within session stability, arena 1
        var = 'pitch_stability';
        v2b = clumaa.(var)(pidx & clumaa.partn==pp,1); % pitch stability, arena 1
        v2s = shuffs_pitch(:,2); % pitch stability shuffle, arena 1
        r = clumaa.repetition_score(pidx & clumaa.partn==1,1); % repetition score, arena 1
        cutoff = prctile(r,99);
        v3 = v2a(clumaa.repetition_score(pidx & clumaa.partn==pp,1)>cutoff);

        xi = -1:0.01:1;
        bw = 0.05;
        f1 = ksdensity(v2a(:),xi(:),"Bandwidth",bw);
        f2 = ksdensity(v2b(:),xi(:),"Bandwidth",bw);
        f3 = ksdensity(v2s(:),xi(:),"Bandwidth",bw);
        f4 = ksdensity(v3(:),xi(:),"Bandwidth",bw);

        % main plot    
        alph = 0.7;          
        a1 = area(xi,f1,'FaceColor',plot_set{1,2},'EdgeColor','k','FaceAlpha',alph); hold on;
        a2 = area(xi,f2,'FaceColor','m','EdgeColor','k','FaceAlpha',alph);       
        a3 = area(xi,f3,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--'); hold on;       
        a4 = area(xi,f4,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','-'); hold on;       

        % axis settings
        ax2.XLim = [-0.5 1]; 
        ax2.FontSize = 8;
        ytickformat('%.1f')
        box off
        xlabel(sprintf('Correlation ({\\itr})'))    
        ylabel(sprintf('PDF'))  

        [~,leg] = legendflex([a1 a2 a3 a4],...
            {sprintf('Session half 1 vs 2'),...
            sprintf('Nose up vs down'),...
            sprintf('Nose up vs down shuffle'),sprintf('Leftward vs rightward (terrain cells only)')},...
            'anchor',{'sw','sw'},'ncol',1,'box','off','buffer',[-10,-100],'xscale',.5,'fontsize',9); 
        
        leg(5).Children.FaceAlpha = a1.FaceAlpha;
        leg(6).Children.FaceAlpha = a2.FaceAlpha;
        leg(7).Children.FaceAlpha = a3.FaceAlpha;
        leg(8).Children.FaceAlpha = a4.FaceAlpha;            

        leg(5).Children.FaceColor = a1.FaceColor;
        leg(6).Children.FaceColor = a2.FaceColor;
        leg(7).Children.FaceColor = a3.FaceColor;
        leg(8).Children.FaceColor = a4.FaceColor;

%% >>>>>>>>>> Between session correlation results (Cohen's d values)
    % create axis
    ax = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+50 ynow 50 130]);
        % ah = add_panel_title('b',sprintf(''),'yoffset',0,'xoffset',-20,'width',400);
        e1 = computeEffectSize(v2a,v2b); % within vs pitch
        g1 = e1.cliff_delta;
        e3 = computeEffectSize(v2a,v3); % within vs pitch (terrain only)
        g3 = e3.cliff_delta;          
        e2 = computeEffectSize(v2b,v2s); % pitch vs shuffle
        g2 = e2.cliff_delta;

        stem(1:3,[g1 g3 g2],'filled','ko'); hold on;

        % axis settings
        ax.YLim = [-0.1 1];
        ax.XLim = [0.5 3.5];
        ylabel(sprintf('Cliff''s Delta'))   
        ax.XTick = 1:3;
        ax.XTickLabel = {};
        ytickformat('%.1f');
        ax.FontSize = 8;
        box off
        ax.XColor = 'none';

        % additional plots
        line(ax.XLim,[0 0],'Color',[.5 .5 .5]);

        % stats
        vcoef = 1;
        fsiz = 12;
        dat = {v2a v2b;v2a v3;v2b v2s};
        for ii = 1:3
            [z,p,~,~,~] = permutationZTest(dat{ii,1},dat{ii,2},'method','cliff');
            if p<.001
                text(ii,ax.YLim(2).*vcoef,'***','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<0.01
                text(ii,ax.YLim(2).*vcoef,'**','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            elseif p<.05
                text(ii,ax.YLim(2).*vcoef,'*','FontSize',fsiz,'HorizontalAlignment','center','VerticalAlignment','bottom');
            end
        end

        axi = axes('Units','pixels','Position',[ax.Position(1) ax.Position(2)-ax.Position(3) ax.Position(3) ax.Position(3)]);
            axi.Color = 'none';
            axi.XLim = ax.XLim;
            axi.YLim = ax.XLim;
            rsiz = 0.6;
            rectangle('Position',[1-rsiz/2 1.75 rsiz rsiz],'FaceColor',[a1.FaceColor a1.FaceAlpha],'EdgeColor',a1.EdgeColor,'LineStyle',a1.LineStyle);
            rectangle('Position',[1-rsiz/2 0.6 rsiz rsiz],'FaceColor',[a2.FaceColor a2.FaceAlpha],'EdgeColor',a2.EdgeColor,'LineStyle',a2.LineStyle);
            ty = 1.1;
            text(1,ty,'vs','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize)

            rectangle('Position',[2-rsiz/2 1.75 rsiz rsiz],'FaceColor',[a1.FaceColor a1.FaceAlpha],'EdgeColor',a1.EdgeColor,'LineStyle',a1.LineStyle);
            rectangle('Position',[2-rsiz/2 0.6 rsiz rsiz],'FaceColor',[a4.FaceColor a4.FaceAlpha],'EdgeColor',a4.EdgeColor,'LineStyle',a4.LineStyle);
            text(2,ty,'vs','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize)

            rectangle('Position',[3-rsiz/2 1.75 rsiz rsiz],'FaceColor',[a2.FaceColor a2.FaceAlpha],'EdgeColor',a2.EdgeColor,'LineStyle',a2.LineStyle);
            rectangle('Position',[3-rsiz/2 0.6 rsiz rsiz],'FaceColor',[a3.FaceColor a3.FaceAlpha],'EdgeColor',a3.EdgeColor,'LineStyle',a3.LineStyle);
            text(3,ty,'vs','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',ax.FontSize)

            axis off

%% >>>>>>>>>> Relationship azimuth stability repetition score
    xnow = xnow+380;
    ynow = ynow;

    % create axis
    ax2 = axes('Units','pixels','Position',[xnow ynow 160 120]);  
        ah = add_panel_title('D',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400,'fontsize',fs);                
    
        var = 'pitch_stability';
        v1 = clumaa.(var)(pidx & clumaa.partn==2,1); % pitch stability, hills
        var = 'repetition_score';
        v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % pitch stability, hills

        % main plot
        v = [v2(:) v1(:)];
        d = NaN(length(v(:,1)),1);
        v_norm = v ./ max(v,[],1,'omitmissing');
        sigma = 0.04;
        for j = 1:size(v_norm,1) % loop over every point
            D = sqrt(sum(bsxfun(@minus,v_norm(j,:),v_norm).^2,2));
            d(j) = sum( normpdf(D,0,sigma),'all','omitmissing');          % generates filtered distance measure
        end
        scatter(v2,v1,30,d,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5); hold on;

        % axis settings
        xlabel(sprintf('Repetition score'))            
        ylabel(sprintf('Tilt stability (r)'))    
        % ax.XTick = [0.001 0.01 0.1 1 10 100];
        % ax.YTick = [0.001 0.01 0.1 1 10 100];        
        ax2.XLim = [-0.5 1.5]; 
        ax2.YLim = [0 1];         
        ax2.FontSize = 8;
        % ax.YScale = 'log';
        % ax.XScale = 'log';
        ytickformat('%.1f')
        xtickformat('%.1f')        
        rf = refline;
        set(rf,'Color','k');
        grid on

        % % additional plots
        % line([ax.XLim(1) 100],[ax.XLim(1) 100],'Color',[.5 .5 .5],'LineStyle','--')
        % 
        % stats
        [r,p] = corr(v2,v1,'rows','pairwise','type','Spearman');
        text(1,0,sprintf('r = %.2f',r),'FontSize',8,'HorizontalAlignment','right','Units','normalized','VerticalAlignment','bottom','Color','k')
        % [r,p] = corr(fs(:,1),fs(:,3),'rows','pairwise','type','Pearson');    
        % text(0,1.02,sprintf('Arena 1 vs Arena 2: r = %.1f',r),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','#012030')

% % stats comparing directionality in terrain cells vs non terrain cells
% var = 'repetition_score';
% r = clumaa.(var)(pidx & clumaa.partn==1,1); % repetition score, hills
% cutoff = prctile(r,99);
% [ds,gs] = vectorDATAGROUP([],v1(v2<cutoff),v1(v2>cutoff));
% [p,a,s] = anova1(ds,gs,'off')
%                 keyboard

        % horizontal colormap
        axc = axes('Units','pixels','Position',[xnow ynow-60 100 10]);
            x = linspace(ax2.CLim(1),ax2.CLim(2),100);
            imagesc(x,'XData',ax2.CLim,'YData',[0 1]);
            colormap(ax2.Colormap);
            axis on
            axc.XTick = [];
            axc.YTick = [];
            axc.FontSize = 7;
            text(0,0.5,sprintf('%.1f ',0),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
            text(1,0.5,sprintf(' Max'),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',axc.FontSize,'Units','normalized')
            text(0.5,1,sprintf('Density'),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',axc.FontSize,'Units','normalized') 


% %% >>>>>>>>>> 3DHD correlation results
%     xnow = 40;
%     ynow = 100;
%     fname = [config.data_out_dir 'PIT_correlations_3dhd_shuffles.mat'];
%     load(fname,'shuffs_3dhd');
% 
%     % create axis
%     ax2 = axes('Units','pixels','Position',[xnow ynow 220 110]);  
%         ah = add_panel_title('c',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400);   
% 
%         var = 'hd3d_stability';
%         v1 = clumaa.(var)(pidx & clumaa.partn==1,2); % arena 1 vs arena 2
%         v2 = clumaa.(var)(pidx & clumaa.partn==1,1); % arena 1 vs hills
%         v4 = shuffs_3dhd(:,2); % arena 1 vs arena 2 shuffles
%         v5 = shuffs_3dhd(:,1); % arena 1 vs hills shuffles
% 
%         xi = -1:0.01:1;
%         bw = 0.05;
%         f1 = ksdensity(v1(:),xi(:),"Bandwidth",bw);
%         f2 = ksdensity(v2(:),xi(:),"Bandwidth",bw);
%         f4 = ksdensity(v4(:),xi(:),"Bandwidth",bw);
%         f5 = ksdensity(v5(:),xi(:),"Bandwidth",bw);
% 
%         % main plot
%         alph = 0.7;  
%         a1 = area(xi,f1,'FaceColor','#012030','EdgeColor','k','FaceAlpha',alph); hold on;
%         a2 = area(xi,f2,'FaceColor','#45C4B0','EdgeColor','k','FaceAlpha',alph);
%         a4 = area(xi,f4,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle','--'); hold on;       
%         a5 = area(xi,f5,'FaceColor',rgb('white'),'EdgeColor','k','FaceAlpha',0.10,'LineStyle',':'); hold on;       
% 
%         % axis settings
%         ax.XLim = [-0.5 1]; 
%         ax.FontSize = 10;
%         set(ax,'yticklabel',num2str(get(ax,'ytick')','%.1f'))
%         box off
%         xlabel(sprintf('Correlation ({\\itr})'))    
%         ylabel(sprintf('PDF'))    
% 
%         % additional plots
%         line(ax.XLim,[0 0],'Color',[.5 .5 .5]) 
% 
%         [~,leg] = legendflex([a1 a2 a4 a5],...
%             {sprintf('%s vs %s',maze_names{1},maze_names{3}),...
%             sprintf('%s vs %s',maze_names{1},maze_names{2}),...
%             sprintf('%s vs %s shuffle',maze_names{1},maze_names{3}),...
%             sprintf('%s vs %s shuffle',maze_names{1},maze_names{2})},...
%             'anchor',{'e','e'},'ncol',2,'box','off','buffer',[40,87],'xscale',.5,'fontsize',9); 
% 
%         leg(5).Children.FaceAlpha = a1.FaceAlpha;
%         leg(6).Children.FaceAlpha = a2.FaceAlpha;
%         leg(7).Children.FaceAlpha = a4.FaceAlpha;
%         leg(8).Children.FaceAlpha = a5.FaceAlpha;
% 
% %% >>>>>>>>>> Relationship azimuth stability repetition score
%     xnow = 40;
%     ynow = ynow;
% 
%     % create axis
%     ax2 = axes('Units','pixels','Position',[xnow ynow 160 110]);  
%         ah = add_panel_title('d',sprintf(''),'yoffset',-10,'xoffset',-10,'width',400);                
% 
%         var = 'hd3d_stability';
%         v1 = clumaa.(var)(pidx & clumaa.partn==1,1); % pitch stability, arena 1 vs hills
%         var = 'repetition_score';
%         v2 = clumaa.(var)(pidx & clumaa.partn==2,1); % repetition score, hills
% 
%         % main plot
%         scatter(v2,v1,30,'k','filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerFaceAlpha',0.5); hold on;
% 
%         % axis settings
%         xlabel(sprintf('Repetition score'))            
%         ylabel(sprintf('Pitch azimuth stability'))    
%         % ax.XTick = [0.001 0.01 0.1 1 10 100];
%         % ax.YTick = [0.001 0.01 0.1 1 10 100];        
%         ax2.XLim = [-0.5 1.5]; 
%         ax2.YLim = [-1 1];         
%         ax2.FontSize = 8;
%         % ax.YScale = 'log';
%         % ax.XScale = 'log';
%         ytickformat('%.1f')
%         xtickformat('%.1f')        
%         refline
%         grid on
% 
%         % % additional plots
%         % line([ax.XLim(1) 100],[ax.XLim(1) 100],'Color',[.5 .5 .5],'LineStyle','--')
%         % 
%         % stats
%         [r,p] = corr(v2,v1,'rows','pairwise','type','Spearman');
%         text(1,0,sprintf('r = %.2f',r),'FontSize',8,'HorizontalAlignment','right','Units','normalized','VerticalAlignment','bottom','Color','k')
%         % [r,p] = corr(fs(:,1),fs(:,3),'rows','pairwise','type','Pearson');    
%         % text(0,1.02,sprintf('Arena 1 vs Arena 2: r = %.1f',r),'FontSize',8,'HorizontalAlignment','left','Units','normalized','VerticalAlignment','bottom','Color','#012030')
% 
% 
%         return
keyboard
%% >>>>>>>>>> Save the overall figure
    if 1
        fname = [config.fig_dir '\Fig S9.png']; 
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








