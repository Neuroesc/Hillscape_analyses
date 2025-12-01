function get_manual_cell_type
%NAME short desc. 
% long desc                                  
%
%   out = name(in) process with default settings
%
%   out = name(in,optional) process using optional argument 1
%
%   out = rate_mapper(__,name,value) process with Name-Value pairs 
%
%   parameters include:
%
%   'in'            -   Scalar, positive integer that specifies X, units are in Y.
%
%                       Default value is Z.
%
%   outputs include:
%
%   'out'           -   [Nx2] desc
%
%
%   Class Support
%   -------------
%   The input matrix in must be a real, non-sparse matrix of
%   the following classes: uint8, int8, uint16, int16, uint32, int32,
%   single or double.
%
%   Notes
%   -----
%   1. 
%
%   2. 
%
%   Example
%   ---------
% 
%   See also NAME, NAME

% HISTORY:
% version 1.0.0, Release 17/01/18 Initial release
% version 2.0.0, Release 29/06/19 v2 created from determineCELLTYPE
% version 3.0.0, Release 20/02/23 renamed to get_manual_cell_type for 3D hills data
% version 4.0.0, Release 20/02/23 overhauled for 3D hills data
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
    figvis = 'on';
    global cell_type;         
    outname = 'klustest';
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    stk = dbstack;
    tnow = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',stk.name,tnow))
    
    fig_overall = figure('visible',figvis,'Units','pixels','Position',[20, 150, 1500, 700]); % open/create figure
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    fsiz = 8; % the fontsize for different texts
    ann = annotation('textbox',[0, 1, 1, 0],'string','','FontSize',8,'LineStyle','none','interpreter','none');  

%% >>>>>>>>>> Prepare the data
    sname = [pwd '\' outname '\sdata.mat'];
    disp(sprintf('\t...loading sdata: %s',sname))    
    load(sname,'sdata'); % load saved session data
    pdata = sdata.Properties.CustomProperties.pdata;
    part_config = pdata.part_config;    
    [~,wav,~] = get_spk_for_klustest('Neuralynx',pdata.data_dirs,pdata.tetrodes,0);
    
    pdata = sdata.Properties.CustomProperties.pdata;
    part_config = pdata.part_config;
    nparts = size(part_config,1); 
    
%% ##################################### Run through every cluster
    ucis = unique(sdata.uci); % list of unique cells in sdata  
    
    for uu = 1:length(ucis)
        uci = ucis{uu};
        disp(sprintf('\tCell %d of %d (%.f%%): %s',uu,length(ucis),uu/length(ucis)*100,uci))
        
        % clear figure for next cluster
        clf(fig_overall);     
        delete(ann);
        
        % add an annotation to the figure with some important info
        ann_str = sprintf('Cell: %s, Analysed: %s',uci,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        ann = annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
        
        idx = find(ismember(sdata.uci,uci));
        
%% #################### automated cell typing
        [~,midx] = max(sdata.frate(idx));
        idx = idx(midx);
        
        dat = NaN(1,6); % to hold data
        frate = sdata.frate(idx); % firing rate
        wow = sdata.wave_width(idx); % width of waveform, this is for the channel with the highest mean amplitue
        
        spatz1 = sdata.planar_spatial_info_shuffles(idx,2); % z-scored spatial information (planar)
        spatz2 = sdata.surficial_spatial_info_shuffles(idx,2); % z-scored spatial information (surficial)
        spat1 = sdata.planar_spatial_info_shuffles(idx,1); % actual spatial information (planar)
        spat2 = sdata.surficial_spatial_info_shuffles(idx,1); % actual spatial information (surficial)
        
        gscorez1 = sdata.planar_spatial_info_shuffles(idx,4); % z-scored grid score (planar)
        gscorez2 = sdata.surficial_spatial_info_shuffles(idx,4); % z-scored grid score (surficial)
        gscore1 = sdata.planar_spatial_info_shuffles(idx,3); % actual grid score (planar)
        gscore2 = sdata.surficial_spatial_info_shuffles(idx,3); % actual grid score (surficial)        
        
        rayz1 = sdata.planar_spatial_info_shuffles(idx,6); % z-scored rayleigh (planar)
        rayz2 = sdata.surficial_spatial_info_shuffles(idx,6); % z-scored rayleigh (surficial)
        ray1 = sdata.planar_spatial_info_shuffles(idx,5); % actual rayleigh (planar)
        ray2 = sdata.surficial_spatial_info_shuffles(idx,5); % actual rayleigh (surficial)
                
        rate_min = 0.1;
        rate_max = 10;
        width_cutoff = 0.250;
        spatial_cutoff = 1; % b/s
        grid_cutoff = 0.8; % g-score
        hd_cutoff = 0.2; % rayleigh
        z_cutoff = 5; % z-score cutoff for shuffles
        % 1 = noise
        % 2 = pyramidal
        % 3 = place cell
        % 4 = grid cell
        % 5 = hd cell
        % 6 = interneuron
        cell_type = 1; % noise by default
        cols = {'g','w','w','w','w','w'}; % for plotting later        
        if frate > rate_min % if the cluster's firing is not too low
            if frate < rate_max % if the cell's firing rate does not exceed the maximum
                if wow > width_cutoff % if the cell has a pyramidal waveform
                    cell_type = 2; % pyramidal
                    cols{1} = 'w';
                    cols{2} = 'g';
                    if (spatz1>z_cutoff && spat1>spatial_cutoff) || (spatz2>z_cutoff && spat2>spatial_cutoff) % if the pyramidal cell is spatial
                        cell_type = 3; % place cell
                        cols{1} = 'w';   
                        cols{2} = 'w';                        
                        cols{3} = 'g';
                    end                     
                else
                    cell_type = 6; % low firing interneuron
                    cols{1} = 'w';                    
                    cols{6} = 'g';
                end   
            else
                cell_type = 6; % high firing interneuron
                cols{1} = 'w';                    
                cols{6} = 'g';                
            end
            if (gscorez1>z_cutoff && gscore1>grid_cutoff) || (gscorez2>z_cutoff && gscore2>grid_cutoff) % if the grid score is significant
                cols{1} = 'w';                    
                cols{4} = 'b';                
            end 
            if (rayz1>z_cutoff && ray1>hd_cutoff) || (rayz2>z_cutoff && ray2>hd_cutoff) % if the grid score is significant
                cols{1} = 'w';                    
                cols{5} = 'b';                
            end
        end
        cval1 = cell_type; % automated value
        
%% #################### manual refining
        % pushbuttons for cell typing
        bxpos = 1400;
        ymx = 600;
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{1},'String','Noise','Position',[bxpos ymx 100 50],'Callback',{@change_ctype,1});  
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{2},'String','Pyramidal cell','Position',[bxpos ymx-80 100 50],'Callback',{@change_ctype,2}); 
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{3},'String','Place cell','Position',[bxpos ymx-160 100 50],'Callback',{@change_ctype,3});   
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{6},'String','Interneuron','Position',[bxpos ymx-240 100 50],'Callback',{@change_ctype,6});              
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor','m','String','Agree','Position',[bxpos ymx-320 100 50],'Callback',{@change_ctype,cell_type});              
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{4},'String','Grid cell','Position',[bxpos ymx-400 100 50]);     
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{5},'String','HD cell','Position',[bxpos ymx-480 100 50]);  
        
        xmin = 20;
        xbuff = 250;
        xvec = [xmin xmin+xbuff xmin+2*xbuff xmin+3.16*xbuff xmin+3.9*xbuff+20 xmin+4.8*xbuff];
        ymax = 450;
        ybuff = 200;
        yvec = [ymax ymax-ybuff ymax-2*ybuff ymax-3*ybuff ymax-4*ybuff ymax-5*ybuff];
        pwidth = 220;
        pheight = 190;        
        
        idx = find(ismember(sdata.uci,uci));
        nsess = length(idx);
        for ss = 1:nsess
            pidx = idx(ss);
            part_now = part_config.part_names{sdata.partn(pidx)};
            part_duration = pdata.part_config.part_duration(sdata.partn(pidx));
            ppox = pdata.(part_now).pox;
            ppoy = pdata.(part_now).poy;
            ratemap = sdata.ratemap_planar{pidx};
            amap = sdata.planar_amap{pidx};         
            ppot = pdata.(part_now).pot;         
            pspx = ppox(sdata.spt_pot_index{pidx});
            pspy = ppoy(sdata.spt_pot_index{pidx});   

%% >>>>>>>>>> Positions and spikes
            axps = axes('Units','pixels','Position',[xvec(1),yvec(ss),pwidth,pheight]);
                % plot position data, excluding pieces not included in this part
                % by inserting NaNs between intervals we can plot this as one line, which saves on memory
                dindax = abs([0; diff(ppot)])>0.1;
                pos_plot = [ppox ppoy];
                pos_plot(dindax,:) = NaN;
                plot(pos_plot(:,1),pos_plot(:,2),'Color',[.5 .5 .5 .5]); hold on;

                % plot spikes after position data so they are all on top     
                plot(pspx,pspy,'Color',[1 0 0 0.5],'Marker','.','LineStyle','none') 

                % additional settings
                daspect([1 1 1])
                axis xy off tight
                text(0,1.1,sprintf('%d spikes (%.2f Hz), %d seconds (%.1f mins)',numel(pspx),numel(pspx)/part_duration,round(part_duration),part_duration/60),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
        
%% >>>>>>>>>> Firing rate map       
            axrt = axes('Units','pixels','Position',[xvec(2),yvec(ss),pwidth,pheight]);
                im = imagesc(ratemap,'alphadata',~isnan(ratemap));
                daspect([1 1 1])
                caxis([0 max([0.1 max(ratemap(:),[],'omitnan')])])       
                colormap(axrt,flipud(viridis));
                axis xy off

                sinfo = sdata.planar_spatial_info_shuffles(pidx,:);            
                text(0,1.05,sprintf('SI: %.2f, z: %.2f',sinfo(1),sinfo(2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')

                axp = get(gca,'Position');
                ccb = colorbar;
                set(gca,'Position',axp);
                set(ccb,'Position',get(ccb,'Position')+[0 0 -0.004 0])
                title(ccb,'Hz','FontSize',fsiz)        
                set(ccb,'yticklabel',num2str(get(ccb,'ytick')','%.1f')) % change Ytick labels to have the same number of decimal places
                
%% >>>>>>>>>> Grid autocorrelogram              
            axgs = axes('Units','pixels','Position',[xvec(3),yvec(ss),pwidth,pheight]);
                imc = imagesc(amap,'alphadata',~isnan(amap));
                daspect([1 1 1])
                caxis([-0.2 1])
                colormap(axgs,turbo);        
                axis xy off

                text(0,1.05,sprintf('G: %.2f, z: %.2f',sinfo(3),sinfo(4)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')

                axp = get(gca,'Position');
                ccb = colorbar;
                set(gca,'Position',axp);
                set(ccb,'Position',get(ccb,'Position')+[0 0 -0.004 0])
                title(ccb,'Hz','FontSize',fsiz)        
                set(ccb,'yticklabel',num2str(get(ccb,'ytick')','%.1f')) % change Ytick labels to have the same number of decimal places

%% >>>>>>>>>> Waveforms      
            tt = sdata.tetrode(pidx);
            wav_now = wav{tt}; % waveforms for this tetrode

            clus = pdata.clusters;
            clu = clus{tt}; % clusters on this tetrode
            cnow = sdata.cluster(pidx); 
            spt = pdata.spike_times{tt}; % spike times for this tetrode
            part_times = part_config.part_times{sdata.partn(pidx)}; 
            sindax = logical(sum(spt' > part_times(:,1) & spt' < part_times(:,2) & (clu==cnow)',1)); % spikes                    
            
            axwav = axes('Units','pixels','Position',[xvec(4),yvec(ss)+30,180,180]);
                nch = 4;
                mxs = NaN(nch,1);
                wav_means = NaN(nch,size(wav_now,1));
                wav_stds = NaN(nch,size(wav_now,1));  
                wav_now_clus = cell(1,4);
                for ww = 1:nch % for every channel
                    wav_now_clus{ww} = squeeze(wav_now(:,ww,sindax))'; % should be [spikes x samples]
                    wav_means(ww,:) = mean(wav_now_clus{ww},1,'omitnan');
                    wav_stds(ww,:) = std(double(wav_now_clus{ww}),[],1,'omitnan') ./ sqrt(size(wav_now_clus{ww},1));
                    mxs(ww,1) = max(wav_means(ww,:),[],'omitnan');
                end
                [~,widx] = sort(mxs,'descend'); % sort from largest > smallest waveform

                % plot the waveform(s) with the highest amplitude
                wnow = squeeze(wav_now_clus{widx(1)});
                max_plot_waves = 250;
                if max_plot_waves>0 % if we want to limit the number of waveforms shown
                    max_plot_waves_now = min([size(wnow,1) max_plot_waves]); % if there are less than max_plot_waves waveforms, plot them all
                    rindx = randperm(size(wnow,1),max_plot_waves_now);
                    wnow = wnow(rindx,:);
                end

                wavtime = pdata.wavtime; 
                if isempty(wnow) % if there are no spikes
                    wnow = NaN(size(wavtime'));
                end

                plot(wavtime,wnow','Color','k'); hold on;
                [hl,hp] = boundedline(wavtime,wav_means(widx(1),:),wav_stds(widx(1),:),'-k');
                set(hl,'Color','r','LineStyle','-','LineWidth',1); 
                set(hp,'FaceColor','b','FaceAlpha',.5);     
                ax = gca;
                ax.XLim = [-0.25 0.75];
                ax.XTick = [];
                ax.YDir = 'normal';
                box on
                grid on
                xlabel('Time (ms)')
                ylabel(sprintf('Amplitude (%cV)',char(956)))
                text(1,1.05,sprintf('Ch%d, Peak: %.1f%cV, Width: %.2fms',widx(1),mxs(widx(1)),char(956),sdata.wave_width(pidx)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','right')
                if max_plot_waves>0 % if we want to limit the number of waveforms shown
                    text(0.99,1,sprintf('N = %d',max_plot_waves_now),'Units','normalized','FontSize',8,'HorizontalAlignment','right','VerticalAl','top')           
                end               

%% >>>>>>>>>> Spike autocorrelation - 25ms refractory period             
            axref = axes('Units','pixels','Position',[xvec(5),yvec(ss)+30,180,180]);
                autoc = sdata.autocorr_25{pidx};
                tlag = pdata.autocorr_25_xvalues;

                bar(tlag,autoc,1,'k','EdgeColor','k'); hold on;
                bar(tlag(abs(tlag)<2),autoc(abs(tlag)<2),1,'r','EdgeColor','r'); hold on;

                ax = gca;            
                ax.XLim = [-20 20];
                ax.XTick = [];      
                ax.YTick = [];
                text(0,1.05,sprintf('RPV: %d (%.2f%%)',sdata.autocorr_25_info(pidx,1),sdata.autocorr_25_info(pidx,2)*100),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
                set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))        
                grid on;

%% >>>>>>>>>> Spike autocorrelation - 500ms theta modulation           
            axref = axes('Units','pixels','Position',[xvec(6),yvec(ss)+30,180,180]);
                autoc = sdata.autocorr_500{pidx}(:,1);
                tlag = pdata.autocorr_500_xvalues;

                bar(tlag,autoc,1,'k'); hold on;

                ax = gca;
                ax.XLim = [-500 500];
                ax.XTick = [];      
                ax.YTick = [];
                text(.5,1.05,sprintf('Theta index: %.1f, frequency: %.2fHz',sdata.autocorr_500_info(pidx,1),sdata.autocorr_500_info(pidx,2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
                set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))        
                grid on;        
        
        end

%% #################### get users response and add it to sdata
        uiwait(fig_overall); % at this point start waiting for a cell type button to be pressed

        if ~any(ismember(sdata.Properties.VariableNames,'cell_type')) % if the column(s) do not exist yet
            sdata.cell_type = NaN(size(sdata,1),2); % preallocate
        end
        cval2 = cell_type; % curated value
        
        sdata.cell_type(idx,:) = repmat([cell_type cval1==cval2],size(idx,1),1);  % cell type, agreement between user and automation   

    end
    close(fig_overall)

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    sdata = rmprop(sdata,'pdata');
    sdata = addprop(sdata,{'pdata'},{'table'});
    sdata.Properties.CustomProperties.pdata = pdata;
    save([pwd '\' pdata.outname '\sdata.mat'],'sdata'); % save session data
    analysis_log({'cell_type'},1,'version',{'v4.0.0'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('cell typing has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',[pwd '\' pdata.outname],' &'');">','klustest folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dumb pushbutton control function
function change_ctype(~,~,value) 
    global cell_type; 
    cell_type = value; 
    uiresume(gcbf);
end









































