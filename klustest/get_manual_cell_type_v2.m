function get_manual_cell_type_v2(data_dir,rnow,dnow)
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
    sname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    disp(sprintf('Loading cluma: %s',sname))    
    load(sname,'cluma'); % load saved session data
    pdata = cluma.Properties.CustomProperties.pdata;
    part_names = {'arena1','hills','arena2'};
    nparts = numel(unique(pdata.pos.session));
    disp(sprintf('\t...%d sessions',nparts));   
    disp(sprintf('\t...recording date: %s',dnow));            
    disp(sprintf('\t...done')) 
    tets = unique(cluma.ele);
    wav = cell(max(tets),1);
    clu = cell(max(tets),1);
    spt = cell(max(tets),1); 
    for tt = 1:length(tets)
        fname = [data_dir '\' rnow '\' dnow '\cluma.t' num2str(tets(tt)) '.spikes'];
        if exist(fname,'file')
            dat = CLUMA_load_merged_file(fname,'.spikes');
            wav{tets(tt)} = dat.wav;
        end
        fname = [data_dir '\' rnow '\' dnow '\cluma.t' num2str(tets(tt)) '.clusters'];
        if exist(fname,'file')
            dat = CLUMA_load_merged_file(fname,'.clusters');
            clu{tets(tt)} = dat.cut;
        end
        fname = [data_dir '\' rnow '\' dnow '\cluma.t' num2str(tets(tt)) '.spt'];
        if exist(fname,'file')
            dat = CLUMA_load_merged_file(fname,'.spt');
            spt{tets(tt)} = dat.spt;
        end
    end

%% >>>>>>>>>> Add width of waveform
    if ~any(ismember(cluma.Properties.VariableNames,'wave_width')) % if the column(s) do not exist yet
        cluma.wave_width = NaN(size(cluma,1),1); % preallocate
        cluma.wave_amps = NaN(size(cluma,1),1); % preallocate
    
        ucis = unique(cluma.uci); % list of unique cells in sdata  
        for uu = 1:length(ucis)
            for pp = 1:nparts
                idx = find(ismember(cluma.uci,ucis{uu}) & cluma.partn==pp);
                if isempty(idx)
                    continue
                end
                if cluma.n_spikes(idx)==0
                    continue
                end
                cluster_now = cluma.clu(idx);            
                tt = cluma.ele(idx);
                waves = wav{tt}; % waveforms for this tetrode
                clus = clu{tt}; % clusters for this tetrode
                spt_now = spt{tt}; % spike times for this tetrode
                part_times = pdata.session_times(pp,:);
    
                y_mat = waves(:,:,clus==cluster_now & spt_now>part_times(1) & spt_now<part_times(2)); % cut the waveforms    
                meanz = cell(2,4);                
                for ww = 1:size(y_mat,2) % for every channel
                    wnow = squeeze(y_mat(:,ww,:));
                    if isempty(wnow)
                        continue
                    end
                    wav_mean = mean(wnow,2,'omitnan');
                    meanz{1,ww} = wav_mean;
                    wav_std = std(wnow,[],2,'omitnan');
                    meanz{2,ww} = wav_std;                    
                end
                mxs = cell2mat( cellfun(@max,meanz(1,:),'UniformOutput',false) );
                [~,widx] = sort(mxs,'descend'); % sort from largest > smallest waveform
                max_wav_means = meanz(1,widx);
    
                % use the waveform with the largest amplitude to calculate wave characteristics
                wnow = max_wav_means{1};
                [~,maxvalt] = max(wnow);
                wnow(1:maxvalt) = NaN;
                [~,minvalt] = min(wnow,[],'omitnan');
                wi = pdata.waveform_xi;
                width_ms = wi(minvalt) - wi(maxvalt); % waveform width in ms
    
                % accumulate data
                cluma.wave_width(idx) = single(width_ms);
                cluma.wave_amps(idx) = single(mxs(1));
            end
        end
    end

%% ##################################### Run through every cluster
    ucis = unique(cluma.uci); % list of unique cells in sdata  
    for uu = 1:length(ucis)
        uci = ucis{uu};
        disp(sprintf('\tCell %d of %d (%.f%%): %s',uu,length(ucis),uu/length(ucis)*100,uci))
        
        % clear figure for next cluster
        clf(fig_overall);     
        delete(ann);
        
        % add an annotation to the figure with some important info
        ann_str = sprintf('Cell: %s, Analysed: %s',uci,datestr(now,'yyyy-mm-dd-HH-MM-SS'));
        ann = annotation('textbox',[0, 1, 1, 0],'string',ann_str,'FontSize',12,'LineStyle','none','interpreter','none');  
        
        idx = find(ismember(cluma.uci,uci));
        
        %% add waveform characteristics
        mwaves = cluma.mean_waveforms(idx,:);
        pos = pdata.pos;

%% #################### automated cell typing
        [~,midx] = max(cluma.f_rate_hz(idx));
        idx = idx(midx);
        
        dat = NaN(1,6); % to hold data
        frate = cluma.f_rate_hz(idx); % firing rate
        wow = cluma.wave_width(idx); % width of waveform, this is for the channel with the highest mean amplitue
        
        spatz1 = cluma.planar_spatial_info_shuffles(idx,2); % z-scored spatial information (planar)
        spatz2 = cluma.surficial_spatial_info_shuffles(idx,2); % z-scored spatial information (surficial)
        spat1 = cluma.planar_spatial_info_shuffles(idx,1); % actual spatial information (planar)
        spat2 = cluma.surficial_spatial_info_shuffles(idx,1); % actual spatial information (surficial)
        
        gscorez1 = cluma.planar_spatial_info_shuffles(idx,4); % z-scored grid score (planar)
        gscorez2 = cluma.surficial_spatial_info_shuffles(idx,4); % z-scored grid score (surficial)
        gscore1 = cluma.planar_spatial_info_shuffles(idx,3); % actual grid score (planar)
        gscore2 = cluma.surficial_spatial_info_shuffles(idx,3); % actual grid score (surficial)        
        
        rayz1 = cluma.planar_spatial_info_shuffles(idx,6); % z-scored rayleigh (planar)
        rayz2 = cluma.surficial_spatial_info_shuffles(idx,6); % z-scored rayleigh (surficial)
        ray1 = cluma.planar_spatial_info_shuffles(idx,5); % actual rayleigh (planar)
        ray2 = cluma.surficial_spatial_info_shuffles(idx,5); % actual rayleigh (surficial)
                
        part_times = pdata.session_times(midx,:);
        tt = cluma.ele(idx);        
        spt_now = spt{tt}; % spike times for this tetrode    
        clu_now = clu{tt}; % clusters for this tetrode        
        sindax = logical(sum(spt_now' > part_times(1) & spt_now' < part_times(2) & (clu_now==cluma.clu(idx))',1)); % spikes          
        spt2 = spt_now(sindax);               
        rp = 0.002; % refractory period in ms
        rpv = sum(diff(spt2(:))<rp);
        rpvprop = rpv./numel(spt2).*100;

        rate_min = 0.1;
        rate_max = 10;
        width_cutoff = 0.250;
        spatial_cutoff = 1; % b/s
        grid_cutoff = 0.8; % g-score
        hd_cutoff = 0.2; % rayleigh
        z_cutoff = 5; % z-score cutoff for shuffles
        rpv_cutoff = 5; % percentage cutoff for refractory period violations
        % 1 = noise
        % 2 = pyramidal
        % 3 = place cell
        % 4 = grid cell
        % 5 = hd cell
        % 6 = interneuron
        cell_type = 1; % noise by default
        cols = {'g','w','w','w','w','w','w'}; % for plotting later        
        if frate > rate_min && rpvprop < rpv_cutoff % if the cluster's firing is not too low and there are not too many RPVs
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
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{7},'String','Spatial cell','Position',[bxpos ymx-240 100 50],'Callback',{@change_ctype,7});           
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{6},'String','Interneuron','Position',[bxpos ymx-320 100 50],'Callback',{@change_ctype,6});              
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor','m','String','Agree','Position',[bxpos ymx-400 100 50],'Callback',{@change_ctype,cell_type});              
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{4},'String','Grid cell','Position',[bxpos ymx-480 100 50]);     
        uicontrol(fig_overall,'Style','pushbutton','Backgroundcolor',cols{5},'String','HD cell','Position',[bxpos ymx-560 100 50]);  
        
        xmin = 20;
        xbuff = 250;
        xvec = [xmin xmin+xbuff xmin+2*xbuff xmin+3.16*xbuff xmin+3.9*xbuff+20 xmin+4.8*xbuff];
        ymax = 450;
        ybuff = 200;
        yvec = [ymax ymax-ybuff ymax-2*ybuff ymax-3*ybuff ymax-4*ybuff ymax-5*ybuff];
        pwidth = 220;
        pheight = 190;        
        
        for pp = 1:nparts
            idx = find(ismember(cluma.uci,uci) & cluma.partn==pp);
            if isempty(idx)
                continue
            end
            if isempty(cluma.spike_index{idx})
                continue
            end

            part_now = part_names{pp};
            part_duration = cluma.session_duration_s(idx);

            session_times = pdata.session_times;
            pidx = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part

            ppox = pos.pox_planar(pidx);
            ppoy = pos.poy_planar(pidx);
            ppot = pos.pot(pidx);
            pspx = pos.pox_planar(cluma.spike_index{idx});
            pspy = pos.poy_planar(cluma.spike_index{idx});
            ratemap = cluma.ratemap_planar{idx};
            amap = cluma.planar_amap{idx};           

%% >>>>>>>>>> Positions and spikes
            axps = axes('Units','pixels','Position',[xvec(1),yvec(pp),pwidth,pheight]);
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
            axrt = axes('Units','pixels','Position',[xvec(2),yvec(pp),pwidth,pheight]);
                im = imagesc(ratemap,'alphadata',~isnan(ratemap));
                daspect([1 1 1])
                caxis([0 max([0.1 max(ratemap(:),[],'omitnan')])])       
                colormap(axrt,turbo);
                axis xy off

                sinfo = cluma.planar_spatial_info_shuffles(idx,:);            
                text(0,1.05,sprintf('SI: %.2f, z: %.2f',sinfo(1),sinfo(2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')

                axp = get(gca,'Position');
                ccb = colorbar;
                set(gca,'Position',axp);
                set(ccb,'Position',get(ccb,'Position')+[0 0 -0.004 0])
                title(ccb,'Hz','FontSize',fsiz)        
                set(ccb,'yticklabel',num2str(get(ccb,'ytick')','%.1f')) % change Ytick labels to have the same number of decimal places
                
%% >>>>>>>>>> Grid autocorrelogram              
            axgs = axes('Units','pixels','Position',[xvec(3),yvec(pp),pwidth,pheight]);
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
            tt = cluma.ele(idx);
            wav_now = wav{tt}; % waveforms for this tetrode
            clu_now = clu{tt}; % clusters for this tetrode
            spt_now = spt{tt}; % spike times for this tetrode
            part_times = pdata.session_times(pp,:);

            sindax = logical(sum(spt_now' > part_times(1) & spt_now' < part_times(2) & (clu_now==cluma.clu(idx))',1)); % spikes                    
            
            axwav = axes('Units','pixels','Position',[xvec(4),yvec(pp)+30,180,180]);
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

                wavtime = pdata.waveform_xi; 
                if isempty(wnow) % if there are no spikes
                    wnow = NaN(size(wavtime'));
                end

                plot(wavtime,wnow','Color','k'); hold on;
                [hl,hp] = boundedline(wavtime,wav_means(widx(1),:),wav_stds(widx(1),:),'-k');
                set(hl,'Color','r','LineStyle','-','LineWidth',1); 
                set(hp,'FaceColor','b','FaceAlpha',.5);     
                ax = gca;
                ax.XLim = [-0.25 0.75];
                ax.YDir = 'normal';
                box on
                grid on
                if pp<nparts
                    ax.XTick = [];  
                else
                    xlabel('Time (ms)')
                end                
                ylabel(sprintf('Amplitude (%cV)',char(956)))
                text(1,1.05,sprintf('Ch%d, Peak: %.1f%cV, Width: %.2f%cs',widx(1),mxs(widx(1)),char(956),cluma.wave_width(idx).*1000,char(956)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','right')
                if max_plot_waves>0 % if we want to limit the number of waveforms shown
                    text(0.99,1,sprintf('N = %d',max_plot_waves_now),'Units','normalized','FontSize',8,'HorizontalAlignment','right','VerticalAl','top')           
                end               

%% >>>>>>>>>> Spike autocorrelation - 25ms refractory period             
            axref = axes('Units','pixels','Position',[xvec(5),yvec(pp)+30,180,180]);               
                [autoc,tlag,e] = CLUMA_spike_auto(spt_now(sindax),'bin_size',1,'win_size',30,'method','correlogram');

                bar(tlag,autoc,1,'k','EdgeColor','k'); hold on;
                bar(tlag(abs(tlag)<2),autoc(abs(tlag)<2),1,'r','EdgeColor','r'); hold on;

                ax = gca;            
                ax.XLim = [-20 20];
                if pp<nparts
                    ax.XTick = [];    
                end
                ax.YTick = [];
                text(0,1.05,sprintf('RPV: %d (%.2f%%)',rpv,rpvprop),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','left')
                set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))        
                grid on;

%% >>>>>>>>>> Spike autocorrelation - 500ms theta modulation           
            axref = axes('Units','pixels','Position',[xvec(6),yvec(pp)+30,180,180]);
                [autoc,tlag,e] = CLUMA_spike_auto(spt_now(sindax),'bin_size',10,'win_size',500,'method','correlogram');

                bar(tlag,autoc,1,'k'); hold on;

                ax = gca;
                ax.XLim = [-500 500];
                if pp<nparts
                    ax.XTick = [];    
                end
                ax.YTick = [];
                % text(.5,1.05,sprintf('Theta index: %.1f, frequency: %.2fHz',sdata.autocorr_500_info(pidx,1),sdata.autocorr_500_info(pidx,2)),'Units','normalized','FontSize',fsiz,'HorizontalAlignment','center')
                set(ax,'yticklabel',num2str(get(gca,'ytick')','%.3f'))        
                grid on;        
        
        end

%% #################### get users response and add it to sdata
        uiwait(fig_overall); % at this point start waiting for a cell type button to be pressed

        if ~any(ismember(cluma.Properties.VariableNames,'cell_type')) % if the column(s) do not exist yet
            cluma.cell_type = NaN(size(cluma,1),2); % preallocate
        end
        cval2 = cell_type; % curated value

        cluma.cell_type(ismember(cluma.uci,uci),:) = repmat([cval1 cval2],sum(ismember(cluma.uci,uci)),1);  % automated cell type, manual cell type 

    end
    close(fig_overall)

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    cname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    save(cname,'cluma'); % save session data
    analysis_log({'cell_type_v2'},1,'version',{'v2.0.0'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('cell typing has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dumb pushbutton control function
function change_ctype(~,~,value) 
    global cell_type; 
    cell_type = value; 
    uiresume(gcbf);
end









































