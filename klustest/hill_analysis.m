function hill_analysis
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short description
% long description
%
% USAGE:
%       [out] = template(in,in2)
%
% INPUT:
%       in - input 1
%       in2 - input 2
%
% OUTPUT:
%       p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 17/02/23 Initial release, adapted from GIT_hilltest
% version 1.0.1, Release 17/02/23 Most position data transforms already done
% version 1.0.2, Release 17/02/23 Most mapping already done
% version 2.0.0, Release 18/02/23 Switched to rate_mapper for firing rate maps
% version 3.0.0, Release 18/02/23 Switched to get_grid_score for grid scores
% version 4.0.0, Release 19/02/23 Created savelli_shuffle for determining spatial significance
% version 5.0.0, Release 19/02/23 Created fit_hill_frame for determining maze boundaries
% version 5.1.0, Release 19/02/23 Extensive bug testing of savelli_shuffle
% version 5.1.1, Release 19/02/23 Added transform of maze frame for surficial data
% version 6.0.0, Release 19/02/23 Dwell maps and position projections are added to sdata (first cell only)

%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
    %% Map settings
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
    
    % HD analysis settings
    mapset.hd_type      = 'histogram'; % (default 'histogram') enter 'density' for a circular kolmogorov smirnov density estimate plot, enter 'histogram' for the traditional histogram polar plot
    mapset.hd_bins      = 64; % (default 64) the number of bins to use when computing HD plot
    mapset.hd_boxcar    = 3; % (defualt 3) the number of bins over which to compute the HD histogram boxcar average        
    
    outname = 'klustest';
    override.maze_frame = 0; % 1 = do not load a maze frame from file
    skipfigs = 1; % 1 = skip making a figure if it exists already    
    save_figs = 1; % 1 = make and save figures (VERY time consuming, takes about 90% of klustest's time)
    fast_figures = 1; % 1 = 4-5x faster figure saving, but at a slightly lower quality
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    stk = dbstack;
    tnow = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',stk.name,tnow))
    
%% >>>>>>>>>> Prepare the data
    sname = [pwd '\' outname '\sdata.mat'];
    disp(sprintf('Loading sdata: %s',sname))    
    load(sname,'sdata'); % load saved session data
    pdata = sdata.Properties.CustomProperties.pdata;
    part_config = pdata.part_config;
    nparts = size(part_config,1);
    disp(sprintf('\t...%d sessions',size(pdata.sessions,1)));   
    disp(sprintf('\t...recording date: %s',pdata.date));            
    disp(sprintf('\t...done'))    

%% >>>>>>>>>> Analyse trajectory and add 3D HD
    disp(sprintf('Loading 3D trajectory...'))       
    for pp = 1:nparts % for every part     
        part_now = part_config.part_names{pp};
        if pp==1
            disp(sprintf('\t\t%s',part_now))                   
        else
            disp(sprintf('\b | %s',part_now))       
        end
        
        % cut the position data to include only this part
        part_times = part_config.part_times{pp};                
        pos = pdata.pos;
        pot = pos.pot;
        pindax = logical(sum(pot' >= part_times(:,1) & pot' <= part_times(:,2),1));
        
        pdata.(part_now).pot = pos.pot(pindax,1); % pos time for this part            
        pdata.(part_now).pox = pos.pox(pindax,1); % pos x for this part
        pdata.(part_now).poy = pos.poy(pindax,1); % pos y for this part
        pdata.(part_now).poz = pos.poz(pindax,1); % pos y for this part            
        pdata.(part_now).yaw = pos.poh(pindax,1); % yaw HD for this part
        pdata.(part_now).pit = pos.pitch(pindax,1); % pitch HD for this part
        pdata.(part_now).rol = pos.roll(pindax,1); % roll HD for this part
    end
    disp(sprintf('\t...done'))    
    
%% >>>>>>>>>> Add maze boundaries
    disp(sprintf('Calculating maze boundaries...'))    
    mname = [pwd '\' pdata.outname '\maze_boundaries.mat'];    
    if exist(mname,'file') && ~override.maze_frame && 1
        disp(sprintf('\t\t...loading from file'))                                    
        load(mname,'lx','ly','lz','-mat');             
    else
        [lx,ly,lz] = fit_hill_frame(pdata,0); % fit the maze frame to the data for all 3 parts simultaneously
        disp(sprintf('\t\t...saving to file'))                            
        save(mname,'lx','ly','lz','-mat');     
    end
    pdata.maze_frame = single([lx(:) ly(:) lz(:)]);      
    disp(sprintf('\t...done'))    
    
%% >>>>>>>>>> Add position data projections  
    disp(sprintf('Calculating trajectory projections...'))    

    for pp = 1:nparts % for every part  
        part_now = part_config.part_names{pp};
        if pp==1
            disp(sprintf('\t\t%s',part_now))                   
        else
            disp(sprintf('\b | %s',part_now))       
        end
        
        ppox = double( pdata.(part_now).pox ); % pos x for this part
        ppoy = double( pdata.(part_now).poy ); % pos y for this part
        ppoz = double( pdata.(part_now).poz ); % pos z for this part
        
        switch part_now
            case {'arena1','arena2'}
                % if this is an arena session
                pos_planar = [ppox,ppoy,ones(size(ppox))];                    
                pdata.(part_now).pos_planar = single(pos_planar); % surficial and planar are the same in this maze                
                
            case {'hills'}
                % if this is a hills session
                % k (frequency) = (2*pi)/period

%                 if ~isfield(pdata.(part_now),'pos_curve')
                    % fit curve to x,z data (side view to capture hill shape)
                    [x,y] = prepareCurveData( ppox, -ppoz );
%                     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%                     opts.Display = 'Off';
%                     opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
%                     opts.StartPoint = [-24.9070167618834 -729.888034703216 421.720294911524 -24.9070167618834 -308.167739791691 421.720294911524 -24.9070167618834 113.552555119833 421.720294911524 -24.9070167618834 535.272850031357 421.720294911524 -24.9070167618834 956.993144942882 421.720294911524 -24.9070167618834 1378.71343985441 421.720294911524];
%                     [fitresult, ~] = fit( x,y, fittype( 'gauss6' ), opts ); % Fit model to data.
%                     ppoz_curved = -fitresult(ppox); % flattened Z positions

                    ft = fittype( 'fourier1' );
                    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                    opts.Display = 'Off';
                    opts.StartPoint = [0 0 0 0.0643996589160774];
                    [fitresult, gof] = fit( x, y, ft, opts );                    
                    ppoz_curved = -fitresult(ppox); % flattened Z positions
                    
                    % surface projection (onto curved surface, probably only useful for plotting)
                    pos_flat_curve = [ppox,ppoy,ppoz_curved];
                    pdata.(part_now).pos_curve = single(pos_flat_curve);
                    pdata.(part_now).curve_gof = gof;
%                 end

%                 if ~isfield(pdata.(part_now),'pos_planar')
                    % planar projection (down onto flat surface)
                    pos_planar = [ppox,ppoy,ones(size(ppox))];
                    pdata.(part_now).pos_planar = single(pos_planar);
%                 end

%                 if ~isfield(pdata.(part_now),'pos_surficial')
                    % surficial projection (onto curve and then flattened)
                    newx = min(lx(:)) : 1 : max(lx(:)); % x points along edge of maze (1mm resolution)
                    newy = -fitresult(newx); % approximate y coordinates we expect on the maze along edge of maze
                    points = [newx(:) newy(:)];
                    ds = cumsum( [0; sqrt(sum((points(1:end-1,:)-points(2:end,:)).^2,2))] ); % calculate multidimensional distance from start of curve
                    ppox_flat = interp1(newx,ds,ppox,'linear','extrap')+min(lx(:)); % find new x coordinates based on distance along curve
                    ppoz_flat = ones(size(ppox)); % new z coordinates will just be bottom of maze
                    pos_surficial = [ppox_flat,ppoy,ppoz_flat];
                    pdata.(part_now).pos_surficial = single(pos_surficial);
                    
                    lx_surf = interp1(newx,ds,lx,'linear','extrap'); % find new x coordinates based on distance along curve
                    lx_surf = lx_surf + min(lx(:));
                    save(mname,'lx_surf','-mat','-append');   
                    pdata.maze_frame = single([lx(:) ly(:) lz(:) lx_surf(:)]);

                    % create figure
                    figure('visible','off','Units','pixels','Position',[10, 10, 850, 950]); % open/create figure
                    ax = axes('Units','pixels','Position',[50,690,700,260]);
                        plot(ppox,ppoz,'k'); hold on;
                        axis xy
                        daspect([1 1 1])
                        xi = min(ppox) : 1 : max(ppox);
                        yi = -fitresult(xi);
                        plot(xi,yi,'r');
                    
                    ax = axes('Units','pixels','Position',[50,420,700,260]);
                        ppox = ppox-min(lx);
                        ppox_flat = ppox_flat-min(lx_surf);
                        lx = lx-min(lx);
                        lx_surf = lx_surf-min(lx_surf);
                        xlims = [min([lx(:);lx_surf(:)])-20 max([lx(:);lx_surf(:)])+20];
                        ylims = [min([ly(:);ly(:)])-20 max([ly(:);ly(:)])+20];
                    
                        plot(ppox,ppoy,'k'); hold on;
                        plot(lx,ly,'r');
                        axis xy
                        ax.XLim = xlims;
                        ax.YLim = ylims;
                        daspect([1 1 1])                    
                        text(0,1.05,sprintf('Planar > Maze width: %.2f mm',range(lx)),'Units','normalized','FontSize',10,'HorizontalAlignment','left')
                    
                    ax = axes('Units','pixels','Position',[50,90,700,260]);
                        plot(ppox_flat,ppoy,'k'); hold on;
                        plot(lx_surf,ly,'r');
                        axis xy
                        ax.XLim = xlims;  
                        ax.YLim = ylims;                        
                        daspect([1 1 1])                     
                        text(0,1.05,sprintf('Surficial > Maze width: %.2f mm (%.1f x)',range(lx_surf),range(lx_surf)/range(lx)),'Units','normalized','FontSize',10,'HorizontalAlignment','left')
                    
                    fname = [pwd '\' pdata.outname '\maze_curve.png'];        
                    [~,~,~] = mkdir([pwd '\' pdata.outname '\']);
                    print(gcf,'-dpng','-r250',fname)
                    close(gcf); 
%                 end                   

            otherwise
                keyboard % do something
        end
    end
    disp(sprintf('\t...done'))    

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Run through clusters
    disp(sprintf('Analysing clusters...'))

    ucis = unique(sdata.uci); % list of unique cells in sdata    
    for uu = 1:length(ucis)
        uci = ucis{uu};
        disp(sprintf('\tCell %d of %d (%.f%%): %s',uu,length(ucis),uu/length(ucis)*100,uci))

        for pp = 1:nparts % for every part               
            part_now = part_config.part_names{pp};
            disp(sprintf('\t\t%s',part_now))                   
            
            idx = find( ismember(sdata.uci,uci) & sdata.partn==pp );

%% >>>>>>>>>> Planar firing rate maps            
            % planar (vertical projection) ratemap     
            % make this for every session
            disp(sprintf('\b >>| planar map'))                               
            rmset = mapset;
            pos = pdata.(part_now).pos_planar(:,1:2) .* 10; % in mm
            lx = pdata.maze_frame(:,1) .* 10; % in mm
            ly = pdata.maze_frame(:,2) .* 10; % in mm
            lz = pdata.maze_frame(:,3) .* 10; % in mm
            
            spk = pos(sdata.spt_pot_index{idx},:);
            rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
            rmset.srate = pdata.pos_srate;
            if isfield(pdata.(part_now),'planar_speedlift') % if a previously computed dwellmap exists, use this to save time
                speedlift = pdata.(part_now).planar_speedlift{1};
                [ratemap,dwellmap, ~, ~, ~] = rate_mapper(pos,spk,rmset,speedlift);             
            else
                [ratemap,dwellmap, ~, ~,speedlift] = rate_mapper(pos,spk,rmset);              
                pdata.(part_now).planar_speedlift = { single(speedlift) };
            end
            sdata.ratemap_planar(idx) = { single(ratemap) };
            pdata.(part_now).dwellmap_planar = single(dwellmap);
                
            % add dwell map to sdata
            if uu==1 % if this is the first cell of the session
                sdata.dwellmap_planar(idx) = { pdata.(part_now).dwellmap_planar }; % add surficial dwellmap            
            end            
            
%% >>>>>>>>>> Planar autocorrelation   
            disp(sprintf('\b | planar autocorr'))                               
            automap = ndautoCORR(ratemap,ratemap,50); % spatial autocorrelogram
            [~,gdata] = get_grid_score(automap,mapset.binsize/10,'method','allen'); % grid score and stats

            if ~any(ismember(sdata.Properties.VariableNames,'planar_amap')) % if the column(s) do not exist yet
                sdata.planar_amap = cell(size(sdata,1),1); % preallocate
            end
            sdata.planar_amap(idx,1) = { single(automap) };             
            if ~any(ismember(sdata.Properties.VariableNames,'planar_field_info')) % if the column(s) do not exist yet
                sdata.planar_field_info = NaN(size(sdata,1),6); % preallocate
            end
            sdata.planar_field_info(idx,:) = single([gdata.radius gdata.majaxislength gdata.minaxislength gdata.height gdata.width gdata.field_orientation]); 

%% >>>>>>>>>> Planar SI, rayleigh and grid score shuffle
            disp(sprintf('\b | planar shuffle'))                               
            pox = pdata.(part_now).pos_planar(:,1) .* 10; % in mm
            poy = pdata.(part_now).pos_planar(:,2) .* 10; % in mm
            pot = pdata.(part_now).pot;
            poh = pdata.(part_now).yaw; % in radians
            spt = sdata.spike_times{idx};            
            spindx = sdata.spt_pot_index{idx};
            sinfo = savelli_shuffle(pox,poy,pot,poh,spt,spindx,rmset,100);                     
            si = mean(sinfo.si_boot,'omitnan'); % spatial info (bootstrapped)
            si_z = (si - mean(sinfo.si_shuff,'omitnan')) ./ std(sinfo.si_shuff,'omitnan'); % spatial info (z-scored against shuffle)
            gs = mean(sinfo.g_boot,'omitnan'); % grid score (bootstrapped)
            gs_z = (gs - mean(sinfo.g_shuff,'omitnan')) ./ std(sinfo.g_shuff,'omitnan'); % grid score (z-scored against shuffle)
            rv = mean(sinfo.r_boot,'omitnan'); % rayleigh vector (bootstrapped)
            rv_z = (rv - mean(sinfo.r_shuff,'omitnan')) ./ std(sinfo.r_shuff,'omitnan'); % rayleigh vector (z-scored against shuffle)    
            if ~any(ismember(sdata.Properties.VariableNames,'planar_spatial_info_shuffles')) % if the column(s) do not exist yet
                sdata.planar_spatial_info_shuffles = NaN(size(sdata,1),6); % preallocate
            end               
            sdata.planar_spatial_info_shuffles(idx,:) = single([si si_z gs gs_z rv rv_z]);
            
%% >>>>>>>>>> Surficial firing rate maps                        
            % surficial (surface projection) ratemap    
            % make this for hilly sessions only
            if strcmp(part_now,'hills')   
                disp(sprintf('\b | surficial map'))                               
                
                rmset = mapset;
                pos = pdata.(part_now).pos_surficial(:,1:2) .* 10; % in mm
                lx = pdata.maze_frame(:,4) .* 10; % in mm
                ly = pdata.maze_frame(:,2) .* 10; % in mm
                lz = pdata.maze_frame(:,3) .* 10; % in mm

                spk = pos(sdata.spt_pot_index{idx},:);
                rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
                rmset.srate = pdata.pos_srate;              
                if isfield(pdata.(part_now),'surficial_speedlift') % if a previously computed dwellmap exists, use this to save time
                    speedlift = pdata.(part_now).surficial_speedlift{1};
                    [ratemap,dwellmap, ~, ~, ~] = rate_mapper(pos,spk,rmset,speedlift);                     
                else
                    [ratemap,dwellmap, ~, ~,speedlift] = rate_mapper(pos,spk,rmset);              
                    pdata.(part_now).surficial_speedlift = { single(speedlift) };
                end
                sdata.ratemap_surficial(idx) = { single(ratemap) };
                pdata.(part_now).dwellmap_surficial = single(dwellmap); 
                
                % add surficial positions to sdata
                if uu==1 % if this is the first cell of the session
                    if ~any(ismember(sdata.Properties.VariableNames,'pos_surficial')) % if the column does not exist yet
                        sdata.pos_surficial = cell(size(sdata,1),1); % preallocate
                    end
                    sdata.pos_surficial(idx) = { single(pdata.(part_now).pos_surficial) }; % add surficial positions
                    if ~any(ismember(sdata.Properties.VariableNames,'dwellmap_surficial')) % if the column does not exist yet
                        sdata.dwellmap_surficial = cell(size(sdata,1),1); % preallocate
                    end
                    sdata.dwellmap_surficial(idx) = { pdata.(part_now).dwellmap_surficial }; % add surficial dwellmap            
                end    
                
%% >>>>>>>>>> Surficial autocorrelation   
                disp(sprintf('\b | surficial autocorr'))                               

                automap = ndautoCORR(ratemap,ratemap,50); % spatial autocorrelogram
                [~,gdata] = get_grid_score(automap,mapset.binsize/10,'method','allen'); % grid score and stats

                if ~any(ismember(sdata.Properties.VariableNames,'surficial_amap')) % if the column(s) do not exist yet
                    sdata.surficial_amap = cell(size(sdata,1),1); % preallocate
                end
                sdata.surficial_amap(idx,1) = { single(automap) }; 
                if ~any(ismember(sdata.Properties.VariableNames,'surficial_field_info')) % if the column(s) do not exist yet
                    sdata.surficial_field_info = NaN(size(sdata,1),6); % preallocate
                end
                sdata.surficial_field_info(idx,:) = single([gdata.radius gdata.majaxislength gdata.minaxislength gdata.height gdata.width gdata.field_orientation]); 

%% >>>>>>>>>> Surficial SI, rayleigh and grid score shuffle
                disp(sprintf('\b | surficial shuffle'))                               

                pox = pdata.(part_now).pos_surficial(:,1) .* 10; % in mm
                poy = pdata.(part_now).pos_surficial(:,2) .* 10; % in mm
                pot = pdata.(part_now).pot;
                poh = pdata.(part_now).yaw; % in radians
                spt = sdata.spike_times{idx};            
                spindx = sdata.spt_pot_index{idx};
                sinfo = savelli_shuffle(pox,poy,pot,poh,spt,spindx,rmset,100);                  
                si = mean(sinfo.si_boot,'omitnan'); % spatial info (bootstrapped)
                si_z = (si - mean(sinfo.si_shuff,'omitnan')) ./ std(sinfo.si_shuff,'omitnan'); % spatial info (z-scored against shuffle)
                gs = mean(sinfo.g_boot,'omitnan'); % grid score (bootstrapped)
                gs_z = (gs - mean(sinfo.g_shuff,'omitnan')) ./ std(sinfo.g_shuff,'omitnan'); % grid score (z-scored against shuffle)
                rv = mean(sinfo.r_boot,'omitnan'); % rayleigh vector (bootstrapped)
                rv_z = (rv - mean(sinfo.r_shuff,'omitnan')) ./ std(sinfo.r_shuff,'omitnan'); % rayleigh vector (z-scored against shuffle)   
                if ~any(ismember(sdata.Properties.VariableNames,'surficial_spatial_info_shuffles')) % if the column(s) do not exist yet
                    sdata.surficial_spatial_info_shuffles = NaN(size(sdata,1),6); % preallocate
                end                  
                sdata.surficial_spatial_info_shuffles(idx,:) = single([si si_z gs gs_z rv rv_z]);   
            end
        end
        
%% >>>>>>>>>> Figure for this cell
        fname = [pwd '\' pdata.outname '\part_figures\' uci '_projections.png'];
        if ~save_figs % if there were no spikes in this part or we don't want figures at all
            continue % don't make the figure  
        elseif exist(fname,'file') && skipfigs % if the figure exists and we don't want to overwrite it
            continue % don't make the figure             
        else
            %putvar(pdata,sdata,uci)
            klustfig_hills(pdata,sdata,uci,fast_figures,'off');
        end 
    end

%% >>>>>>>>>> Plot all place cells
    for pp = 1:nparts
        part_now = part_config.part_names{pp};
        fname = [pwd '\' pdata.outname '\' part_now '_all_place_cells.png'];
        if ~save_figs % if there were no spikes in this part or we don't want figures at all
            continue % don't make the figure  
        elseif exist(fname,'file') && skipfigs % if the figure exists and we don't want to overwrite it
            continue % don't make the figure             
        else
            klustfig_all_pcells(pdata,sdata,pp,fast_figures,'off')
        end     
    end
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    sdata = rmprop(sdata,'pdata');
    sdata = addprop(sdata,{'pdata'},{'table'});
    sdata.Properties.CustomProperties.pdata = pdata;
    save([pwd '\' pdata.outname '\sdata.mat'],'sdata'); % save session data
    analysis_log({'hill_analysis'},1,'version',{'v6.0.0'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('hill_analysis has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',[pwd '\' pdata.outname],' &'');">','klustest folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');























