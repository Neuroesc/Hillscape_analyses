function hill_analysis_v2(data_dir,rnow,dnow)
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
    
    override.maze_frame = 0; % 1 = do not load a maze frame from file
    skipfigs = 1; % 1 = skip making a figure if it exists already    
    save_figs = 1; % 1 = make and save figures (VERY time consuming, takes about 90% of klustest's time)
    fast_figures = 1; % 1 = 4-5x faster figure saving, but at a slightly lower quality
    frate_cutoff = 1;
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    stk = dbstack;
    tnow = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',stk.name,tnow))
    
%% >>>>>>>>>> Prepare the data
    sname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    disp(sprintf('Loading cluma: %s',sname))    
    load(sname,'cluma'); % load saved session data
    pdata = cluma.Properties.CustomProperties.pdata;
    part_names = {'arena1','hills','arena2'};
    nparts = size(pdata.session_times,1);
    disp(sprintf('\t...%d sessions',nparts));   
    disp(sprintf('\t...recording date: %s',dnow));            
    disp(sprintf('\t...done'))     
    
%% >>>>>>>>>> Add maze boundaries
    % disp(sprintf('Calculating maze boundaries...'))    
    % mname = [data_dir '\' rnow '\' dnow '\3Danalysis\maze_boundaries.mat']; 
    % fname = [data_dir '\' rnow '\' dnow '\3Danalysis\maze_boundaries.png'];       
    % if exist(mname,'file') && ~override.maze_frame && 1
    %     disp(sprintf('\t\t...loading from file'))                                    
    %     load(mname,'lx','ly','lz','-mat');             
    % else
    %     [lx,ly,lz] = fit_hill_frame_v2(pdata,0,mname,fname); % fit the maze frame to the data for all 3 parts simultaneously    
    % end
    disp(sprintf('\t...done'))    
    pos = pdata.pos;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Run through clusters
    disp(sprintf('Analysing clusters...'))

    ucis = unique(cluma.uci); % list of unique cells in sdata    
    for uu = 1:length(ucis)
        uci = ucis{uu};
        disp(sprintf('\tCell %d of %d (%.f%%): %s',uu,length(ucis),uu/length(ucis)*100,uci))

        for pp = 1:nparts % for every part               
            part_now = part_names{pp};
            disp(sprintf('\t\t%s',part_now))                   
            
            idx = find( ismember(cluma.uci,uci) & cluma.partn==pp );
            session_times = pdata.session_times;
            pidx = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part

            mid_time = median(pos.pot(pidx));
            pidx1 = pos.pot > session_times(pp,1) & pos.pot < mid_time; % index for position data in this part, first half
            pidx2 = pos.pot > mid_time & pos.pot < session_times(pp,2); % index for position data in this part, second half

%% >>>>>>>>>> Planar firing rate maps            
            % planar (vertical projection) ratemap     
            % make this for every session
            disp(sprintf('\b >>| planar map'))                               
            rmset = mapset;
            lx = pdata.maze_frame(:,1); % in mm
            ly = pdata.maze_frame(:,2); % in mm
            lz = pdata.maze_frame(:,3); % in mm

            spt = cluma.spike_times_s{idx};
            sidx1 = spt > session_times(pp,1) & spt < mid_time; % index for spike data in this part, first half
            sidx2 = spt > mid_time & spt < session_times(pp,2); % index for spike data in this part, second half
            spike_index = cluma.spike_index{idx};

            % first half map
            pos_map = [pos.pox_planar(pidx1) pos.poy_planar(pidx1)];
            spk_map = [pos.pox_planar(spike_index(sidx1)) pos.poy_planar(spike_index(sidx1))];
            rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
            rmset.srate = 50;
            if isfield(pdata,part_now) && isfield(pdata.(part_now),'planar_speedlift_half1') % if a previously computed dwellmap exists, use this to save time
                speedlift = pdata.(part_now).planar_speedlift_half1{1};
                [ratemap1,~, ~, ~, ~] = rate_mapper(pos_map,spk_map,rmset,speedlift);             
            else
                [ratemap1,~, ~, ~,speedlift] = rate_mapper(pos_map,spk_map,rmset);              
                pdata.(part_now).planar_speedlift_half1 = { single(speedlift) };
            end

            % second half map
            pos_map = [pos.pox_planar(pidx2) pos.poy_planar(pidx2)];
            spk_map = [pos.pox_planar(spike_index(sidx2)) pos.poy_planar(spike_index(sidx2))];
            if isfield(pdata,part_now) && isfield(pdata.(part_now),'planar_speedlift_half2') % if a previously computed dwellmap exists, use this to save time
                speedlift = pdata.(part_now).planar_speedlift_half2{1};
                [ratemap2,~, ~, ~, ~] = rate_mapper(pos_map,spk_map,rmset,speedlift);             
            else
                [ratemap2,~, ~, ~,speedlift] = rate_mapper(pos_map,spk_map,rmset);              
                pdata.(part_now).planar_speedlift_half2 = { single(speedlift) };
            end
            if ~any(ismember(cluma.Properties.VariableNames,'ratemap_planar_half')) % if the column(s) do not exist yet
                cluma.ratemap_planar_half = cell(size(cluma,1),2); % preallocate
            end            
            cluma.ratemap_planar_half(idx,:) = { single(ratemap1) single(ratemap2) };

            % correlation
            if ~any(ismember(cluma.Properties.VariableNames,'within_stability')) % if the column(s) do not exist yet
                cluma.within_stability = NaN(size(cluma,1),5); % preallocate
            end    
            if max(ratemap1(:),[],'omitnan')>=frate_cutoff || max(ratemap2(:),[],'omitnan')>=frate_cutoff
                cluma.within_stability(idx,1) = corr(ratemap1(:),ratemap2(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills
            end       
         

            %% Extra within session correlations
            % correlations limited to edges
            e_idx = zeros(size(ratemap1));
            e_idx(:,1) = 1;
            e_idx(:,end) = 1;
            e_idx(1,:) = 1;
            e_idx(end,:) = 1;
            e_dist = bwdist(e_idx);
            e_log = e_dist<5;

            m1e = ratemap1(e_log);
            m2e = ratemap2(e_log);
            if max(m1e(:),[],'omitnan')>=frate_cutoff || max(m2e(:),[],'omitnan')>=frate_cutoff            
                cluma.within_stability(idx,3) = corr(m1e(:),m2e(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills (troughs only)
            end

            % correlations away from edges
            e_log = e_dist>10;
            m1m = ratemap1(e_log);
            m2m = ratemap2(e_log);
            if max(m1m(:),[],'omitnan')>=frate_cutoff || max(m2m(:),[],'omitnan')>=frate_cutoff            
                cluma.within_stability(idx,4) = corr(m1m(:),m2m(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills (troughs only)
            end

            % correlations limited to corners
            c_idx = zeros(size(ratemap1));
            c_idx(1,1) = 1;
            c_idx(1,end) = 1;
            c_idx(end,1) = 1;
            c_idx(end,end) = 1;
            c_dist = bwdist(c_idx);
            c_log = c_dist<15;

            m1c = ratemap1(c_log);
            m2c = ratemap2(c_log);
            if max(m1c(:),[],'omitnan')>=frate_cutoff || max(m2c(:),[],'omitnan')>=frate_cutoff            
                cluma.within_stability(idx,5) = corr(m1c(:),m2c(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills (troughs only)
            end

            if 0
                figure
                subplot(1,2,1)
                imagesc(ratemap1)
                daspect([1 1 1])
                
                subplot(1,2,2)
                imagesc(ratemap2)
                daspect([1 1 1])
                keyboard
                close(gcf)
            end

%% >>>>>>>>>> Surficial firing rate maps                        
            % surficial (surface projection) ratemap    
            % make this for hilly sessions only
            if strcmp(part_now,'hills')   
                disp(sprintf('\b | surficial map'))                               
                
                rmset = mapset;
                lx = pdata.maze_frame(:,4); % in mm
                ly = pdata.maze_frame(:,2); % in mm
                lz = pdata.maze_frame(:,3); % in mm

                % first half map
                pos_map = [pos.pox_surficial(pidx1) pos.pox_surficial(pidx1)];
                spk_map = [pos.pox_surficial(spike_index(sidx1)) pos.pox_surficial(spike_index(sidx1))];
                rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
                rmset.srate = 50;
                if isfield(pdata,part_now) && isfield(pdata.(part_now),'surficial_speedlift_half1') % if a previously computed dwellmap exists, use this to save time
                    speedlift = pdata.(part_now).surficial_speedlift_half1{1};
                    [ratemap1,~, ~, ~, ~] = rate_mapper(pos_map,spk_map,rmset,speedlift);             
                else
                    [ratemap1,~, ~, ~,speedlift] = rate_mapper(pos_map,spk_map,rmset);              
                    pdata.(part_now).surficial_speedlift_half1 = { single(speedlift) };
                end
    
                % second half map
                pos_map = [pos.pox_surficial(pidx2) pos.pox_surficial(pidx2)];
                spk_map = [pos.pox_surficial(spike_index(sidx2)) pos.pox_surficial(spike_index(sidx2))];
                if isfield(pdata,part_now) && isfield(pdata.(part_now),'surficial_speedlift_half2') % if a previously computed dwellmap exists, use this to save time
                    speedlift = pdata.(part_now).surficial_speedlift_half2{1};
                    [ratemap2,~, ~, ~, ~] = rate_mapper(pos_map,spk_map,rmset,speedlift);             
                else
                    [ratemap2,~, ~, ~,speedlift] = rate_mapper(pos_map,spk_map,rmset);              
                    pdata.(part_now).surficial_speedlift_half2 = { single(speedlift) };
                end
                if ~any(ismember(cluma.Properties.VariableNames,'ratemap_surficial_half')) % if the column(s) do not exist yet
                    cluma.ratemap_surficial_half = cell(size(cluma,1),2); % preallocate
                end                 
                cluma.ratemap_surficial_half(idx,:) = { single(ratemap1) single(ratemap2) };
    
                % correlation
                if max(ratemap1(:),[],'omitnan')>=frate_cutoff || max(ratemap2(:),[],'omitnan')>=frate_cutoff
                    cluma.within_stability(idx,2) = corr(ratemap1(:),ratemap2(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills
                end                    

            end
        end
    end
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    cluma = rmprop(cluma,'pdata');
    cluma = addprop(cluma,{'pdata'},{'table'});
    pdata.pos = pos;
    pdata.mapset = mapset;    
    cluma.Properties.CustomProperties.pdata = pdata;
    cname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    save(cname,'cluma'); % save session data
    analysis_log({'stability_analysis'},1,'version',{'v1.1.0'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('stability_analysis has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');























