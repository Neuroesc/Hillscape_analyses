function [clumaa] = PIT_correlations(config,pidx,clumaa,posdata,overwrite)
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
    frate_cutoff = 1;
    iti = 1000; % N of shuffles
    if ~exist('overwrite','var') || isempty(overwrite) || isnan(overwrite)
        overwrite = 0;
    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>>>>>>>>>>>> Between session stability
    disp(sprintf('\t...between session correlations'))
    if ~any(ismember(clumaa.Properties.VariableNames,'between_session_stability')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.between_session_stability = NaN(size(clumaa,1),12);
        
        for ii = 1:length(ucis)
            % normal correlations
            m1 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1}; % arena 1
            m2 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==2}; % hills
            m2s = clumaa.ratemap_surficial{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==2}; % hills, surficial map                    
            m3 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==3}; % arena 2
            idx = ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1;
    
            clumaa.between_session_stability(idx,1) = get_correlation(m1,m2,frate_cutoff); % arena 1 vs hills
            clumaa.between_session_stability(idx,2) = get_correlation(m1,m3,frate_cutoff); % arena 1 vs arena 2
    
            % correlations limited to overlapping bins
            msk = get_map_mask(m1,'overlapping'); % sub function to get mask
            m1b = m1(msk);
            m2b = m2(msk);
            m3b = m3(msk);
            clumaa.between_session_stability(idx,3) = get_correlation(m1b,m2b,frate_cutoff); % arena 1 vs hills (troughs only)
            clumaa.between_session_stability(idx,4) = get_correlation(m1b,m3b,frate_cutoff); % arena 1 vs arena 2 (troughs only)

            % correlations limited to edges
            msk = get_map_mask(m1,'edges'); % sub function to get mask
            m1e = m1(msk);
            m2e = m2(msk);
            m3e = m3(msk);
            clumaa.between_session_stability(idx,5) = get_correlation(m1e,m2e,frate_cutoff); % arena 1 vs hills (edges only)
            clumaa.between_session_stability(idx,6) = get_correlation(m1e,m3e,frate_cutoff); % arena 1 vs arena 2 (edges only)

            % correlations limited to center
            msk = get_map_mask(m1,'center'); % sub function to get mask
            m1e = m1(msk);
            m2e = m2(msk);
            m3e = m3(msk);
            clumaa.between_session_stability(idx,7) = get_correlation(m1e,m2e,frate_cutoff); % arena 1 vs hills (center only)
            clumaa.between_session_stability(idx,8) = get_correlation(m1e,m3e,frate_cutoff); % arena 1 vs hills (center only)

            % correlations limited to corners
            msk = get_map_mask(m1,'corners'); % sub function to get mask
            m1c = m1(msk);
            m2c = m2(msk);
            m3c = m3(msk);
            clumaa.between_session_stability(idx,9) = get_correlation(m1c,m2c,frate_cutoff); % arena 1 vs hills (corners only)
            clumaa.between_session_stability(idx,10) = get_correlation(m1c,m3c,frate_cutoff); % arena 1 vs hills (corners only)

            % correlations aligned by maze ends
            m1a = zeros(size(m2s));
            m1a(:,1:size(m1,2)) = m1;
            m1b = zeros(size(m2s));
            m1b(:,(size(m1b,2)-size(m1,2))+1:end) = m1;
            clumaa.between_session_stability(idx,11) = get_correlation(m1a,m2s,frate_cutoff); % arena 1 vs hills (aligned left)
            clumaa.between_session_stability(idx,12) = get_correlation(m1b,m2s,frate_cutoff); % arena 1 vs hills (aligned right)
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>> Between session shuffles
    shuff = 1;
    rng(999); % for reproducibility    
    disp(sprintf('\t...between session shuffles'))
    if shuff
        fname = [config.data_out_dir 'PIT_correlations_between_shuffles.mat'];
        if ~exist(fname,'file') || overwrite
            ucis = unique(clumaa.uci(pidx));
            shuffs = NaN(iti,12); % preallocate
            for jj = 1:iti
                uix = randperm(numel(ucis),2); % 2 random place cells
    
                % normal correlations
                m1 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{uix(1)}) & clumaa.partn==1}; % arena 1, cell 1
                m2 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==1}; % arena 1, cell 2
                m3 = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==2}; % hills, cell 2

                shuffs(jj,1) = get_correlation(m1,m2,frate_cutoff); % arena 1 vs hills
                shuffs(jj,2) = get_correlation(m1,m3,frate_cutoff); % arena 1 vs arena 2
        
                % correlations limited to overlapping bins
                msk = get_map_mask(m1,'overlap'); % sub function to get mask
                m1b = m1(msk);
                m2b = m2(msk);
                m3b = m3(msk);
                shuffs(jj,3) = get_correlation(m1b,m2b,frate_cutoff); % arena 1 vs hills (troughs only)
                shuffs(jj,4) = get_correlation(m1b,m3b,frate_cutoff); % arena 1 vs arena 2 (troughs only)

                % correlations limited to edges
                msk = get_map_mask(m1,'edges'); % sub function to get mask
                m1e = m1(msk);
                m2e = m2(msk);
                m3e = m3(msk);
                shuffs(jj,5) = get_correlation(m1e,m2e,frate_cutoff); % arena 1 vs hills (edges only)
                shuffs(jj,6) = get_correlation(m1e,m3e,frate_cutoff); % arena 1 vs arena 2 (edges only)
    
                % correlations limited to center
                msk = get_map_mask(m1,'center'); % sub function to get mask
                m1e = m1(msk);
                m2e = m2(msk);
                m3e = m3(msk);
                shuffs(jj,7) = get_correlation(m1e,m2e,frate_cutoff); % arena 1 vs hills (center only)
                shuffs(jj,8) = get_correlation(m1e,m3e,frate_cutoff); % arena 1 vs arena 2 (center only)
    
                % correlations limited to corners
                msk = get_map_mask(m1,'corners'); % sub function to get mask
                m1c = m1(msk);
                m2c = m2(msk);
                m3c = m3(msk);
                shuffs(jj,9) = get_correlation(m1c,m2c,frate_cutoff); % arena 1 vs hills (corners only)
                shuffs(jj,10) = get_correlation(m1c,m3c,frate_cutoff); % arena 1 vs arena 2 (corners only)
    
                % correlations aligned by maze ends
                m1a = zeros(size(m2s));
                m1a(:,1:size(m1,2)) = m1;
                m1b = zeros(size(m2s));
                m1b(:,(size(m1b,2)-size(m1,2))+1:end) = m1;
                shuffs(jj,11) = get_correlation(m1a,m2s,frate_cutoff); % arena 1 vs hills (aligned left)
                shuffs(jj,12) = get_correlation(m1b,m2s,frate_cutoff); % arena 1 vs hills (aligned right)
            end
            save(fname,'shuffs');
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>> Within session stability
%% Make half session firing rate maps first
    % Map settings
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

    disp(sprintf('\t...within session half-maps'))
    if ~any(ismember(clumaa.Properties.VariableNames,'ratemap_planar_half')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.ratemap_planar_half = cell(size(clumaa,1),2); % preallocate
        clumaa.ratemap_surficial_half = cell(size(clumaa,1),2); % preallocate

        for ii = 1:length(ucis) % for every cell
            for pp = 1:3 % for every part               
                idx = find( ismember(clumaa.uci,ucis{ii}) & clumaa.partn==pp );
                pos_idx = clumaa.pos_idx(idx);
                pos = posdata.pos{posdata.pos_idx==pos_idx};
                mf = posdata.maze_frame{posdata.pos_idx==pos_idx}; 

                session_times = posdata.session_times{posdata.pos_idx==pos_idx};
                pidxn = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part

                mid_time = median(pos.pot(pidxn));
                pidx1 = pos.pot > session_times(pp,1) & pos.pot < mid_time; % index for position data in this part, first half
                pidx2 = pos.pot > mid_time & pos.pot < session_times(pp,2); % index for position data in this part, second half

                %% Planar firing rate maps            
                rmset = mapset;
                lx = mf(:,1); % in mm
                ly = mf(:,2); % in mm

                spt = clumaa.spike_times_s{idx};
                sidx1 = spt > session_times(pp,1) & spt < mid_time; % index for spike data in this part, first half
                sidx2 = spt > mid_time & spt < session_times(pp,2); % index for spike data in this part, second half
                spike_index = clumaa.spike_index{idx};

                % first half map
                pos_map = [pos.pox_planar(pidx1) pos.poy_planar(pidx1)];
                spk_map = [pos.pox_planar(spike_index(sidx1)) pos.poy_planar(spike_index(sidx1))];
                rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
                rmset.srate = 50;
                [ratemap1,~, ~, ~,~] = rate_mapper(pos_map,spk_map,rmset);              

                % second half map
                pos_map = [pos.pox_planar(pidx2) pos.poy_planar(pidx2)];
                spk_map = [pos.pox_planar(spike_index(sidx2)) pos.poy_planar(spike_index(sidx2))];
                [ratemap2,~, ~, ~,~] = rate_mapper(pos_map,spk_map,rmset);              

                % accumulate maps           
                clumaa.ratemap_planar_half(idx,:) = { single(ratemap1) single(ratemap2) };

                %% Surficial firing rate maps                        
                % surficial (surface projection) ratemap    
                if pp==2                  
                    rmset = mapset;
                    lx = mf(:,4); % in mm
                    ly = mf(:,2); % in mm
    
                    % first half map
                    pos_map = [pos.pox_surficial(pidx1) pos.pox_surficial(pidx1)];
                    spk_map = [pos.pox_surficial(spike_index(sidx1)) pos.pox_surficial(spike_index(sidx1))];
                    rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
                    rmset.srate = 50;
                    [ratemap1,~, ~, ~,~] = rate_mapper(pos_map,spk_map,rmset);              
        
                    % second half map
                    pos_map = [pos.pox_surficial(pidx2) pos.pox_surficial(pidx2)];
                    spk_map = [pos.pox_surficial(spike_index(sidx2)) pos.pox_surficial(spike_index(sidx2))];
                    [ratemap2,~, ~, ~,~] = rate_mapper(pos_map,spk_map,rmset);              
    
                    % accumulate maps               
                    clumaa.ratemap_surficial_half(idx,:) = { single(ratemap1) single(ratemap2) };
                end
            end
        end
    end
    
%% Calculate correlations
    disp(sprintf('\t...within session correlations'))
    if ~any(ismember(clumaa.Properties.VariableNames,'within_session_stability')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.within_session_stability = NaN(size(clumaa,1),6);
        
        for ii = 1:length(ucis) % for every cell
            for pp = 1:3 % for every part   
                idx = ismember(clumaa.uci,ucis{ii}) & clumaa.partn==pp;
                m1a = clumaa.ratemap_planar_half{idx,1}; % planar half 1
                m1b = clumaa.ratemap_planar_half{idx,2}; % planar half 2
        
                % normal correlation
                clumaa.within_session_stability(idx,1) = get_correlation(m1a,m1b,frate_cutoff); % arena 1 vs hills

                % correlations limited to overlapping bins
                msk = get_map_mask(m1a,'overlapping'); % sub function to get mask
                m2a = m1a(msk);
                m2b = m1b(msk);
                clumaa.within_session_stability(idx,3) = get_correlation(m2a,m2b,frate_cutoff); % arena 1 vs hills (troughs only)
    
                % correlations limited to edges
                msk = get_map_mask(m1a,'edges'); % sub function to get mask
                m3a = m1a(msk);
                m3b = m1b(msk);
                clumaa.within_session_stability(idx,4) = get_correlation(m3a,m3b,frate_cutoff); % arena 1 vs hills (troughs only)
    
                % correlations limited to center
                msk = get_map_mask(m1a,'center'); % sub function to get mask
                m4a = m1a(msk);
                m4b = m1b(msk);
                clumaa.within_session_stability(idx,5) = get_correlation(m4a,m4b,frate_cutoff); % arena 1 vs hills (troughs only)
    
                % correlations limited to corners
                msk = get_map_mask(m1a,'corners'); % sub function to get mask
                m5a = m1a(msk);
                m5b = m1b(msk);
                clumaa.within_session_stability(idx,6) = get_correlation(m5a,m5b,frate_cutoff); % arena 1 vs hills (troughs only)
            end
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>> Within session shuffles
    %% shuffles
    shuff = 1;
    rng(999); % for reproducibility
    disp(sprintf('\t...within session shuffles'))
    if shuff
        fname = [config.data_out_dir 'PIT_correlations_within_shuffles.mat'];
        if ~exist(fname,'file') || overwrite        
            ucis = unique(clumaa.uci(pidx));
            shuffs_within = NaN(iti,4,2); % preallocate
            for jj = 1:iti
                uix = randperm(numel(ucis),2); % 2 random place cells
    
                for pp = 1:2 % for each part
                    % 1st and 2nd half
                    m1a = clumaa.ratemap_planar_half{(ismember(clumaa.uci,ucis{uix(1)}) & clumaa.partn==pp),1}; % maze half 1, cell 1 
                    m1b = clumaa.ratemap_planar_half{(ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==pp),2}; % maze half 2, cell 2
                    shuffs_within(jj,1,pp) = get_correlation(m1a,m1b,frate_cutoff);
        
                    % correlations limited to edges
                    msk = get_map_mask(m1a,'edges'); % sub function to get mask
                    m1e = m1a(msk);
                    m2e = m1b(msk);
                    shuffs_within(jj,2,pp) = get_correlation(m1e,m2e,frate_cutoff); % arena 1 vs hills (troughs only)
        
                    % correlations away from edges
                    msk = get_map_mask(m1a,'center'); % sub function to get mask
                    m1m = m1a(msk);
                    m2m = m1b(msk);
                    shuffs_within(jj,3,pp) = get_correlation(m1m,m2m,frate_cutoff); % arena 1 vs hills (troughs only)
        
                    % correlations limited to corners
                    msk = get_map_mask(m1a,'corners'); % sub function to get mask
                    m1c = m1a(msk);
                    m2c = m1b(msk);
                    shuffs_within(jj,4,pp) = get_correlation(m1c,m2c,frate_cutoff); % arena 1 vs hills (troughs only)
                end
            end
            save(fname,'shuffs_within');
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>> Pitch movement remapping
%% Make half session firing rate maps first
    % Map settings
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

    disp(sprintf('\t...pitch half maps'))
    if ~any(ismember(clumaa.Properties.VariableNames,'ratemap_pitch_half')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.ratemap_pitch_half = cell(size(clumaa,1),2); % preallocate

        loopout = looper(length(ucis));
        for ii = 1:length(ucis) % for every cell
            for pp = 1:3 % for every part     
                idx = find( ismember(clumaa.uci,ucis{ii}) & clumaa.partn==pp );
                pos_idx = clumaa.pos_idx(idx);
                pos = posdata.pos{posdata.pos_idx==pos_idx};
                mf = posdata.maze_frame{posdata.pos_idx==pos_idx}; 

                session_times = posdata.session_times{posdata.pos_idx==pos_idx};
                pidxn = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part

                %% Planar firing rate maps            
                rmset = mapset;
                lx = mf(:,1); % in mm
                ly = mf(:,2); % in mm

                spt = clumaa.spike_times_s{idx};   
                sidx = clumaa.spike_index{idx}(spt > session_times(pp,1) & spt < session_times(pp,2)); % index for spikes in this part              
                spike_pitch = pos.pitch(sidx);

                % first half map
                pos_map1 = [pos.pox_planar(pidxn & pos.pitch<0) pos.poy_planar(pidxn & pos.pitch<0)];
                spk_map1 = [pos.pox_planar(sidx(spike_pitch<0)) pos.poy_planar(sidx(spike_pitch<0))];
                rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
                rmset.srate = 50;
                [ratemap1,~, ~, ~,~] = rate_mapper(pos_map1,spk_map1,rmset);              

                % second half map
                pos_map2 = [pos.pox_planar(pidxn & pos.pitch>0) pos.poy_planar(pidxn & pos.pitch>0)];
                spk_map2 = [pos.pox_planar(sidx(spike_pitch>0)) pos.poy_planar(sidx(spike_pitch>0))];
                [ratemap2,~, ~, ~,~] = rate_mapper(pos_map2,spk_map2,rmset);              

                % accumulate maps           
                clumaa.ratemap_pitch_half(idx,:) = { single(ratemap1) single(ratemap2) };

                % figure
                % subplot(3,2,1)
                % plot(pos.pox_planar(pidxn),pos.poy_planar(pidxn),'k'); hold on;
                % plot(pos.pox_planar(sidx),pos.poy_planar(sidx),'r.','MarkerSize',20); hold on;
                % daspect([1 1 1])
                % 
                % subplot(3,2,2)
                % imagesc(clumaa.ratemap_planar{idx})
                % daspect([1 1 1])
                % 
                % subplot(3,2,3)
                % plot(pos_map1(:,1),pos_map1(:,2),'Color',[.5 .5 .5],'Marker','o'); hold on;
                % plot(spk_map1(:,1),spk_map1(:,2),'r.','MarkerSize',20); hold on;
                % daspect([1 1 1])
                % 
                % subplot(3,2,4)
                % imagesc(ratemap1)
                % daspect([1 1 1])
                % 
                % subplot(3,2,5)
                % plot(pos_map2(:,1),pos_map2(:,2),'Color',[.5 .5 .5],'Marker','o'); hold on;
                % plot(spk_map2(:,1),spk_map2(:,2),'r.','MarkerSize',20); hold on;
                % daspect([1 1 1])
                % 
                % subplot(3,2,6)
                % imagesc(ratemap2)
                % daspect([1 1 1])
                % 
                % keyboard
            end
            loopout = looper(loopout);            
        end
    end
    
%% Calculate pitch correlations
    disp(sprintf('\t...pitch map correlations'))
    if ~any(ismember(clumaa.Properties.VariableNames,'pitch_stability')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.pitch_stability = NaN(size(clumaa,1),1);
        
        loopout = looper(length(ucis));        
        for ii = 1:length(ucis) % for every cell
            for pp = 1:3 % for every part   
                idx = ismember(clumaa.uci,ucis{ii}) & clumaa.partn==pp;
                m1a = clumaa.ratemap_pitch_half{idx,1}; % planar half 1
                m1b = clumaa.ratemap_pitch_half{idx,2}; % planar half 2
        
                % normal correlation
                clumaa.pitch_stability(idx,1) = get_correlation(m1a,m1b,frate_cutoff); % arena 1 vs hills
            end
            loopout = looper(loopout);                        
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>> pitch shuffles
    %% shuffles
    shuff = 1;
    rng(999); % for reproducibility
    disp(sprintf('\t...pitch shuffles'))
    if shuff
        fname = [config.data_out_dir 'PIT_correlations_pitch_shuffles.mat'];
        if ~exist(fname,'file') || overwrite        
            ucis = unique(clumaa.uci(pidx));
            shuffs_pitch = NaN(iti,2); % preallocate
            for jj = 1:iti
                uix = randperm(numel(ucis),2); % 2 random place cells
    
                for pp = 1:2 % for each part
                    % pitch up and down maps
                    m1a = clumaa.ratemap_pitch_half{(ismember(clumaa.uci,ucis{uix(1)}) & clumaa.partn==pp),1}; % pitch up map, cell 1 
                    m1b = clumaa.ratemap_pitch_half{(ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==pp),2}; % pitch down map, cell 2
                    shuffs_pitch(jj,pp) = get_correlation(m1a,m1b,frate_cutoff);
                end
            end
            save(fname,'shuffs_pitch');
        end
    end
    disp(sprintf('\t...done'))


%% >>>>>>>>>> Azimuth movement remapping
%% Make half session firing rate maps first
    % Map settings
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

    disp(sprintf('\t...azimuth half maps'))
    if ~any(ismember(clumaa.Properties.VariableNames,'ratemap_azimuth_half')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.ratemap_azimuth_half = cell(size(clumaa,1),2); % preallocate

        loopout = looper(length(ucis));
        for ii = 1:length(ucis) % for every cell
            for pp = 1:3 % for every part     
                idx = find( ismember(clumaa.uci,ucis{ii}) & clumaa.partn==pp );
                pos_idx = clumaa.pos_idx(idx);
                pos = posdata.pos{posdata.pos_idx==pos_idx};
                mf = posdata.maze_frame{posdata.pos_idx==pos_idx}; 

                session_times = posdata.session_times{posdata.pos_idx==pos_idx};
                pidxn = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part

                %% Planar firing rate maps            
                rmset = mapset;
                lx = mf(:,1); % in mm
                ly = mf(:,2); % in mm

                spt = clumaa.spike_times_s{idx};   
                sidx = clumaa.spike_index{idx}(spt > session_times(pp,1) & spt < session_times(pp,2)); % index for spikes in this part              
                spike_azimuth = pos.yaw(sidx);

                % first half map
                idx1 = abs(pos.yaw)>90;
                idx2 = ~idx1;  
                sidx1 = abs(spike_azimuth)>90;
                sidx2 = abs(spike_azimuth)<90;                 
                pos_map1 = [pos.pox_planar(pidxn & idx1) pos.poy_planar(pidxn & idx1)];
                spk_map1 = [pos.pox_planar(sidx(sidx1)) pos.poy_planar(sidx(sidx1))];
                rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
                rmset.srate = 50;
                [ratemap1,~, ~, ~,~] = rate_mapper(pos_map1,spk_map1,rmset);              

                % second half map
                pos_map2 = [pos.pox_planar(pidxn & idx2) pos.poy_planar(pidxn & idx2)];
                spk_map2 = [pos.pox_planar(sidx(sidx2)) pos.poy_planar(sidx(sidx2))];
                [ratemap2,~, ~, ~,~] = rate_mapper(pos_map2,spk_map2,rmset);              

                % accumulate maps           
                clumaa.ratemap_azimuth_half(idx,:) = { single(ratemap1) single(ratemap2) };

% figure
% subplot(1,3,1)
% cline(pos.pox_planar(pidxn),pos.poy_planar(pidxn),pos.yaw(pidxn))
% colorbar
% set(gca,'CLim',[-180 180])
% c = cmocean('balance');
% colormap(gca,[c;flipud(c)])
% 
% subplot(1,3,2)
% cline(pos.pox_planar(pidxn & idx1),pos.poy_planar(pidxn & idx1),pos.yaw(pidxn & idx1))
% colorbar
% set(gca,'CLim',[-180 180])
% colormap(gca,[c;flipud(c)])
% 
% subplot(1,3,3)
% cline(pos.pox_planar(pidxn & idx2),pos.poy_planar(pidxn & idx2),pos.yaw(pidxn & idx2))
% colorbar
% set(gca,'CLim',[-180 180])
% colormap(gca,[c;flipud(c)])

                % figure
                % subplot(3,2,1)
                % plot(pos.pox_planar(pidxn),pos.poy_planar(pidxn),'k'); hold on;
                % plot(pos.pox_planar(sidx),pos.poy_planar(sidx),'r.','MarkerSize',20); hold on;
                % daspect([1 1 1])
                % 
                % subplot(3,2,2)
                % imagesc(clumaa.ratemap_planar{idx})
                % daspect([1 1 1])
                % 
                % subplot(3,2,3)
                % plot(pos_map1(:,1),pos_map1(:,2),'Color',[.5 .5 .5],'Marker','.'); hold on;
                % % plot(spk_map1(:,1),spk_map1(:,2),'r.','MarkerSize',20); hold on;
                % daspect([1 1 1])
                % 
                % subplot(3,2,4)
                % imagesc(ratemap1)
                % daspect([1 1 1])
                % 
                % subplot(3,2,5)
                % plot(pos_map2(:,1),pos_map2(:,2),'Color',[.5 .5 .5],'Marker','.'); hold on;
                % % plot(spk_map2(:,1),spk_map2(:,2),'r.','MarkerSize',20); hold on;
                % daspect([1 1 1])
                % 
                % subplot(3,2,6)
                % imagesc(ratemap2)
                % daspect([1 1 1])
                % 
                % keyboard
            end
            loopout = looper(loopout);            
        end
    end
    
%% Calculate azimuth correlations
    disp(sprintf('\t...azimuth map correlations'))
    if ~any(ismember(clumaa.Properties.VariableNames,'azimuth_stability')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.azimuth_stability = NaN(size(clumaa,1),1);
        
        loopout = looper(length(ucis));        
        for ii = 1:length(ucis) % for every cell
            for pp = 1:3 % for every part   
                idx = ismember(clumaa.uci,ucis{ii}) & clumaa.partn==pp;
                m1a = clumaa.ratemap_azimuth_half{idx,1}; % planar half 1
                m1b = clumaa.ratemap_azimuth_half{idx,2}; % planar half 2
        
                % normal correlation
                clumaa.azimuth_stability(idx,1) = get_correlation(m1a,m1b,frate_cutoff); % arena 1 vs hills
            end
            loopout = looper(loopout);                        
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>> azimuth shuffles
    %% shuffles
    shuff = 1;
    rng(999); % for reproducibility
    disp(sprintf('\t...azimuth shuffles'))
    if shuff
        fname = [config.data_out_dir 'PIT_correlations_azimuth_shuffles.mat'];
        if ~exist(fname,'file') || overwrite        
            ucis = unique(clumaa.uci(pidx));
            shuffs_azimuth = NaN(iti,2); % preallocate
            for jj = 1:iti
                uix = randperm(numel(ucis),2); % 2 random place cells
    
                for pp = 1:2 % for each part
                    % pitch up and down maps
                    m1a = clumaa.ratemap_azimuth_half{(ismember(clumaa.uci,ucis{uix(1)}) & clumaa.partn==pp),1}; % pitch up map, cell 1 
                    m1b = clumaa.ratemap_azimuth_half{(ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==pp),2}; % pitch down map, cell 2
                    shuffs_azimuth(jj,pp) = get_correlation(m1a,m1b,frate_cutoff);
                end
            end
            save(fname,'shuffs_azimuth');
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>>>>>>>>>>>> 3D HD stability
    disp(sprintf('\t...3DHD correlations'))
    if ~any(ismember(clumaa.Properties.VariableNames,'hd3d_stability')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.hd3d_stability = NaN(size(clumaa,1),2);
        
        for ii = 1:length(ucis)
            % normal correlations
            m1 = clumaa.hd_3d_ratemap{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1}; % arena 1
            m2 = clumaa.hd_3d_ratemap{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==2}; % hills
            m3 = clumaa.hd_3d_ratemap{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==3}; % arena 2
            idx = ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1;
    
            clumaa.hd3d_stability(idx,1) = get_correlation(m1,m2,frate_cutoff); % arena 1 vs hills
            clumaa.hd3d_stability(idx,2) = get_correlation(m1,m3,frate_cutoff); % arena 1 vs arena 2
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>> 3D HD shuffles
    shuff = 1;
    rng(999); % for reproducibility    
    disp(sprintf('\t...3dHD shuffles'))
    if shuff
        fname = [config.data_out_dir 'PIT_correlations_3dhd_shuffles.mat'];
        if ~exist(fname,'file') || overwrite
            ucis = unique(clumaa.uci(pidx));
            shuffs_3dhd = NaN(iti,2); % preallocate
            for jj = 1:iti
                uix = randperm(numel(ucis),2); % 2 random place cells
    
                % normal correlations
                m1 = clumaa.hd_3d_ratemap{ismember(clumaa.uci,ucis{uix(1)}) & clumaa.partn==1}; % arena 1, cell 1
                m2 = clumaa.hd_3d_ratemap{ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==1}; % arena 1, cell 2
                m3 = clumaa.hd_3d_ratemap{ismember(clumaa.uci,ucis{uix(2)}) & clumaa.partn==2}; % hills, cell 2

                shuffs_3dhd(jj,1) = get_correlation(m1,m2,frate_cutoff); % arena 1 vs hills
                shuffs_3dhd(jj,2) = get_correlation(m1,m3,frate_cutoff); % arena 1 vs arena 2
            end
            save(fname,'shuffs_3dhd');
        end
    end
    disp(sprintf('\t...done'))

%% >>>>>>>>>>>>>>>>>>>> Reflected map similarity
    disp(sprintf('\t...reflected correlations'))
    if ~any(ismember(clumaa.Properties.VariableNames,'reflected_stability')) || overwrite
        ucis = unique(clumaa.uci);

        % preallocate
        clumaa.reflected_stability = NaN(size(clumaa,1),2);
        
        for ii = 1:length(ucis)
            m1a = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1}; % arena 1
            m2a = fliplr(m1a); % arena 1, flipped horizontally
            m1b = clumaa.ratemap_planar{ismember(clumaa.uci,ucis{ii}) & clumaa.partn==2}; % hills
            m2b = fliplr(m1b); % hills, flipped horizontally
            
            idx = ismember(clumaa.uci,ucis{ii}) & clumaa.partn==1;
    
            clumaa.reflected_stability(idx,1) = get_correlation(m1a,m2a,frate_cutoff); % arena 1
            clumaa.reflected_stability(idx,2) = get_correlation(m1b,m2b,frate_cutoff); % hills

        end
    end
    disp(sprintf('\t...done'))

end % end function

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Sub function for map correlations
function r = get_correlation(m1,m2,frate_cutoff)
    r = NaN;
    if ~isempty(m1) && ~isempty(m2)
        if max(m1(:),[],'omitnan')>=frate_cutoff || max(m2(:),[],'omitnan')>=frate_cutoff       
            r = corr(m1(:),m2(:),'rows','pairwise','type','Pearson'); % arena 1 vs hills (troughs only)
        else
            r = NaN;
        end
    end
end























