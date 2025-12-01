function PIT_behaviour(data_dir,rnow,dnow)
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
    
    %% Map settings    
    pdata.mapset.srate = 50;
    mapset = pdata.mapset;
    pos = pdata.pos;

    % settings for anisotropy calculation
    dist_cutoff = 120;
    dist_cutoff_bins = dist_cutoff ./ mapset.binsize;
    min_run_speed = 5; % cm/s

%% >>>>>>>>>> Behaviour anisotropy
    session_times = pdata.session_times;
    pos = pdata.pos;

    for pp = 1:size(session_times,1)      
        pidx = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part
                    
        % create dummy ratemap just to get coordinate conversion
        rmset = mapset;
        spk = [NaN NaN];
        if pp==2
            p3 = [pos.pox_surficial(pidx) pos.poy_surficial(pidx)];  
            lx = pdata.maze_frame(:,4); % in mm
        else
            p3 = [pos.pox_planar(pidx) pos.poy_planar(pidx)];   
            lx = pdata.maze_frame(:,1); % in mm                
        end
        ly = pdata.maze_frame(:,2); % in mm        
        rmset.maplims = [min(lx,[],'omitnan') min(ly,[],'omitnan') max(lx,[],'omitnan') max(ly,[],'omitnan')]; % in mm
        
        % ratemap
        [~,dwellmap, ~, outputs,~] = rate_mapper(p3,spk,rmset);     

        % running speed & direction
        p3b = circshift(p3,-1,1); % shift data one time point backward
        [podn,ds] = cart2pol(p3(:,1)-p3b(:,1),p3(:,2)-p3b(:,2)); 
        podn = rad2deg(podn); % movement direction                                
        ds(end) = NaN; % last sample cannot be computed
        wsize = 3;
        dist_sum = movsum(ds,wsize,'omitnan','Endpoints','shrink'); % moving sum of distance over window of wsize samples
        tsize = wsize .* (1/mapset.srate); % the length (s) of each window
        povn = (dist_sum./10) ./ tsize; % speed in cm per second = distance in cm / time in seconds

        % map positions to dwellmap
        mapXY = outputs.points_to_map(p3);

        % overall anisotropy score
        mov_y_idx = (podn>45 & podn<135) | (podn<-45 & podn>-135);
        mov_x_idx = (podn>-45 & podn<45) | podn>135 | podn<-135;
        time_y = sum( povn>min_run_speed & mov_y_idx );
        time_x = sum( povn>min_run_speed & mov_x_idx );
        anisotropy_score = (time_y - time_x) ./ (time_y + time_x);
        pdata.(part_names{pp}).anisotropy_score = single(anisotropy_score);

        % anisotropy map
        anisotropy_map = NaN(size(dwellmap));
        for bb = 1:numel(dwellmap)
            [r,c] = ind2sub(size(dwellmap),bb);
            box_idx = mapXY(:,1)>(c-dist_cutoff_bins) & mapXY(:,1)<(c+dist_cutoff_bins) & mapXY(:,2)>(r-dist_cutoff_bins) & mapXY(:,2)<(r+dist_cutoff_bins) & povn>min_run_speed;

            mov_y_idx = (podn>45 & podn<135) | (podn<-45 & podn>-135);
            mov_x_idx = (podn>-45 & podn<45) | podn>135 | podn<-135;

            time_y = sum( box_idx & mov_y_idx );
            time_x = sum( box_idx & mov_x_idx );
            anisotropy_score = (time_y - time_x) ./ (time_y + time_x);
            anisotropy_map(bb) = anisotropy_score;

            if 0
                figure
                subplot(2,2,1)
                imagesc(dwellmap); hold on;
                axis xy
                daspect([1 1 1])
                plot(mapXY(:,1),mapXY(:,2),'w')
                
                subplot(2,2,2)
                imagesc(dwellmap); hold on;
                axis xy
                daspect([1 1 1])
                plot(mapXY(:,1),mapXY(:,2),'w')
                plot(mapXY(box_idx,1),mapXY(box_idx,2),'r')
                
                subplot(2,2,3)
                plot(mapXY(box_idx,1),mapXY(box_idx,2),'r')
                axis xy
                daspect([1 1 1])
                
                subplot(2,2,4)
                plot_now = [mapXY podn];
                plot_now(~box_idx,:) = repmat([NaN NaN NaN],sum(~box_idx),1);
                cline(plot_now(:,1),plot_now(:,2),podn(:));
                axis xy
                daspect([1 1 1])
                keyboard
            end
        end
        anisotropy_map = imgaussfilt(anisotropy_map,1);
        pdata.(part_names{pp}).anisotropy_map = single(anisotropy_map);

        if 0
            figure
            subplot(2,2,1)
            imagesc(dwellmap); hold on;
            axis xy
            daspect([1 1 1])
            plot(mapXY(:,1),mapXY(:,2),'w')
            
            subplot(2,2,2)
            cline(mapXY(:,1),mapXY(:,2),podn)
            axis xy
            daspect([1 1 1])
            
            subplot(2,2,3)
            imagesc(anisotropy_map); hold on;
            axis xy
            daspect([1 1 1])
            colorbar
            ax = gca;
            ax.CLim = [-1 1];
            plot(mapXY(:,1),mapXY(:,2),'Color',[.5 .5 .5 .5])
            
            subplot(2,2,4)
            imagesc(imgaussfilt(anisotropy_map,1)); hold on;
            axis xy
            daspect([1 1 1])
            colorbar
            ax = gca;
            ax.CLim = [-1 1];
        end
    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    cluma = rmprop(cluma,'pdata');
    cluma = addprop(cluma,{'pdata'},{'table'});
    cluma.Properties.CustomProperties.pdata = pdata;
    cname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    save(cname,'cluma'); % save session data
    analysis_log({'PIT_behaviour'},1,'version',{'v1.0.0'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('%s has finished. It took %0.3f seconds or %0.3f minutes',stk.name,toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');








 


















