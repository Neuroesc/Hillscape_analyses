function hill_repeat_analysis_v2(data_dir,rnow,dnow)
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
% version 1.0.0, Release 21/02/23 Initial release

%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK  
    outname = 'klustest';
    skipfigs = 1; % 1 = skip making a figure if it exists already    
    save_figs = 1; % 1 = make and save figures (VERY time consuming, takes about 90% of klustest's time)
    fast_figures = 1; % 1 = 4-5x faster figure saving, but at a slightly lower quality
    hill_spacing = 950; % mm, spacing between peaks
    cutoff_dist = 120; % mm, area around bin locations to combine
    % mapset.binsize      = 32; % (mm) firing rate map bin size

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
    mapset = pdata.mapset;
    pos = pdata.pos;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Run through clusters
    disp(sprintf('Analysing clusters...'))
    ucis = unique(cluma.uci); % list of unique cells in sdata    
    
    for pp = 1:nparts % for every part
        part_now = part_names{pp};
        disp(sprintf('\t%s',part_now))      
            
        % session_times = pdata.session_times;
        % pidx = pos.pot > session_times(pp,1) && pos.pot < session_times(pp,2); % index for position data in this part

        loopout = looper(length(ucis));
        for uu = 1:length(ucis)
            uci = ucis{uu};
            idx = find( ismember(cluma.uci,uci) & cluma.partn==pp );

%% >>>>>>>>>> Get (planar) firing rate map etc
            % rmap = cluma.ratemap_planar{idx};
            amap = cluma.planar_amap{idx,1};
            % dmap = pdata.(part_now).dwellmap_planar;

            spacing_bins = hill_spacing/mapset.binsize;
            field_index = ceil(sort(unique([size(amap,2)/2:spacing_bins:size(amap,2)-spacing_bins/2,size(amap,2)/2:-spacing_bins:spacing_bins/2])));
            field_index(3) = []; % remove central field
            non_field_index = ceil(sort(unique([size(amap,2)/2+spacing_bins/2:spacing_bins:size(amap,2)-spacing_bins/2,size(amap,2)/2-spacing_bins/2:-spacing_bins:spacing_bins/2])));
            
            amap_idx = zeros(1,size(amap,2));
            amap_idx(field_index) = 1;
            dist_idx = bwdist(amap_idx,'euclidean');
            dist_log_fields = dist_idx < cutoff_dist/mapset.binsize;
            
            amap_idx = zeros(1,size(amap,2));
            amap_idx(non_field_index) = 1;
            dist_idx = bwdist(amap_idx,'euclidean');
            dist_log_nfields = dist_idx < cutoff_dist/mapset.binsize; 
            
            % get the repeat score
            % this one is specific to the horizontal axis
            central_row = amap( ceil(size(amap,1)/2),: );
            repeat_score1 = mean(central_row(dist_log_fields),'omitnan') - mean(central_row(dist_log_nfields),'omitnan');
            
            % this one averages vertically and will include offset fields
            fmat = repmat(dist_log_fields,size(amap,1),1);
            fnmat = repmat(dist_log_nfields,size(amap,1),1);            
            repeat_score2 = mean(amap(fmat),'omitnan') - mean(amap(fnmat),'omitnan');
            
            if ~any(ismember(cluma.Properties.VariableNames,'repetition_score')) % if the column(s) do not exist yet
                cluma.repetition_score = NaN(size(cluma,1),2); % preallocate
            end            
            cluma.repetition_score(idx,:) = [repeat_score1 repeat_score2];
            if ~any(ismember(cluma.Properties.VariableNames,'amap_cent')) % if the column(s) do not exist yet
                cluma.amap_cent = cell(size(cluma,1),1); % preallocate
            end            
            cluma.amap_cent(idx) = { single(central_row) };    
            
            loopout = looper(loopout);            
        end
        
        disp(sprintf('\t...saving figure'))              
        fname = [data_dir '\' rnow '\' dnow '\3Danalysis\part_figures\' part_now '_repetition.png'];
        [~,~,~] = mkdir([data_dir '\' rnow '\' dnow '\3Danalysis\part_figures']);        
        if (exist(fname,'file') && skipfigs) || ~save_figs  % if the figure exists and we don't want to overwrite it
            % do nothing
        else
            klustfig_repetition(pdata,cluma,rnow,dnow,pp,fast_figures,'off',fname,part_now);
        end    
        disp(sprintf('\t...done'))                      
    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    cluma = rmprop(cluma,'pdata');
    cluma = addprop(cluma,{'pdata'},{'table'});
    cluma.Properties.CustomProperties.pdata = pdata;
    cname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    save(cname,'cluma'); % save session data
    analysis_log({'hill_repeat_analysis_v2'},1,'version',{'v2.0.0'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('%s has finished. It took %0.3f seconds or %0.3f minutes',stk.name,toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');























