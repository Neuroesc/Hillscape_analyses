function clumaa = run_fun_for_audit(config)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% get_dat_for_audit  gather all data for an expriment
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
% See also: klustest analysis_log

% HISTORY:
% version 1.0.0, Release 16/02/23 Initial release
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
    % p = inputParser;
    % addRequired(p,'data_dir',@(x) isstring(x) || ischar(x)); 
    % addRequired(p,'data_out_dir',@(x) isstring(x) || ischar(x));     
    % addRequired(p,'to_skip',@(x) iscell(x)); 
    % addParameter(p,'rats_to_skip',cell(1,1),@(x) iscell(x));     
    % addParameter(p,'f',cell(1,1),@(x) iscell(x)); 
    % parse(p,data_dir,to_skip,varargin{:});
    % config = p.Results;
    start_dir = pwd;

    if ~isfield(config,'to_skip') || isempty(config.to_skip)
        config.to_skip = cell(1,2);
    end
    if ~isfield(config,'rats_to_skip') || isempty(config.rats_to_skip)
        config.rats_to_skip = cell(1,2);
    end    
    warning('off','MATLAB:table:RowsAddedExistingVars');

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
    disp(sprintf('\tdata directory: %s',config.data_dir)) % display this
    dnames = dir(config.data_dir);    
    fnames = {dnames.name}';
    fnames = fnames([dnames.isdir]');
    rat_names = fnames( ~contains(fnames,{'Unused','unused','.','..'}) ); % the rat directories   
    % rat_names = {'RG6'}
    disp(sprintf('\t\t...%d rats found',length(rat_names))) % display this
    
    clumaa = table;
    caa = 0;
    posdata = table;
    p_idx = 1;
    for rr = 1:length(rat_names) % for every rat directory
        rnow = rat_names{rr}; % get its name
        disp(sprintf('\trat %s...',rnow)) % display this
        
        dnames = dir([config.data_dir '\' rnow]);
        fnames = {dnames.name}';        
        fnames = fnames([dnames.isdir]');
        dat_names = fnames( ~contains(fnames,{'Unused','unused','.','..'}) ); % the rat directories
        disp(sprintf('\t\t...%d recordings found',length(dat_names))) % display this
    
        for dd = 1:length(dat_names) % for every recording directory
            dnow = dat_names{dd};
            disp(sprintf('\t\t\t%s',dnow)) % display this
            
            if any( strcmp(config.to_skip(:,1),rnow) & strcmp(config.to_skip(:,2),dnow) )
                disp(sprintf('\t\t\tincluded in skip list... skipping')) % display this   
                continue
            elseif any( strcmp(config.rats_to_skip,rnow) )
                disp(sprintf('\t\t\trat included in skip list... skipping')) % display this   
                continue                
            else
                cd([config.data_dir '\' rnow '\' dnow]); % move to the directory
            end         

%% >>>>>>>>>> 3D trajectory reconstruction
            if 0
                % GIT_reconstruction('resync',1,'repos',1,'retri',1,'remerge',1);
                CLUMA_3d_reconstruction;
            end

%% >>>>>>>>>> load sdata, useful for minor functions that can't do this on their own
            if 0
                cname = [config.data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
                load(cname,'cluma')
                pdata = cluma.Properties.CustomProperties.pdata;

                % add UCI & partn
                cluma.uci = cell(size(cluma,1),1);
                cluma.partn = NaN(size(cluma,1),1);
                for ii = 1:size(cluma,1)
                    cluma.uci{ii} = [rnow '_' dnow '_' cluma.electrode{ii}];

                    if contains(cluma.session_name{ii},'a_arena')
                        cluma.partn(ii) = 1;
                    elseif contains(cluma.session_name{ii},'b_hills')
                        cluma.partn(ii) = 2;
                    elseif contains(cluma.session_name{ii},'c_arena')
                        cluma.partn(ii) = 3;
                    end
                end
                save(cname,'cluma'); % save session data                
            end
            
%% >>>>>>>>>> add 3D trajectory to pdata
            if 0
                cname = [config.data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
                load(cname,'cluma')
                pdata = cluma.Properties.CustomProperties.pdata;

                CLUMA_3d_reconstruction;
                dat = CLUMA_load_merged_file([config.data_dir '\' rnow '\' dnow '\cluma.positions'],'positions');
                pdata.old_pos = pdata.pos;
                pdata.pos = dat.pos;
                dat = CLUMA_load_merged_file([config.data_dir '\' rnow '\' dnow '\cluma.positions3'],'positions3');
                pdata.pos.poz = dat.pos.poz;

                cluma = rmprop(cluma,'pdata');
                cluma = addprop(cluma,{'pdata'},{'table'});
                cluma.Properties.CustomProperties.pdata = pdata;
                cname = [config.data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
                save(cname,'cluma'); % save session data                
            end
            if 0 % short snippet just to add new 3D HD data to pdata
                cname = [config.data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
                load(cname,'cluma')
                pdata = cluma.Properties.CustomProperties.pdata;

                CLUMA_3d_reconstruction;
                dat = CLUMA_load_merged_file([config.data_dir '\' rnow '\' dnow '\cluma.positions3'],'positions3');
                pdata.pos.rz = dat.pos.rz;
                pdata.pos.gz = dat.pos.gz;
                pdata.pos.bz = dat.pos.bz;

                cluma = rmprop(cluma,'pdata');
                cluma = addprop(cluma,{'pdata'},{'table'});
                cluma.Properties.CustomProperties.pdata = pdata;
                cname = [config.data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
                save(cname,'cluma'); % save session data                
            end

%% >>>>>>>>>> fit frame to maze data
            overwrite_maze_frame = 0;
            if 0
                [~,~,~] = mkdir([config.data_dir '\' rnow '\' dnow '\3Danalysis']);
                mname = [config.data_dir '\' rnow '\' dnow '\3Danalysis\maze_boundaries.mat']; 
                fname = [config.data_dir '\' rnow '\' dnow '\3Danalysis\maze_boundaries.png'];                 
                if ~exist(mname,'file') || overwrite_maze_frame
                    [lx,ly,lz] = fit_hill_frame_v2(pdata,0,mname,fname); % fit the maze frame to the data for all 3 parts simultaneously   
                else
                    load(mname,'lx','ly','lz','-mat');             
                end
                pdata.maze_frame = single([lx(:) ly(:) lz(:)]);   

                cluma = rmprop(cluma,'pdata');
                cluma = addprop(cluma,{'pdata'},{'table'});
                cluma.Properties.CustomProperties.pdata = pdata;
                cname = [config.data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
                save(cname,'cluma'); % save session data                 
            end

%% >>>>>>>>>> hill analysis
            overwrite_hill_analysis = 0;
            if 0
                alog_name = [config.data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed, skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if any( strcmp(alog.analysis,{'hill_analysis'}) ) && ~overwrite_hill_analysis
                        disp(sprintf('\t\t\talready analysed with hill_analysis... skipping')) % display this  
                    else
                        hill_analysis_v2(config.data_dir,rnow,dnow);
                    end
                end                
            end

%% >>>>>>>>>> stability analysis
            overwrite_stab_analysis = 0;
            if 0
                alog_name = [config.data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed, skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if any( strcmp(alog.analysis,{'stability_analysis'}) ) && ~overwrite_stab_analysis
                        disp(sprintf('\t\t\talready analysed with hill_analysis... skipping')) % display this  
                    else
                        PIT_stability_analysis(config.data_dir,rnow,dnow);
                    end
                end                
            end

%% >>>>>>>>>> hill repetition function
            overwrite_rep_analysis = 0;
            if 0
                alog_name = [data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed, skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if any( strcmp(alog.analysis,{'hill_repeat_analysis'}) ) && ~overwrite_rep_analysis
                        disp(sprintf('\t\t\talready analysed with hill_repeat_analysis... skipping')) % display this  
                    else
                        hill_repeat_analysis_v2(data_dir,rnow,dnow);
                    end
                end                
            end    
          
%% >>>>>>>>>> manually categorise cells
            overwrite_cell_typing = 0;
            if 0
                alog_name = [config.data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed, skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if any( strcmp(alog.analysis,{'cell_type_v2'}) ) && ~overwrite_cell_typing
                        disp(sprintf('\t\t\talready analysed with get_manual_cell_type... skipping')) % display this  
                    elseif ~any( strcmp(alog.analysis,{'hill_analysis'}) )
                        disp(sprintf('\t\t\thill_analysis function run yet... skipping')) % display this  
                    else
                        get_manual_cell_type_v2(data_dir,rnow,dnow);
                    end
                end                
            end 

%% >>>>>>>>>> Behaviour analysis
            overwrite_behave = 0;
            if 0
                alog_name = [config.data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed, skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if any( strcmp(alog.analysis,{'PIT_behaviour'}) ) && ~overwrite_behave
                        disp(sprintf('\t\t\talready analysed with PIT_behaviour... skipping')) % display this  
                    else
                        PIT_behaviour(config.data_dir,rnow,dnow);
                    end
                end                
            end 

%% >>>>>>>>>> Place field analysis
            overwrite_fields = 0;
            if 0
                alog_name = [config.data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed, skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if any( strcmp(alog.analysis,{'hill_field_analysis'}) ) && ~overwrite_fields
                        disp(sprintf('\t\t\talready analysed with hill_field_analysis... skipping')) % display this  
                    else
                        hill_field_analysis(config.data_dir,rnow,dnow);
                    end
                end                
            end 

%% >>>>>>>>>> plot all place cells (do nothing else)
            if 0
                alog_name = [config.data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed, skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if ~any( strcmp(alog.analysis,{'hill_analysis'}) )
                        disp(sprintf('\t\t\thill_analysis function run yet... skipping')) % display this  
                    else
                        sname = [config.data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
                        load(sname,'cluma'); % load saved session data                        
                        pdata = sdata.Properties.CustomProperties.pdata;                        
                        nparts = numel(unique(pdata.pos.session));

                        for pp = 1:nparts
                            part_now = pdata.part_config.part_names{pp};
                            % klustfig_all_pcells
                            fname = [config.data_dir '\' rnow '\' dnow '\3Danalysis\part_figures\' part_now '_all_place_cells.png'];
                            if exist(fname,'file') && 0 % if the figure exists and we don't want to overwrite it
                                continue % don't make the figure             
                            else                                
                                klustfig_all_pcells(pdata,cluma,rnow,dnow,pp,1,'off',0,fname)
                            end 

                            % klustfig_repetition
                            fname = [config.data_dir '\' rnow '\' dnow '\3Danalysis\part_figures\' part_now '_repetition.png'];
                            if exist(fname,'file') && 0 % if the figure exists and we don't want to overwrite it
                                continue % don't make the figure             
                            else                                
                                klustfig_repetition(pdata,cluma,rnow,dnow,pp,1,'off',fname,part_now);
                            end 

                            % klustfig_hills
                            for uu = 1:length(ucis)
                                fname = [config.data_dir '\' rnow '\' dnow '\3Danalysis\part_figures\' uci '_projections.png'];
                                if exist(fname,'file') && 0 % if the figure exists and we don't want to overwrite it
                                    continue % don't make the figure             
                                else                                
                                    klustfig_hills(pdata,cluma,uci,fast_figures,'off',fname);
                                end 
                            end                            
                        end                          
                    end
                end                
            end 

%% >>>>>>>>>> pitch modulation analysis
            if 0
                alog_name = [config.data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed, skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if any( strcmp(alog.analysis,{'pitch_tuning'}) ) && 1
                        disp(sprintf('\t\t\talready analysed with pitch_modulation_analysis... skipping')) % display this  
                    else
                        PIT_pitch_tuning(config.data_dir,rnow,dnow);
                    end
                end                
            end              
            
%% >>>>>>>>>> concatenate loaded sdata tables
            if 0
                alog_name = [config.data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
                if ~exist(alog_name,'file')
                    disp(sprintf('\t\t\tunanalysed... skipping')) % display this                
                else
                    alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                    if 0%~all( ismember({'klustest','3Dreconstruction','hill_analysis','cell_type'},alog.analysis) )
                        disp(sprintf('\t\t\tnot all required functions run yet... skipping')) % display this  
                    else
                        load([config.data_dir '\' rnow '\' dnow '\cluma\cluma.mat'],'cluma')
                        cluma.pos_idx = ones(size(cluma,1),1).*p_idx; % preallocate
                        clumaa = vertcat(clumaa, cluma); % concatenate cluma
                        
                        pdata = cluma.Properties.CustomProperties.pdata;

                        % keyboard
                        posdata.pos_idx(p_idx) = p_idx;
                        posdata.session_times(p_idx) = { pdata.session_times };
                        posdata.maze_frame(p_idx) = { pdata.maze_frame };
                        posdata.pos(p_idx) = { single(table2array(pdata.pos)) };
                        posdata.pos(p_idx) = { pdata.pos };   
                        posdata.anisotropy_map(p_idx,:) = { pdata.arena1.anisotropy_map pdata.hills.anisotropy_map pdata.arena2.anisotropy_map};
                        posdata.anisotropy_score(p_idx,:) = [pdata.arena1.anisotropy_score pdata.hills.anisotropy_score pdata.arena2.anisotropy_score];
                        p_idx = p_idx+1;
                        caa = 1;
                    end
                end
            end
            %%
            
        end
    end
    cd(start_dir);

    % save concatenated data
    if caa
        clumaa = rmprop(clumaa,{'pdata'});                
        clumaa = addprop(clumaa,{'posdata'},{'table'});        
        clumaa.Properties.CustomProperties.posdata = posdata;
        cname = [config.data_out_dir '\clumaa.mat'];
        save(cname,'clumaa'); % save session data
    end






























