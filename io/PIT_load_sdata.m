function PIT_load_sdata
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% PIT_load_sdata  load raw data into an sdata table
% Run through raw data files stored in individual rat directories, load them all
% and add them to a single summary table
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
% version 1.0.0, Release 08/11/22 Code conception
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    if any(strcmp(evalin('base','who'),'sdata')) && ~overwrite_sdata && maintain_sdata % if we want to use an mtint currently held in memory
        mtint = evalin('base','sdata');        
        m = whos('sdata');     
        disp(sprintf('\t...using sdata held in memory (%.1fMb)',m.bytes./1e+6));
    elseif ~exist(config.sname,'file') || config.overwrite_sdata
        disp(sprintf('\t...building sdata'));    
        
        %% code to create sdata
        rmstring = {'Unused','unused','.','..','media','testing'}; % directories we want to ignore consistently
        dirs1 = {config.place_raw_data_dir};

        %% for every major directory
        sdata_big = table;
        pdata_big = table;
        for aa = 1:length(top_dirs) % for every major directory
            cd(dirs1{aa}); % move to this folder    
            disp(sprintf('%s directory...',dirs1{aa})) % display this

            % get a list of all folders (should be rat folders)
            rnames = dir;
            rnames = {rnames.name}';
            rnames(~ismember(rnames,rmstring)); % remove stuff we want to ignore, like '.' or 'Unused' etc    

            %% for every rat directory
            for rr = 1:length(rnames) % for every rat directory
                cd([dir1 '\' rnames{rr}]); % move to its folder
                disp(sprintf('\t%s directory...',[dir1 '\' rnames{rr}])) % display this
                dnames = dir; % list the folders in this rats directory
                dnames = {dnames.name}'; % get their names
                dnames(~ismember(dnames,rmstring)); % remove stuff we want to ignore, like '.' or 'Unused' etc    

                %% for every date directory
                for dd = 1:length(dnames) % for every date directory
                    cd([dir1 '\' rnames{rr} '\' dnames{dd}]); % move to it
                    disp(sprintf('\t\t...%s',dnow)) % display this
                    
                    if exist('000_ready_for_audit.txt','file')

                        out_name = 'kwiktint';
                        load([pwd '\klustest\' out_name '\' out_name '_sdata.mat'],'sdata'); % load this session's sdata
                        load([pwd '\klustest\' out_name '\' out_name '_pdata.mat'],'pdata'); % load this session's pdata          

                        sdata_big = [sdata_big; sdata]; % concatenate tables
                        pdata_big = [pdata_big; pdata]; % concatenate tables

                        keyboard
                    end
                end
            end
        end
        
        
        
        
        
        
        
        
        
        
        

        info = whos('sdata');
        siz = info.bytes / 1000000;
        disp(sprintf('\t...saving sdata (%.1fMb)',siz)) % I have tried my best to keep the mtint >50Mb for most recordings
        save(config.sname,'sdata','-v7.3');
    else
        m = dir(config.sname);
        disp(sprintf('\t...loading sdata (%.1fMb)',m.bytes./1e+6));
        load(config.sname,'sdata');
    end














































