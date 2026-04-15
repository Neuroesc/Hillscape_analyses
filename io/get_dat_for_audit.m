function sdata_acc = get_dat_for_audit(data_dir,to_skip,varargin)
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
    p = inputParser;
    addRequired(p,'data_dir',@(x) isstring(x) || ischar(x)); 
    addParameter(p,'to_skip',cell(1,2),@(x) iscell(x));        
    parse(p,data_dir,varargin{:});
    config = p.Results;
    to_skip = config.to_skip;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
    disp(sprintf('\tdata directory: %s',data_dir)) % display this
    dnames = dir(data_dir);    
    fnames = {dnames.name}';
    fnames = fnames([dnames.isdir]');
    rat_names = fnames( ~contains(fnames,{'Unused','unused','.','..'}) ); % the rat directories
    disp(sprintf('\t\t...%d rats found',length(rat_names))) % display this
    
    sdata_acc = table;
    for rr = 1:length(rat_names) % for every rat directory
        rnow = rat_names{rr}; % get its name
        disp(sprintf('\trat %s...',rnow)) % display this
        
        dnames = dir([data_dir '\' rnow]);
        fnames = {dnames.name}';        
        fnames = fnames([dnames.isdir]');
        dat_names = fnames( ~contains(fnames,{'Unused','unused','.','..'}) ); % the rat directories
        disp(sprintf('\t\t...%d recordings found',length(dat_names))) % display this
    
        for dd = 1:length(dat_names) % for every recording directory
            dnow = dat_names{dd};
            disp(sprintf('\t\t\t%s',dnow)) % display this
            
            if any( strcmp(to_skip(:,1),rnow) && strcmp(to_skip(:,2),dnow) )
                disp(sprintf('\b: included in skip list')) % display this                
            end            

            alog_name = [data_dir '\' rnow '\' dnow '\_analysis_log.txt'];
            if ~exist(alog_name,'file')
                disp(sprintf('\b: unanalysed, skipping')) % display this                
            else
                alog = readtable(alog_name,'FileType','text','ReadVariableNames',1,'VariableNamingRule','preserve','Delimiter','\t');
                if ~any( strcmp(alog.analysis,'klustest') )
                    if ~any( strcmp(alog.analysis,'3Dreconstruction') )
                        disp(sprintf('\b: (WARNING: 3Dreconstruction not found)')) % display this  
                    end
                    disp(sprintf('\b: klustest not run yet, skipping')) % display this  
                else
                    load([data_dir '\' rnow '\' dnow '\klustest\sdata.mat')
                    sdata_acc = [sdata_acc; sdata];
                end
            end
        end
    end


































