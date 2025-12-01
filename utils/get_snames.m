function [snames,cname,nsess,fnames] = get_snames(varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% get_snames  get session names based on available .set files
% This function just gathers some session names by asking the user to select .set files
% and outputs some useful forms of this info
%
% USAGE:
%       [snames,cname,nsess,fnames] = get_snames() process with default settings
% 
%       [snames,cname,nsess,fnames] = get_snames(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'data_dir'      -   (default = current directory) String, directory top search for .set files
% 
%       'mode'          -   (default = 1) Scalar, 1 = find all .set files, otherwise dialog box will ask user to select files
%
%       'subs'          -   (default = 0) Scalar, 1 = include subdirectories, otherwise only top directory is included
%
% OUTPUT:
%       snames   - output as a vector
%       cname   - output as a vector
%       nsess   - output as a vector
%       fnames   - output as a vector
%
% EXAMPLES:
%       % run function in current directory, find all .set files
%       [snames,cname,nsess,fnames] = get_snames;
%
% See also: gittest uipickfiles fileparts

% HISTORY:
% version 1.0.0, Release 24/11/16 created so this functionality can be shared between 2D and 3D functions
% version 2.0.0, Release 13/10/21 renamed get_snames for non-planar grid cell project
% version 2.0.1, Release 13/10/21 added variable input arguments
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
    p = inputParser;
    addParameter(p,'data_dir',pwd,@(x) isstring(x) | ischar(x));       
    addParameter(p,'mode',1,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'subs',0,@(x) isnumeric(x) && isscalar(x));    
    parse(p,varargin{:});
    config = p.Results;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    if config.mode==1
        if config.subs==1
            snames = dir([config.data_dir '\**/*.set']);            
        else
            snames = dir([config.data_dir '\*.set']);
        end
        
        if ~numel(snames)
            error('ERROR: no detected .set files, check working directory matches data directory... exiting');
        end
        snames = {snames(:).name};
        snames = snames(:);
    else
        snames = uipickfiles('FilterSpec','*.set','Output','cell','Prompt','Select the sessions to analyse...');
    end
    nsess = numel(snames);

    % sort out parameters for later and reduce .set file names
    fnames = cell(length(snames),1);
    if length(snames) == 1
        nnow = snames{1};
        [~,nme,~] = fileparts(nnow);
        cname = nme;  
        snames{1} = nme;    
        fnames{1} = [config.data_dir '\' nme];
    else
        for ff = 1:length(snames)
            nnow = snames{ff};
            [~,nme,~] = fileparts(nnow);
            snames{ff} = nme;
            fnames{ff} = [config.data_dir '\' nme];
            if ff == 1 % the first filename
                cname = nme;
            elseif ff == length(snames) % the last filename
                cname = [cname '_' nme];
            else % middle filenames
                cname = [cname '_' nme '_'];
            end
        end
    end
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function just gathers some session names by asking the user to select .set files
%   it outputs some useful forms of this info
%   [snames,cname,nsess] = getSNAMES
%
%%%%%%%% Outputs
%   snames = cell array of selected file names, filename with no extension or path
%   cname = a new filename generated using all the input snames concatenated
%   nsess = the number of sessions/set files selected
%   fnames = a cell array of full file names
%
%%%%%%%% Comments 
%   24/11/16 created so this functionality can be shared between 2D and 3D functions
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial settings
if ~exist('get_all','var') || isempty(get_all)
    get_all = 0;
end % if ~exist('get_all','var') || ismepty(get_all)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find sessions




