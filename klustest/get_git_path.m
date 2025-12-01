function [rpox,rpoy,rpoz,gpox,gpoy,gpoz,bpox,bpoy,bpoz,pot] = get_git_path(data_dir,snames,recon_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function just loads all of the 3D trajectory info produced by lattice for a session
%   given a set of filenames, an mtint (to check position data lengths) and some settings (optional)
%   [pox,poy,poz] = getLATTICEpath(snames,mtint,config3)
%
%%%%%%%% Inputs
%   snames = cell array of file names, no extensions
%   mtint = an mtint structure formed by klustest
%   config 3 = (optional) settings structure given by klustest3, otherwise defaults are used
%
%%%%%%%% Outputs
%   pox,poy,poz,pot = position x, y, z, t respectively
%
%%%%%%%% Comments
%   10/11/17 created so klustest can take advantage of reconstructed trajectories
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
samp_rate_hz = 50; % synchronised dacqTrack data should be 50Hz
pos_tb = 1 / samp_rate_hz;    
pixel_ratio = 1000;
leds = 1; % typically we only use 1 LED


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the data
rpox = []; 
rpoy = []; 
rpoz = []; % initialise variables
gpox = []; 
gpoy = []; 
gpoz = []; % initialise variables
bpox = []; 
bpoy = []; 
bpoz = []; % initialise variables
pot = [];

nsess = length(snames);
for pp = 1:nsess % for each session
    snow = snames{pp}; % get its filename
    fcheck = [data_dir '\' recon_name '\' snow '_merged.txt'];
    if exist(fcheck,'file') % check it exists
        cdata = readtable(fcheck,'Delimiter','\t');
    else
        error('Predetermined 3D trajectory not found: %s\n ...exiting',fcheck)
    end
    
    rpox = [rpox(:); double( cdata.rpox(:) )];
    rpoy = [rpoy(:); double( -cdata.rpoy(:) )];
    rpoz = [rpoz(:); double( -cdata.rpoz(:) )];

    gpox = [gpox(:); double( cdata.gpox(:) )];
    gpoy = [gpoy(:); double( -cdata.gpoy(:) )];
    gpoz = [gpoz(:); double( -cdata.gpoz(:) )];        

    bpox = [bpox(:); double( cdata.bpox(:) )];
    bpoy = [bpoy(:); double( -cdata.bpoy(:) )];
    bpoz = [bpoz(:); double( -cdata.bpoz(:) )];  

    pot = [pot(:); double( cdata.pot(:) )]; 
end







