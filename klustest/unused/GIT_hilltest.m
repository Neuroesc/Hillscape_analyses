function GIT_hilltest(tetrodes,clusters,cname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%klustest  cluster analysis function
% This function utilises many minor functions to both analyse and plot cluster data generated 
% in Tint. If requested, it will output a figure for each cluster, the cluster space of each
% tetrode, the cross-correlations of every cluster on a tetrode and a session data structure (sdata.mat)
% It will also generate an mtint file (mtint.mat) containing all the tetrode and cluster info.
% klustest(tetrodes,clusters,cname)
%
% USAGE:
%         klustest(tetrodes,clusters,cname)
%
% INPUT:
%         tetrodes - (default = 1:16) the tetrodes to run on in a vector format (i.e. [1 2 3 4 5 6] would run on tetrodes 1 to 6)
%         clusters - (default = 0) the clusters to run on, set this to 0 to run on all clusters
%         rname - the rat name/number as a string - this will be used in the sdata structure so its important to give this
%         cname - (optional) function will automatically look for klustakwiked files named 'kwiktint', but you can change this here (i.e if you chose a custom output name in kwiktint)
%
% EXAMPLES:
%
% See also: KWIKTINT

% HISTORY:
% version 01.0.0, Release 05/08/16 created from an older version
% version 01.1.0, Release 08/08/16 modified readDACQDATA to output waveforms and added post processing of this info
% version 01.2.0, Release 08/08/16 readDACQDATA working, added postprocessing of mtint
% version 01.3.0, Release 09/08/16 main plots working
% version 01.4.0, Release 10/08/16 options added for different styles of plot
% version 02.0.0, Release 11/08/16 added theta autocorrelation sine estimate (O'Mara 2015)
% version 02.1.0, Release 12/08/16 added refractory period violations (Navratilova and McNaughton (2016)
% version 02.1.1, Release 13/08/16 added saving data, figures
% version 02.2.0, Release 14/08/16 added cell type identification
% version 03.0.0, Release 15/08/16 cluster quality assessment added in clusterQUALITY
% version 03.1.0, Release 16/08/16 added cluster space figures
% version 03.1.1, Release 17/08/16 fixed phase plot and phase map
% version 03.2.0, Release 18/08/16 modified clusterQUALITY to deal better with missing channels
% version 03.3.0, Release 18/08/16 added the option to ignore position data, fixed bug with waveform plot
% version 03.3.1, Release 19/08/16 make it so none of the position figures are made if there is no position data, this should be faster
% version 03.4.0, Release 22/08/16 fixed issues with cluster space plot, concentrate on first feature, added legend
% version 03.5.0, Release 23/08/16 added cluster space subplot to cell figure, for Ele, this uses the first clustering feature, the 1st and 2nd highest amplitude channel
% version 04.0.0, Release 25/06/16 added tetrode input to readalldacqdata, this prevents it trying to open unnecessary files, fixed bug in the detection of dead/lfp channels
% version 04.1.0, Release 25/08/16 fixed a minor bug in histograms when there is only one spike
% version 05.0.0, Release 22/11/16 adapted from KlusterAnalysisTINT
% version 05.1.0, Release 30/03/17 started editing the newer version with more emphasis on sdata structure
% version 06.0.0, Release 30/03/17 adapted to new partitioning system using a part_config cell array, function can now split based on any inputs
% version 07.0.0, Release 31/03/17 started moving figure work into subfunctions to contain them, each uses the sdata structure
% version 07.1.0, Release 01/04/17 finished most of the figure subfunctions
% version 08.0.0, Release 02/04/17 added kluspartfig, which generates a summary figure containing all the parts, and correlations between them
% version 09.0.0, Release 05/04/17 fixed multiple problems in clusterQUALITY, including wrong lratios and reliance on kwiktint files, now using clusterQUALITY_v2
% version 10.0.0, Release 13/04/17 added ability to extract parent filenames from .cut file, meaning the function only needs the kwiktint out name
% version 10.1.0, Release 19/04/17 fixed head direction plot and analyses
% version 10.2.0, Release 20/04/17 changed from using mapDATA to mapDATA_v2, latter allows for maps with same limits and converted position data necessary for overdispersion
% version 10.3.0, Release 20/04/17 added overdispersion calculation (my own)
% version 10.4.0, Release 05/05/17 added comments, added easier naming of part_config fields
% version 10.5.0, Release 10/05/17 fixed bugs with interval times in part_config, added figKEYS to handle interval times, added exceptions to figCLUST and overdispersion to handle overlapping trials
% version 10.6.0, Release 12/05/17 fixed bug in getDACQDATA where keypress times were not accumulated in time, fixed bug in partSESS where incorrect interval vector was being loaded
% version 11.0.0, Release 12/05/17 added ability to send emails when completed
% version 11.1.0, Release 17/05/17 added possibility to calculate shuffled spatial measures
% version 11.2.0, Release 19/05/17 replaced read_key with saveKEY, this automatically saves a text file version of the .inp files as it loads them
% version 12.0.0, Release 19/05/17 replaced getPARTconfig with getPARTconfig_v2, this saves a much nicer text file which is easier to edit
% version 12.1.0, Release 07/06/17 fixed a problem in getDACQDATA where it would try to run on all tetrodes even if a subset is specified
% version 12.2.0, Release 22/06/17 fixed error in time allocation
% version 12.3.0, Release 10/07/17 changed part_config specification so that it is done in a loop, added information required by contest and klustest3
% version 13.0.0, Release 25/07/17 added variable pixel ratio capability (big job!) getDACQDATA now uses getDACQDATAHEADERS instead of the horible key_value functions
% version 13.1.0, Release 01/11/17 fixed bug where LFP and session duration did not match
% version 14.0.0, Release 01/11/17 replaced celltype with getCELLTYPE, replaced mapDATA_v3 with mapDATA4, replaced GridAnalaysis with gridSCORE2, replaced AutoCorr with ndautoCORR
% version 15.0.0, Release 04/08/18 overhauled, code cleaned up, data maintained in cm throughout
% version 15.1.0, Release 10/08/18 recoded overdispersion and speed score analysis
% version 16.0.0, Release 04/04/19 moved part_config guts to later in the function, replaced part_config structure to be a table, now saved as .mat
% version 16.0.1, Release 04/04/19 replaced getTRODES and getCNAMES with a combined version called getTRODES
% version 16.0.2, Release 05/04/19 moved pixel ratio handling to postprocessDACQDATA
% version 16.0.3, Release 07/04/19 part indexing is now built logically instead of in a loop (should be faster)
% version 16.0.4, Release 08/04/19 updated spatial mapping analyses
% version 16.1.0, Release 09/04/19 moved HD bining to mapHD, theta curve fitting to getTHETAintrinsic, refractory analysis to getRPVcount, speed analysis to getSPEEDmod      
% version 16.1.1, Release 09/04/19 changed getCELLTYPE for getCELLTYPE2 which outputs a binary approach to cell typing
% version 16.2.0, Release 10/04/19 moved cluster quality to getDACQDATA
% version 16.2.1, Release 10/04/19 better handling of waveforms
% version 16.2.1, Release 13/04/19 figure functions finished, figPART and figCLUST
% version 16.3.0, Release 15/04/19 Finished converting sdata to table format, variables have to be preallocated by subfunction addToTable
% version 16.3.1, Release 15/04/19 progress is now displayed using fileExchange function 'ProgressBar'
% version 16.3.2, Release 15/04/19 improved handling of keypresses in figKEY, added manual definition of start/end keys, part_config now contains files required for each part
% version 17.0.0, Release 25/06/19 klustest3 version 2 created from this code to bring it up to date
% version 18.0.0, Release 18/10/21 renamed GIT_klustest3 for non-planar grid cell project, tailored for this project
% version 18.0.0, Release 18/10/21 getDACQDATA replaced with get_dacq_data, postprocessDACQDATA with postprocess_dacq_data, getLATTICEpath with get_git_path, makeLATTICE_v2 with fit_maze_frame

%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
    % I can extract the rat name from the file path because I know my folder structure, if this is not the case for you, comment out the next line
    pname = pwd; segs = strsplit(pname,'\'); rname = segs{end-1};
%     rname = '999';

    % deal with the other variables
    inps = {'tetrodes','clusters','cname','rname'};
    vals = {'[]','0','''kwiktint''','''000'''};
    for ff = 1:length(inps)
        if ~exist(inps{ff},'var')
            eval([inps{ff} '=' vals{ff} ';']);
        end
    end
    recon_name = '3Dreconstruction';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################################################################################## %% INITIAL SETTINGS / INPUTS
    % Specify partition configuration settings
    part_names              = {'arena1','hills','arena2'}; % the names you want the outputs to be saved as, these cannot start with a number: i.e. 'session1'
    % method of partition, corresponding to each of the names above: 1 = combine everything, 2 = take a recording session, 3 = use digital inputs
    part_methods            = [2 2 2]; % i.e. if part_methods=[2 2 2 1], we will have 4 parts, the first 3 correspond to some combination of recording sessions (there can be multiple ones) and the last one will include all data
    % cell array of vectors indicating which intervals to include in each partition: if method = 1 this does nothing, if method = 2 this should specify which recording sessions to include, if method = 3 this should specify which digital input pairs to include (inf = use all)
    part_intervals          = {1 2 3}; % i.e. if part_methods=[2 2 2], then if part_intervals={1 2 3}, rec 1 will go in part 1, rec 2 in part 2 and rec 3 in part 3    OR     if part_methods=[2 2 2 2 2], then if part_intervals={1 [2 4 5] 3}, rec 1 will go in part 1, rec 2,4 and 5 in part 2 and rec 3 in part 3
    pixel_ratio_override    = []; % set to empty [] to ignore, otherwise this must be n recording sessions long
    config.interval_keys    = {}; % what keypresses were used to delineate trial starts and ends {start,end} these can be integers or characters i.e. {1,2} or {'s','e'}, {'1','2'} is the same as {1,2}
    part_dimensions         = [2 3 2];

    % overrides - general settings for overriding normal klustest functionality
    pconfig_override        = 0; % set to 1 if you want to ignore and overwrite an existing part_config, set to 2 to run with current part_config settings without overwriting anything
    maintain_mtint          = 0; % DEBUGGING ONLY set to 1 to save/load mtint in the base workspace, this saves time when running the function mutliple times (for instance in debugging) but should otherwise be set to 0
    mtint_override          = 0; % set this to 1 to force the production of a new mtint file (do this if you reclustered in tint for instance)
    config.trial_override   = 0; % set this to 1 to force the production of a new .key file (do this if you are not happy with the current trials for instance)
    fig_key_off             = 1; % set this to 1 to suppress display of trials/digital keypresses in an interactive UI
    angle_override          = 0; % set this to 1 to recalculate 3D head pose
    frame_override          = 0; % set this to 1 to override any existing maze outlines/frames
    step_save               = 2; % set to 1 and the function will save sdata etc after every tetrode, if the function crashes on the next run it will skip all the completed tetrodes (saving adds a little time, but might save time overall), set to 2 to override an existing stepsave, set to 0 to do no step saving
    config.cname            = cname;
    
    % map settings - settings used when generating 2D spatial firing rate maps    
    config.rmethod          = 'histogram'; % (default 'histogram') the mapping approach to use
    config.map_padd         = 2; % (default 2) the number of bins to pad spatial maps with
    config.bin_size         = 2; % (default 2) bin size in cm for calculating the rate map (make sure data is in cm or pm/sm values are given)
    config.map_sigma        = 4; % (default 1.5) used by nearest and KDE method, sigma of gaussian to use when smoothing traditional dwell and spike maps, or used as bandwidth of KDE
    config.smethod          = 1; % (default 1) used by nearest method, when set to 1 the dwell and spike maps are smoothed before dividing, 2 means the smoothing is applied after dividing, 3 means no smoothing
    config.min_dwell        = 0.01; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    config.g_sigma          = 5; % (default 10) only used by gaussian method - sigma of gaussian used to weight position point and spike point distances, bigger means smoother maps
    config.min_dist         = 5; % (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    config.max_dist         = 20; % (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    config.srate            = 50; % (default 50) sampling rate of data in Hz, used to calculate time

    % map settings - settings used when generating 3D spatial firing rate maps    
%     config3.rmethod         = 'gaussian'; % (default 'nearest') the mapping approach to use, either 'nearest','gaussian','adaptive','KDE'
%     config3.map_padd        = 2; % (default 2) the number of bins to pad spatial maps with
%     config3.bin_size        = 50; % (default 2) bin size in cm for calculating the rate map (make sure data is in cm or pm/sm values are given)
%     config3.map_sigma       = 2; % (default 1.5) used by nearest method, sigma of gaussian to use when smoothing traditional dwell and spike maps [p s] to smooth independently or [x] to smooth ratemap after division, leave empty [] for no smoothing
%     config3.min_dwell       = 1; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
%     config3.g_sigma         = config3.bin_size*1.5; % (default 10) only used by gaussian method - sigma of gaussian used to weight position point and spike point distances, bigger means smoother maps
%     config3.bandwidth       = config3.bin_size*1.5; % bandwidth for KDE mapping method    
%     config3.min_dist        = 100; % (default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
%     config3.max_dist        = config3.bin_size*5; % maximum distance (mm) from voxel to include data
%     config3.srate           = 50; % (default 50) sampling rate of data in Hz, used to calculate time

    % lattice settings - size and shape    
%     config3.lattice_type    = 1;    
%     config3.lat_size        = [970 970 970]; % lattice size in mm
%     config3.of_size         = [1200 1200 200]; % open field size in mm    
    
    % field settings - settings to use when detecting 2D place fields
    config.frcut            = 0.5; % relative minimum firing rate (% of ratemap max) to be considered a field
    config.arcut            = 9; % minimum number of pixels to be considered a field, typically people use 9 contiguous pixels
    config.minfr            = 1; % (Hz) absolute minimum cutoff firing rate to be considered a field

    % field settings - settings to use when detecting 2D place fields
    config3.thresh_method    = 'zscore'; % method for thresholding ratemap, 'proportion','zscore' or 'cutoff'   
    config3.frcut            = 1; % relative minimum firing rate (% of ratemap max) to be considered a field or Z value is method is 'zscore'
    config3.arcut            = 40; % minimum number of pixels to be considered a field, typically people use 9 contiguous pixels
    config3.minfr            = 0; % (Hz) absolute minimum cutoff firing rate to be considered a field
    
    % spike plot settings
    config.over_smooth      = 13; % number of position data points over which to smooth instantaneous firing rate when calculating overdispersion
    config.time_bins        = 2; % (default 2s) time window over which to compute the spike vs time plot

    % HD settings - settings to use when generating head direction maps
    config.hd_type          = 'histogram'; % (default 'histogram') enter 'density' for a circular kolmogorov smirnov density estimate plot, enter 'histogram' for the traditional histogram polar plot
    config.hd_bins          = 64; % (default 64) the number of bins to use when computing HD plot
    config.hd_sigma         = 0.04; % (default 2) the standard deviation of the gaussian kernel used in HD circular density estimate
    config.hd_boxcar        = 3; % (defualt 3) the number of bins over which to compute the HD histogram boxcar average

    % figure settings
    run_figPART            = 1; % [ figPART ] set to 1 for figures showing the activity of a cluster in each part
    run_figCLUS            = 1; % [ figCLUS ] set to 1 for figures showing the activity of a cluster in all parts together

    % optimization settings
    config.wave_asis        = 1; % set to 1 for waveform analysis and plots
    config.temp_asis        = 1; % set to 1 for spike time analysis and plots
    config.fild_asis        = 1; % set to 1 for place field analysis and plots
    config.grid_asis        = 1; % set to 1 for grid cell analysis and plots
    config.head_asis        = 0; % set to 1 for HD cell analysis and plots
    config.ovrd_asis        = 1; % set to 1 for overdispersion analysis and plots
    config.sisi_asis        = 1; % set to 1 for inter-spike interval analysis and plots
    config.autc_asis        = 1; % set to 1 for spike autocorrelation analysis and plots
    config.spph_asis        = 1; % set to 1 for spike theta phase analysis and plots
    config.sped_asis        = 1; % set to 1 for speed modulation analysis and plots
    config.cell_asis        = 0; % set to 1 for cell type analysis and plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Start analysis and get part_config
    % get the function's name (people may rename it from klustest)
    stk = dbstack;
    function_name = stk.name;

    % starting messages
    tic;
    disp('----------------------------------------------------------------------------');
    time_now = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',function_name,time_now))
    if pconfig_override;                    disp(sprintf('\tWARNING: will override part_config...'));                           end
    if maintain_mtint;                      disp(sprintf('\tWARNING: will use/maintain mtint in base workspace...'));           end
    if mtint_override;                      disp(sprintf('\tWARNING: will override mtint...'));                                 end
    if frame_override;                      disp(sprintf('\t...will ask for new maze frames'));                                 end    
    if step_save==1 && ~mtint_override;     disp(sprintf('\t...will save sdata after every cluster'));                          end
    if step_save==2 || mtint_override;      disp(sprintf('\tWARNING: will overwrite any existing sdata'));                      end
    
    disp(sprintf('\t...%s spatial maps with %.1fcm bins (%.1f padding)',config.rmethod,config.bin_size,config.map_padd));
    disp(sprintf('\t...%.1f%% cutoff for fields, %.1f pixels minimum, %.1fHz minimum',config.frcut.*100,config.arcut,config.minfr));
    disp(sprintf('\t...%s HD maps with %d bins',config.hd_type,config.hd_bins));

    % The part_config table contains basic information about how/where/if to separate the data into different
    % partitions or parts. For instance if we record an open field, then some sort of maze, then the open field
    % again, we may want to divide the data into those 3 parts for separate analysis, or we might want to just 
    % lump everything together if we care about some intrinsic property of the cells
    % We want to save this table so that in the future we can just run the function with the same settings
    [~,~,~] = mkdir([pwd '\klustest\' cname]); % create the folder which will hold a lot of klustest data
    part_num = length(part_names);
    Rat = repmat({rname},part_num,1);
    Part = part_names';
    Method = part_methods';
    Intervals = part_intervals';
    Outname = repmat({cname},part_num,1);
    Dimensions = part_dimensions';
    part_config = table(Rat,Part,Method,Intervals,Outname,Dimensions);

    % load or save part_config
    % as mentioned above, the part_config file contains the basic information required to process a dataset, it really
    % is just a permanent copy of the settings outlined manually at the top of klustest. For more information, like
    % the start and end times of each part, or the files used for each one, see the part_config saved in pdata instead
    % partIO here will try to save the part_config in a file named pconfig_name, if this file already exists it will
    % backup the contents using the current date/time and then save the part_config alongside it non destructively
    % Otherwise it will make a new file and save it there. FOr this reason, if you load pconfig_name you will actually find 
    % a structure named 'part_data', the field named 'part_config' is the last part_config setting used, the other fields are
    % older part_configs, where the numbers in the field name reflect the date/time they were overwritten
    pconfig_name = ['klustest\' cname '\' cname '_part_config.mat'];    
    part_config = partIO(pconfig_name,part_config,pconfig_override);
    disp(sprintf('\t...rat name: %s',rname))      

    % start to prepare the pdata (part or session data) and sdata (cell data) tables
    % To save space and make things easier for the user I have divided the data into two files
    % the pdata structure contains data which applies to the whole session, like positions or dwell maps
    % the sdata table contains data for each cluster/part (one per row). The idea being that tables are easier to concatenate
    % to build a full data set
    pdata = struct; % will hold part data
    pdata.part_names = part_names;
    pdata.combined_name = cname;
    pdata.rat = rname;
    pdata.analysed = time_now;
    pdata.directory = pwd;
    pdata.interval_keys = config.interval_keys;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Get session names and tetrodes
    % When sessions are cluster cut together (which you should do with multiple sessions recording the same cells)
    % tint saves the filenames of the individual sessions in the .cut file it makes after cluster cutting with
    % klustakwik. I found extracting this data to be the easiest way to get this information, and it ensures the
    % function will only ever be used to analyse sessions cluster cut in one file
    % As it happens, we can also check to see which tetrodes have a cut file and from this work out which tetrodes
    % we can and can't analyse
    disp(sprintf('Assessing data...'))
    [tets,snames,cutname] = getTRODES(cname); 
    
    disp(sprintf('\t...read %s',cutname));
    nsess = numel(snames);
    disp(sprintf('\t...working on %d sessions: %s',nsess,strjoin(snames,', ')));    
    if ~isempty(tetrodes) % if specific tetrodes were requested
        if ~isempty(setdiff(tetrodes,tets)) % if there were some tetrodes requested that don't have files
            disp(sprintf('\t...WARNING: skipping tetrodes %s, cut file not found',mat2str(setdiff(tetrodes,tets))))            
        end
        tetrodes = intersect(tetrodes,tets);
    else % if the user is running in automatic detection of tetrodes
        tetrodes = tets;
    end
    disp(sprintf('\t...tetrodes: %s accounted for',mat2str((tetrodes(:))')))

    % accumulate
    pdata.tetrodes = tetrodes; % list of tetrodes analysed
    pdata.session_names = snames; % names of the recording sessions (i.e. files) used, cell array
    pdata.sessions = nsess; % the number of recording sessions used
    pdata.part_config = part_config; % a copy of the part_config, this will be extended to include more data that the saved part_config though, table format
    pdata.config = config; % a copy of the configuration struct, so ratemaps etc can be remade with identical settings if need be
    pdata.config3 = config3; % a copy of the configuration3D struct, so ratemaps etc can be remade with identical settings if need be

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Load or generate mtint
    % The first step of the data analysis is just to load all of the Axona file data into a matlab friendly format
    % this is largely unprocessed, raw data. This is stored in a structure called mtint.
    % We create this and will use it a little bit, but mostly it just serves as a foundation for klustest. Our 
    % data and stats will be stored in the pdata and sdata structures. Klustest shouldn't really change the contents
    % of this structure as it serves as a non-destructive reference
    % Because the mtint can be quite large and bulky, we have some time saving options - if one has already been created 
    % we can just load it instead of making a new one. If you re-cluster cut the data you should override any existing
    % mtint though (or delete it manually) otherwise the clusters will not match those in Tint. If we are debugging the function
    % we can just opt to use an mtint stored in the base workspace (so no loading necessary)
    disp(sprintf('Fetching DACQ data (mtint)...'));
    mname = ['klustest\' config.cname '\' config.cname '_mtint.mat']; % the filename of the mtint (one that exists or one we will make)
    
    if any(strcmp(evalin('base','who'),'mtint')) && ~mtint_override && maintain_mtint % if we want to use an mtint currently held in memory
        mtint = evalin('base','mtint');        
        m = whos('mtint');     
        disp(sprintf('\t...using mtint held in memory (%.1fMb)',m.bytes./1e+6));
    elseif ~exist(mname,'file') || mtint_override
        disp(sprintf('\t...building mtint'));    
        mtint = get_dacq_data(config,snames,tetrodes); % my replacement for readalldacqdata
        disp(sprintf('\t...post-processing mtint'));

        % deal with manual override of pixel ratio if necessary
        if numel(pixel_ratio_override) ~= numel(mtint.pos.header)
            disp(sprintf('\tWARNING: number of pixel ratios in part_config (%d) does not equal number of recordings (%d)... ignoring manual values',numel(pixel_ratio_override),length({mtint.pos.header.pixels_per_metre})));
        else
            for pr = 1:numel([mtint.pos.header.pixels_per_metre])
                mtint.pos.header(pr).pixels_per_metre = pixel_ratio_override(pr);
            end
        end    

        % post process position data etc
        % this function adds running speed, head direction data and also smoothes, interpolates the data
        % removes jumps in the data which are too quick to be natural, converts the position data to cm
        % and a bunch of other stuff
        mtint = postprocess_dacq_data(mtint); % my replacement for postprocess_DACQ_data

        % save mtint
        info = whos('mtint');
        siz = info.bytes / 1000000;
        disp(sprintf('\t...saving mtint (%.1fMb)',siz)) % I have tried my best to keep the mtint >50Mb for most recordings
        save(mname,'mtint','-v7.3');
    else
        m = dir(mname);
        disp(sprintf('\t...loading mtint (%.1fMb)',m.bytes./1e+6));
        load(mname,'mtint');
    end

    % basic session info
    duration = mtint.pos.total_duration;
    disp(sprintf('\t...total session time: %ds',duration))
    disp(sprintf('\t...cut file is made up of %d recordings',numel(mtint.header)))

    if maintain_mtint
        assignin('base','mtint',mtint); % leave mtint in base workspace
    end 
    disp(sprintf('\t...done'));

    % accumulate
    pdata.date = mtint.header_date; % date of recording
    pdata.duration = duration; % total duration of the recording (s)
    pdata.recording_times = mtint.pos.trial_duration; % the length in s of each recording
    pdata.session_details = table; % table to hold extra session info in one place
    pdata.session_details.session_names = pdata.session_names; % cell array of session names
    pdata.session_details.session_lengths = mtint.pos.trial_duration(:); % length of each recording (s)
    pdata.session_details.session_ends = cumsum(mtint.pos.trial_duration(:)); % the time (in concatenated klustest time) where this session ends in the data (s)
    pdata.session_details.session_starts = pdata.session_details.session_ends(:) - pdata.session_details.session_lengths(:); % the time (in concatenated klustest time) where this session starts in the data (s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Add interval/part time delineations to part_config (now stored in pdata)
    % We want to divide a long recording session into partitions or parts, these might be whole
    % recording sessions or they might be based on digital keypresses or you might just want the
    % whole thing in one long session. In any case partSESS2 takes the keypresses and the desired
    % method defining each part and gives the start + end time of each so we can separate the data later
    % If you are using keypresses (i.e. intervals) partSESS2 will also run figKEYS which will plot every
    % trial and let you correct keypresses if necessary
    disp(sprintf('Preparing partitions...'))
    pdata.config.trial_override = config.trial_override;
    pdata.config.fig_key_off = fig_key_off;    
    [pdata,nparts] = partSESS2(mtint,pdata); % process the part_config to get the start and end time of each partition

    % basic session info
    disp(sprintf('\t...%d parts',nparts))
    disp(sprintf('\t...part methods: %s',mat2str(pdata.part_config.Method(:))))
    disp(sprintf('\t...part names: %s',strjoin(pdata.part_config.Part(:),', ')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Prepare important data
    disp(sprintf('Extracting initial data...'))

    % get the position data for the whole session
    % for klustest3 this data is collected from multiple cameras tracked using dacqTRack
    % and is alread reconstructed in 3D 
    % get_git_path also flips the Y and Z axes to correct the data
    disp(sprintf('\t...3D positions'))    
    [rpox,rpoy,rpoz,gpox,gpoy,gpoz,bpox,bpoy,bpoz,~] = get_git_path(pwd,snames,recon_name,mtint.pos);
    pot = double(mtint.pos.ts); % extract the time stamp of each position value
    pov = double(mtint.pos.speed(:)); % the running speed throughout the session
    ppm = ones(size(pot)).*1000; % pixels per metre

    % extract the animal's 3D head pose from the tracking LEDs
    disp(sprintf('\t...3D orientation'))        
    pos = [rpox rpoy rpoz gpox gpoy gpoz bpox bpoy bpoz];
    oname = [pwd '\' recon_name '\3Dorientation.txt'];
    %angle_override = 1;
    if ~exist(oname,'file') || angle_override || 0
        rname_temp = 'generic_hpc';
        head_angles = get_3d_orientation(pos,rname_temp); % odat = [X, Y, Z, yaw, pitch, roll]
        writetable(head_angles,oname,'FileType','text','WriteVariableNames',true);
    else
        head_angles = readtable(oname,'FileType','text','ReadVariableNames',true);
    end
    
    % get the position data sampling rate (should be 50hz) or 0.05s
    srate = mtint.pos.header(1).sample_rate_num(1,1); % sampling rate
    sinterval = 1 / srate; % sampling interval

    % display results
    disp(sprintf('\t...positions read: %d',numel(rpox(:))));
    disp(sprintf('\t...tracking colours: %d',size(mtint.pos.led_pos,3)));
    disp(sprintf('\t...median pixel ratio: %dppm',nanmedian(ppm)));
    disp(sprintf('\t...sample rate: %dHz (%.2fs)',srate,sinterval));

    % get the LFP data
    % as explained below in the LFP section, klustest loads the first available LFP data saved in mtint, this corresponds to
    % LFP channel 1 in dacqUSB (i.e. the primary LFP channel displayed during recording)
    % this can be changed here if necessary, but if we try to load a non-existant field of mtint the function will crash
    disp(sprintf('\t...LFP'))        
    Fs = mtint.lfp(1).Fs(1,1);
    lfp = double(mtint.lfp(1).lfp(:,1));
    lfpt = (0:length(lfp)-1)'/Fs; % make a vector for time    

    % filter LFP to get the theta band    
    [b,a] = butter(4,[6 12]/(Fs/2)); % Generate 4th order butterworth filter coefficients for theta band [6 12] Hz
    lftheta = filtfilt(b,a,lfp); % Apply filter to data using zero-phase filtering  
    
    % accumulate data
%     pox = mean([rpox gpox bpox],2,'omitnan'); % rat X center of mass
%     poy = mean([rpoy gpoy bpoy],2,'omitnan'); % rat Y center of mass
%     poz = mean([rpoz gpoz bpoz],2,'omitnan'); % rat Z center of mass
    pox = head_angles.meanRx(:); % X center of the rat's head 
    poy = head_angles.meanRy(:); % Y center of the rat's head 
    poz = head_angles.meanRz(:); % Z center of the rat's head 
    poh = table2array(head_angles(:,4:6)); % yaw, pitch, roll
    
    pdata.pox = single([pox rpox gpox bpox]);
    pdata.poy = single([poy rpoy gpoy bpoy]);
    pdata.poz = single([poz rpoz gpoz bpoz]);    
    pdata.pot = pot; 
    pdata.pov = pov;
    pdata.ppm = ppm; 
    pdata.sampling_rate = srate;
    pdata.sampling_interval = sinterval;

    % display results
    disp(sprintf('\t...LFP samples read: %d',numel(lfp)));
    disp(sprintf('\t...LFP sample rate: %dHz',Fs));    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################################################################################## %% Get maze outlines  
    disp(sprintf('\t...maze boundaries'))    
    for pp = 1:nparts % for every partition
        part_now = part_names{pp}; % the name of the current part
        part_times = pdata.(part_now).times; % the time pairs (intervals) corresponding to this part
        pindax = zeros(size(pox));       
        for ii = 1:length(part_times(:,1)) % for every pair of time values (interval) associated with this part  
            pindax = logical(pindax + (pot > part_times(ii,1) & pot < part_times(ii,2))); % logical array, same length as pox, 1 means position sample falls into this part
        end 
        ppox = pox(pindax); % pos x
        ppoy = poy(pindax); % pos y
        ppoz = poz(pindax); % pos z        

        % fit a frame to this partition data
        lname = ['klustest\' cname '\' part_now '_frame.mat'];
        if exist(lname,'file') && ~frame_override
            load(lname);
        else          
            [lx,ly,lz] = fit_maze_frame(ppox,ppoy,ppoz);
            save(lname,'lx','ly','lz') 
        end 
        pdata.(part_now).lx = lx;
        pdata.(part_now).ly = ly;
        pdata.(part_now).lz = lz;
    end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################################################################################## %% Run through tetrodes and clusters
    disp(sprintf('Analysing tetrode/cluster data...'))
    sdata = table; % will hold all the cluster data
    fdata = table; % will hold all the cluster place field data
    
    if step_save==1 && ~mtint_override
        if exist([pwd '\klustest\' cname '\' cname '_fdata.mat'],'file')
            load([pwd '\klustest\' cname '\' cname '_fdata.mat'],'fdata'); % load field data
            load([pwd '\klustest\' cname '\' cname '_sdata.mat'],'sdata'); % load cluster data
            load([pwd '\klustest\' cname '\' cname '_pdata.mat'],'pdata'); % load session data  
        end
    end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% For every available tetrode
    % The next loop focuses on analysing a single tetrode, there is not a lot of info we really need about tetrodes
    % because the meat of the data is really per cluster or part
    % This is especially true now that cluster quality has been moved to gatDACQDATA, so that it doesn't have to be
    % computed on every run of klustest. However, we will load the waveforms for this tetrode before continuing to the clusters
    for ee = 1:length(tetrodes) 
        tet = tetrodes(ee); % tet = the current tetrode
        disp(sprintf('\tLoading tetrode %d...',tet));
        
        % get a vector of clusters we want to analyse
        if isempty(clusters) || clusters==0
            clus = unique(mtint.tetrode(tet).cut); % find the clusters logged for this tetrode
        else
            clus = clusters;
        end

        % check to see if there are any clusters
        clus_count = numel(clus);   
        disp(sprintf('\t\t...%d spikes detected',mtint.tetrode(tet).nspike_cut));
        disp(sprintf('\t\t\t...%d data clusters detected',sum(clus~=0))); 
        if ~clus_count || ~any(clus) % if there are no clusters, or if there is only a noise cluster
            continue % skip analysis
        end
        disp(sprintf('\t\t\t\t...starting analysis'));

        % get the channel waveforms for this tetrode    
        if config.wave_asis
            disp(sprintf('\t\t\t\t...getting waveforms'));                   
            nspikes = mtint.tetrode(tet).nspike_cut;
            waves = zeros(nspikes,32,4,'single');
            rnow = 1;            
            for ssn = 1:length(pdata.session_names)
                fnamen = pdata.session_names{ssn};
                [~,c1,c2,c3,c4,~] = get_SPIKESV([fnamen,'.',num2str(tet)]);                               
                waves(rnow:rnow+length(c1(:,1))-1,:,1) = c1;
                waves(rnow:rnow+length(c2(:,1))-1,:,2) = c2;
                waves(rnow:rnow+length(c3(:,1))-1,:,3) = c3;
                waves(rnow:rnow+length(c4(:,1))-1,:,4) = c4;                
                rnow = rnow+length(c4(:,1));
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% For every detected cluster    
        % The next loop focuses on analysing a single cluster from this tetrode, almost immediately we will skip it if
        % the cluster is noise (i.e. cluster 0), in some cases it might be interesting to have that information, but as we
        % dump all our noise spikes in cluster 0 this may be extraordinarily large, so it makes more sense to skip it for speed

        % get all the spike times for this tetrode
        spiketime = double(mtint.tetrode(tet).ts);      
        disp(sprintf('\t\t\t\t...analysing clusters'));                   
    
        for cc = 1:length(clus) % for every cluster   
            sdatac = table; % to hold this cluster's data
            fdatac = table; % to hold this cluster's place field data
            tdata = table;
            clu = clus(cc); % clu = the current cluster
            uci = ['r' rname '_' pdata.date '_t' num2str(tet) '_c' num2str(clu)]; % unique cell identifier - a string with [rat name, session date, electrode, cluster];

            % if this is the noise cluster don't continue any further
            if ~clu      
                continue
            end
            disp(sprintf('\t\t\t\t\tcluster %d ',clu));                   
            
            % if we are saving in cluster steps, check if this cluster's data have already been processed/saved
            if step_save==1 && ~isempty(sdata) && ~mtint_override
                if any(sdata.tet(:)==tet & sdata.clu(:)==clu)
                    continue
                end
            end

            % clu_indx is a logical vector, length = N of spikes, true if spike is in cluster clu, false if not
            % mtint.tetrode.cut contains a vector specifying which spikes are in which cluster, each row corresponds to a tetrode
            % so to get the data for one tetrode we need to do: mtint.tetrode(tetrode # we want).cut
            clu_identity = mtint.tetrode(tet).cut; % clu_assign is a vector of numbers, one for each spike, each number corresponds to a cluster                        
            clu_indx = clu_identity==clu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% For each part  
            % The next loop focuses on dividing the data into its different parts and analysing the cell within each of these
            % The very first section deals with extracting the data relevant to this part. I've tried to keep this process
            % logical so as to increase speed, the downside is that if there is a large amount of data this approach may fail
            % due to a lack of memory. However, as the matrices are logical (i.e. 1 byte per value) I imagine this would require
            % an enormous recording session with many intervals
            part_names = pdata.part_config.Part;
            for pp = 1:nparts % for every partition  
                disp(sprintf('\b| part %d ',pp));                   
     
                sdatap = table; % to hold part data for this cluster 
                tdatap = table;
                part_now = part_names{pp}; % the name of the current part
                part_times = pdata.(part_now).times; % the time pairs (intervals) corresponding to this part
            
                % find what data falls into the intervals associated with this part
                % start with positions
                if ~isfield(pdata,part_now)
                    pdata.(part_now) = pdata.(part_now);
                end
                pindax = logical(sum(pot' >= part_times(:,1) & pot' <= part_times(:,2),1));
                ppox = pox(pindax); % pos x
                ppoy = poy(pindax); % pos y
                ppoz = poz(pindax); % pos z                    
                ppot = pot(pindax); % pos time
                ppov = pov(pindax); % pos running speed
                ppoh = poh(pindax,:); % pos HD   
                part_duration = sum(pindax(:))*pdata.sampling_interval; 

                % accumulate
                pdata.(part_now).pox = single(ppox);
                pdata.(part_now).poy = single(ppoy);
                pdata.(part_now).poz = single(ppoz);                    
                pdata.(part_now).pot = single(ppot); 
                pdata.(part_now).poh = single(ppoh); 
                pdata.(part_now).pov = single(ppov);
                pdata.(part_now).duration = part_duration;

                % spike data next
                sindax = logical(sum(spiketime' > part_times(:,1) & spiketime' < part_times(:,2) & clu_indx',1)); % spikes
                nindax = logical(sum(spiketime' > part_times(:,1) & spiketime' < part_times(:,2) & ~clu_identity',1)); % noise 
                pspt = spiketime(sindax);
                sidx = knnsearch(ppot,pspt);
                pspx = ppox(sidx); % spike x
                pspy = ppoy(sidx); % spike y
                pspz = ppoz(sidx); % spike y                
                psph = ppoh(sidx,:); % direction
                pspv = ppov(sidx); % running speed                

                % accumulate, these variables don't need to be preallocated as they should always be filled
                sdatap.rat = {pdata.rat};
                sdatap.date = str2double(pdata.date);
                sdatap.directory = {pwd};                
                sdatap.partn = pp;
                sdatap.part = {part_now};                
                sdatap.uci = {uci};
                sdatap.tet = tet;
                sdatap.clu = clu;
                sdatap.spike_index = {uint32(sidx)}; % this index can be used to get the spike values from the position data values i.e. pdata.(part_now).pox(sdata.spike_index)
                sdatap.spike_times = {double(pspt)}; % the actual spike times
                sdatap.isod = mtint.clu_quality(tet).isolation_distances(clu);
                sdatap.lratio = mtint.clu_quality(tet).lratios(clu);
                sdatap.spikes = numel(pspx);
                sdatap.duration = part_duration;
                sdatap.frate = numel(pspx) / part_duration;
                minspikes = 1; % the minimum number of spikes a cluster has to have before we do the following analyses
                pdata.minspikes = minspikes;
                
                % accumulate for tdata
                tdatap.rat = {pdata.rat};
                tdatap.date = str2double(pdata.date);
                tdatap.directory = {pwd};                
                tdatap.partn = pp;
                tdatap.part = {part_now};                
                tdatap.uci = {uci};
                tdatap.tet = tet;
                tdatap.clu = clu;    
                tdatap.pos = {single(ppox) single(ppoy) single(ppoz) single(ppoh) single(ppov) single(ppot)};
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Waveforms
                % all waveforms for this electrode are loaded above by getSPIKESV, but we want only the waveforms
                % corresponding to this cluster and part, which is what this section is for
                % The waveforms are sadly not saved anywhere because they would take up a lot of space (minimum 100mb)
                % so instead we have to just load them from the spike files every time we run klustest
                % Instead we save the mean +/- SD for each channel as a representation of the cluster
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'waveform_mean',cell(1,4),'waveform_stdv',cell(1,4),'waveform_rms',NaN(1,4),'waveform_snr',NaN(1,4),'waveform_max',NaN(1,4),'waveform_min',NaN(1,4),'waveform_maxt',NaN(1,4),'waveform_mint',NaN(1,4),'waveform_width',NaN(1,4),'waveform_params',NaN(1,2),'waveform_snrs',NaN(1,4));

                % get the waveforms for this cluster   
                if config.wave_asis && numel(pspx)>minspikes
                    waves2 = waves(sindax,:,:);
                    for w = 1:4 % for every recording channel
                        wav = double(waves2(:,:,w));
                        ch = nanmean(wav,1);
                        chs = nanstd(wav,[],1);
                        chrms = nanmean(rms(wav,1));
                        nsrms = nanmean(rms(waves(nindax,:,w),1));
                        [maxval,maxt] = max(ch);
                        [postminval,postmint] = min(ch(maxt:end));
                        postmint = postmint + maxt - 1;
                        width = postmint - maxt;
                        width = width * (1000/50); 

                        sdatap.waveform_mean(1,w) = {ch};
                        sdatap.waveform_stdv(1,w) = {chs};
                        sdatap.waveform_rms(1,w) = chrms;  
                        sdatap.waveform_snr(1,w) = chrms/nsrms; % https://en.wikipedia.org/wiki/Signal-to-noise_ratio                       
                        sdatap.waveform_max(1,w) = maxval;
                        sdatap.waveform_min(1,w) = postminval;   
                        sdatap.waveform_maxt(1,w) = maxt;
                        sdatap.waveform_mint(1,w) = postmint;                         
                        sdatap.waveform_width(1,w) = width;
                        
                        % signal to noise Liu et al. (2014) Quality Metrics of Spike Sorting Using Neighborhood Components Analysis
                        % https://dx.doi.org/10.2174%2F1874120701408010060
                        snrs = 1/size(wav,1) .* nansum((nanmax(wav,[],2)-nanmin(wav,[],2)) ./ (2.*nanstd(wav-nanmean(wav,1),[],2)));
                        sdatap.waveform_snrs(1,w) = snrs;
                    end

                    % get the cell's width of waveform using the waveform with the highest mean amplitude
                    [amp,idx] = max(sdatap.waveform_max(1,:));
                    wow = sdatap.waveform_width(1,idx);
                    sdatap.waveform_params(1,:) = [amp wow]; % accumulate data
                    clear waves2 sindax pindax % clear all the waveform data, we don't need it again and it is quite large
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Ratemap, dwellmap, spatial analyses, grid analyses, overdispersion
                % place field, gridness and overdispersion analyses all require a firing rate map
                % so they are combined here in the same loop where the firing rate map is generated
                % firing rate maps are generated by mapDATA which can utilise a number of different methods
                % it will also accept a previously computed dwellmap to speed up computation
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'ratemap_planar',cell(1,1),'ratemap_surficial',cell(1,1),'spatial_info_bsec',NaN(1,2),'spatial_info_bspike',NaN(1,2),'mutual_info',NaN(1,2),'sparsity',NaN(1,2),'spatial_coherence',NaN(1,2));
                sdatap = addToTable(sdatap,'place_fields',NaN);
                fdatap = table;
                tdatap = addToTable(tdatap,'pos_flat_curve',cell(1,1),'pos_planar',cell(1,1),'pos_surficial',cell(1,1),'dwellmap_planar',cell(1,1),'dwellmap_surficial',cell(1,1));
                
                if part_config.Dimensions(pp)==2 % if this is an arena session
                    pos_planar = [ppox,ppoy,ones(size(ppox))];                    
                    pdata.(part_now).pos_planar = single(pos_planar); % surficial and planar are the same in this maze
                    
                else % if this is a hills session
                    if ~isfield(pdata.(part_now),'pos_flat_curve')
                        % fit curve to x,z data (side view to capture hill shape)
                        [x,y] = prepareCurveData( ppox, -ppoz );
                        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                        opts.Display = 'Off';
                        opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
                        opts.StartPoint = [-24.9070167618834 -729.888034703216 421.720294911524 -24.9070167618834 -308.167739791691 421.720294911524 -24.9070167618834 113.552555119833 421.720294911524 -24.9070167618834 535.272850031357 421.720294911524 -24.9070167618834 956.993144942882 421.720294911524 -24.9070167618834 1378.71343985441 421.720294911524];
                        [fitresult, ~] = fit( x,y, fittype( 'gauss6' ), opts ); % Fit model to data.
                        ppoz_curved = -fitresult(ppox); % flattened Z positions

                        % surface projection (onto curved surface, probably only useful for plotting)
                        pos_flat_curve = [ppox,ppoy,ppoz_curved];
                        pdata.(part_now).pos_flat_curve = single(pos_flat_curve);
                        tdatap.pos_flat_curve = {single(pos_flat_curve)};
                    end

                    if ~isfield(pdata.(part_now),'pos_planar')
                        % planar projection (down onto flat surface)
                        pos_planar = [ppox,ppoy,ones(size(ppox))];
                        pdata.(part_now).pos_planar = single(pos_planar);
                        tdatap.pos_planar = {single(pos_planar)};                        
                    end

                    if ~isfield(pdata.(part_now),'pos_surficial')
                        % surficial projection (onto curve and then flattened)
                        lx = pdata.(part_now).lx; % get bounds of maze
                        newx = min(lx(:))-20 : 1 : max(lx(:))+20; % x points along edge of maze (1mm resolution)
                        newy = fitresult(newx); % approximate y coordinates we expect on the maze along edge of maze
                        points = [newx(:) newy(:)];
                        ds = cumsum( [0; sqrt(sum((points(1:end-1,:)-points(2:end,:)).^2,2))] ); % calculate multidimensional distance from start of curve
                        ppox_flat = interp1(newx,ds,ppox,'linear','extrap'); % find new x coordinates based on distance along curve
                        ppoz_flat = ones(size(ppox)); % new z coordinates will just be bottom of maze
                        pos_surficial = [ppox_flat,ppoy,ppoz_flat];
                        pdata.(part_now).pos_surficial = single(pos_surficial);
                        tdatap.pos_surficial = {single(pos_surficial)};                        
                        
                        lx_flat = interp1(newx,ds,lx,'linear','extrap'); % find new x coordinates based on distance along curve 
                        pdata.(part_now).lx_flat = lx_flat(:);
                    end                   
                end

                % planar (vertical projection) ratemap     
                % make this for every session
                pos = pdata.(part_now).pos_planar(:,1:2) ./ 10;
                spk = pos(sidx,:);
                m = [min(pdata.(part_now).lx(:),[],'omitnan') min(pdata.(part_now).ly(:),[],'omitnan') max(pdata.(part_now).lx(:),[],'omitnan') max(pdata.(part_now).ly(:),[],'omitnan')] ./ 10;
                if isfield(pdata.(part_now),'planar_speedlift') % if a previously computed dwellmap exists, use this to save time
                    [rmap1,dmap1,~,c1] = graphDATA(pos,spk,'method',config.rmethod,'kern','epanechnikov','binsize',config.bin_size,'ssigma',config.map_sigma,'maplims',m,'padding',config.map_padd,'speedlift',pdata.(part_now).planar_speedlift);  
                else
                    [rmap1,dmap1,~,c1] = graphDATA(pos,spk,'method',config.rmethod,'kern','epanechnikov','binsize',config.bin_size,'ssigma',config.map_sigma,'maplims',m,'padding',config.map_padd);
                    pdata.(part_now).planar_speedlift = c1.speedlift;
                end
                sdatap = addToTable(sdatap,'ratemap_planar',{rmap1});
                pdata.(part_now).dwellmap_planar = single(dmap1);
                tdatap.dwellmap_planar = {single(dmap1)};
                
                % surficial (surface projection) ratemap    
                % make this for hilly sessions only
                if part_config.Dimensions(pp)==3
                    pos = pdata.(part_now).pos_surficial(:,1:2) ./ 10;
                    spk = pos(sidx,:);
                    m = [min(pdata.(part_now).lx_flat(:),[],'omitnan') min(pdata.(part_now).ly(:),[],'omitnan') max(pdata.(part_now).lx_flat(:),[],'omitnan') max(pdata.(part_now).ly(:),[],'omitnan')] ./ 10;
                    d = [pdata.(part_now).lx_flat(:) pdata.(part_now).ly(:)]./10;
                    if isfield(pdata.(part_now),'surficial_speedlift') % if a previously computed dwellmap exists, use this to save time
                        [rmap2,dmap2,~,c2] = graphDATA(pos,spk,'method',config.rmethod,'kern','epanechnikov','binsize',config.bin_size,'ssigma',config.map_sigma,'maplims',m,'speedlift',pdata.(part_now).surficial_speedlift);  
                    else
                        [rmap2,dmap2,~,c2] = graphDATA(pos,spk,'method',config.rmethod,'kern','epanechnikov','binsize',config.bin_size,'ssigma',config.map_sigma,'maplims',m);
                        pdata.(part_now).surficial_speedlift = c2.speedlift;
                    end
                    sdatap = addToTable(sdatap,'ratemap_surficial',{rmap2});
                    pdata.(part_now).dwellmap_surficial = single(dmap2);    
                    tdatap.dwellmap_surficial = {single(dmap2)};                    
                end
     
%% ########## %% ratemap analysis
                % This function incorporates a number of different analyses and you can check the function for detailed descriptions of these
                % I will point out that the spatial information content value provided by this function is calculated as in 
                % Skaggs et al. (1996) Theta Phase Precession in Hippocampal Neuronal Populations and the Compression of Temporal Sequences
                % https://doi.org/10.1002/(SICI)1098-1063(1996)6:2%3C149::AID-HIPO6%3E3.0.CO;2-K
                % I note this because the method differs slightly between papers and I have found this one reflects the data the best
                spatm1 = spatialMETRICS(rmap1,dmap1);

                if part_config.Dimensions(pp)==2 % if this is an arena session                    
                    % accumulate data                    
                    sdatap = addToTable(sdatap,'spatial_info_bsec',[spatm1.spatial_information NaN],'spatial_info_bspike',[spatm1.spatial_information_perspike NaN],'mutual_info',[spatm1.mutual_info NaN],'sparsity',[spatm1.sparsity NaN],'spatial_coherence',[spatm1.spatial_coherence NaN]);
                
                else % if this is a hills session
                    spatm2 = spatialMETRICS(rmap2,dmap2);

                    % accumulate data                    
                    sdatap = addToTable(sdatap,'spatial_info_bsec',[spatm1.spatial_information spatm2.spatial_information],'spatial_info_bspike',[spatm1.spatial_information_perspike spatm2.spatial_information_perspike],'mutual_info',[spatm1.mutual_info spatm2.mutual_info],'sparsity',[spatm1.sparsity spatm2.sparsity],'spatial_coherence',[spatm1.spatial_coherence spatm2.spatial_coherence]);                    
                end
                
%% ########## %% place field analysis
                if config.fild_asis
                    % this analysis is pretty standard and is based on what I used before, a good citation would be:
                    % Park, Dvorak and Fenton (2011) Ensemble Place Codes in Hippocampus: CA1, CA3, and Dentate Gyrus Place Cells Have Multiple Place Fields in Large Environments
                    % "A place field was defined as any contiguous set of 9 or more pixels with greater that 0 AP/s firing rate that shared at least one side with another pixel in the field."
                    % https://dx.doi.org/10.1371%2Fjournal.pone.0022349
%                     [pdata,sdatap,fdatap,nfields] = getPFIELDS3(pdata,sdatap,uci,part_now);                    

                    % accumulate data
%                     sdatap = addToTable(sdatap,'place_fields',nfields);
                end

%% ########## %% grid characteristic analysis
                sdatap = addToTable(sdatap,'amap_planar',cell(1,1),'amap_surficial',cell(1,1));
                sdatap = addToTable(sdatap,'amap_planar_mask',cell(1,1),'amap_surficial_mask',cell(1,1));                
                sdatap = addToTable(sdatap,'grid_score',NaN(1,2),'grid_wavelength',NaN(1,2),'grid_radius',NaN(1,2),'grid_orientation',NaN(1,2),'grid_anisotropy',NaN(1,2));
                sdatap = addToTable(sdatap,'field orientation',NaN(1,2));
                        
                if config.grid_asis          
                    % create planar autocorrelation
                    % I'm not entirely sure who to attribute this to, but autocorrelations are pretty ubiquitous now
                    % this function was adapted from xPearson (adapted for speed and to make it work on n-dimensions)
                    amap1 = ndautoCORR(rmap1,rmap1,50);
                    if ~all(isnan(amap1(:)))
                        % autocorrelation analysis
                        % There are a number of different grid score approaches contained here so see the function itself for a description
                        % my preferred method is the default, which is from Langston et al. (2010) Development of the Spatial Representation System in the Rat
                        % In this method rings are cut from the autocorrelation at different distances from the centre and the standard rotation and correlation
                        % is performed on each one. The highest of these grid scores is used, the radius of the ring is the grid spacing. A sine wave
                        % is fitted to the values in the ring and this is used to estimate the positions of the fields. The grid orientation is the
                        % angle from the centre to the first one of these fields, counter-clockwise. The grid field radius is taken as the radius of the centre field.
                        [~,gdata1] = get_grid_score(amap1,config.bin_size,'method','allen');

                        % accumulate data
                        sdatap = addToTable(sdatap,'amap_planar',{amap1});
                        sdatap = addToTable(sdatap,'amap_planar_mask',{gdata1.peaks_mask});

                        if part_config.Dimensions(pp)==2 % if this is an arena session
                            sdatap = addToTable(sdatap,'grid_score',[gdata1.grid_score NaN],'grid_wavelength',[gdata1.wavelength NaN],'grid_radius',[gdata1.radius NaN],'grid_orientation',[gdata1.grid_orientation NaN],'grid_anisotropy',[gdata1.elongation NaN]);
                            sdatap = addToTable(sdatap,'field orientation',[gdata1.field_orientation NaN]);

                        else % if this is a hills session
                            amap2 = ndautoCORR(rmap2,rmap2,50);
                            [~,gdata2] = get_grid_score(amap2,config.bin_size,'method','allen');

                            % accumulate data
                            sdatap = addToTable(sdatap,'amap_surficial',{amap2});
                            sdatap = addToTable(sdatap,'amap_surficial_mask',{gdata1.peaks_mask});                        
                            sdatap = addToTable(sdatap,'grid_score',[gdata1.grid_score gdata2.grid_score],'grid_wavelength',[gdata1.wavelength gdata2.wavelength],'grid_radius',[gdata1.radius gdata2.radius],'grid_orientation',[gdata1.grid_orientation gdata2.grid_orientation],'grid_anisotropy',[gdata1.elongation gdata2.elongation]);                        
                            sdatap = addToTable(sdatap,'field orientation',[gdata1.field_orientation gdata2.field_orientation]);                    
                        end
                    end
                end      

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% %% ################################################################# %% Overdispersion              
%             % Taken from: Fenton et al. (2010) Attention-like modulation of hippocampus place cell discharge
%             % "In the first method, the entire session was divided into 5-sec intervals. For each interval we calculated the expected number of spikes, exp, as [equation in paper]
%             % where ri is the time-averaged rate at location i, and ti is the time spent in location i during the pass. Only intervals during which exp ? 5.0 AP were 
%             % used to calculate overdispersion since the overall firing rate of place cells is ~1.0 AP/sec.
%             % For each selected 5-sec interval, we then calculated z, the normalized standard deviation of obs, the observed number of spikes as [equation in paper]
%             % z measures the deviation of observed discharge from expected in standard deviation units. Overdispersion in turn is the variance of 
%             % the z distribution for a set of passes. The outcome of this calculation was found to be indistinguishable from the somewhat different 
%             % method previously used (Fenton and Muller, 1998)."
%             % https://doi.org/10.1523/JNEUROSCI.5576-09.2010
%             [overz,overd,overr] = getOVERDISPERSE([mapdata.poxnew mapdata.poynew],ppot,pspt,ratemap);
%             % essentially we compare the instantaneous firing of the cell at every position to the firing we would expect at that position according to the firing rate map
%             % the overdispersion z is the deviation of the expected from the observed, the overdispersion value is the standard deviation of this distribution
% 
%             % accumulate data
%             sdatap = addToTable(sdatap,'over_dispersion_z',{overz},'over_dispersion',overd,'over_dispersion_r',overr);                      
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% HD analyses  
            % Normal 2D yaw only HD maps
            sdatap = addToTable(sdatap,'hd_spikemap',cell(1,1),'hd_ratemap',cell(1,1),'hd_rayleigh',NaN,'hd_max',NaN,'hd_mean',NaN,'hd_stdev',NaN);
            if isfield(pdata.(part_now),'hd_dwellmap') % if a previously computed dwellmap exists, use this to save time
                [~,~,hd_spikemap,hd_ratemap,r,mx,mn,sd] = mapHD(pdata.(part_now).hd_dwellmap,ppoh(:,1),psph(:,1),config);   
            else
                [~,hd_dwellmap,hd_spikemap,hd_ratemap,r,mx,mn,sd] = mapHD([],ppoh(:,1),psph(:,1),config);                
                pdata.(part_now).hd_dwellmap = hd_dwellmap;
            end                
            sdatap = addToTable(sdatap,'hd_spikemap',{single(hd_spikemap(:))},'hd_ratemap',{single(hd_ratemap(:))},'hd_rayleigh',r,'hd_max',mx(1),'hd_mean',mn,'hd_stdev',sd);
            tdatap.hd_dwellmap = {single(pdata.(part_now).hd_dwellmap)};                
            
            % HD firing rate maps are now computed by map_hd_3d, like the above it will accept a precomputed dwell time map to speed up computation
            sdatap = addToTable(sdatap,'hd_3d_spikemap',cell(1,1),'hd_3d_ratemap',cell(1,1));      
            if isfield(pdata.(part_now),'hd_3d_dwellmap') % if a previously computed dwellmap exists, use this to save time
                [~,~,~,spikemap,~,ratemap] = map_hd_3d(ppoh,psph,'bins',64,'sigma',10,'dwellmap',pdata.(part_now).hd_3d_dwellmap);
            else
                [~,~,~,spikemap,dwellmap,ratemap] = map_hd_3d(ppoh,psph,'bins',64,'sigma',10);
                pdata.(part_now).hd_3d_dwellmap = dwellmap;
            end                
            sdatap = addToTable(sdatap,'hd_3d_spikemap',{spikemap},'hd_3d_ratemap',{ratemap});
            tdatap.hd_3d_dwellmap = {single(pdata.(part_now).hd_3d_dwellmap)};                
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Spike theta phase analyses   
                % By default any LFP analyses in klustest make use of the first LFP channel contained in
                % the mtint (i.e. LFP channel 1 in dacqUSB)
                % This can be changed at the top of klustest if you want, but it will crash if the LFP file does not
                % exist. Ideally the 'best' LFP channel would be used, where 'best' would be defined as the 
                % channel with the strongest theta or something. But for speed, we are just using the 1st one here
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'theta_phase_mean',NaN,'theta_phase_r',NaN,'theta_phase_max',NaN,'theta_phase_dist',cell(1,1));

                if config.spph_asis && numel(pspx)>minspikes
                    % bin the theta phase data
                    ai = reshape(deg2rad(-180:5:540),[],1); % doing this means we have bins symmetrical around zero
                    
                    % we must calculate the spike phase here precisely for every spike time, rather than indexing into the position thet phase
                    % this is for reasons of resolution (LFP and spikes are sampled at 48kHz positions are only sampled at 50Hz)
                    % This process is fairly well established:
                    % "To obtain a theta phase angle for each spike, LFPs were first bandpass filtered (fourth-order Chebyshev, r = 0.5, MATLAB 
                    % filter and filtfilt routines; 6–10 Hz) before a Hilbert transform was applied to obtain the instantaneous phase angle."
                    % from van der Meer and Redish (2011) Theta Phase Precession in Rat Ventral Striatum Links Place and Reward Information (https://doi.org/10.1523/JNEUROSCI.4869-10.2011)             
                    h = hilbert(lftheta);
                    phase = mod(angle(h),2*pi);  
                    pspp = interp1(lfpt(:),phase,'nearest');
                    spp2 = reshape([pspp; (pspp+2*pi)],[],1);
                    yi = histcounts(spp2,ai);

                    % directional analyses on phase data
                    mx2p = ai(yi == max(yi)); % preferred angle (location of max frate)

                    % accumulate data
                    sdatap = addToTable(sdatap,'theta_phase_mean',circ_mean(pspp),'theta_phase_r',circ_r(pspp(:)),'theta_phase_max',mx2p(1),'theta_phase_dist',{yi});
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% ################################################################# %% Inter-spike interval analyses, theta and bursting analyses, refractory period analyses   
                % ISI analyses are quick and there shouldn't really be any reason to need to disable them,
                % however, it is assumed here that if you don't want refractory period or theta autocorrelogram analyses
                % ISI analyses are not required. The intrinsic theta fit can be seen in the klustest figure output
                % as a red line fitted to the theta spike autocorrelogram
                % This can also be used in combination with the frequency of global (LFP) theta as an indication of phase precession
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'isi_dist',cell(1,1),'isi_fdist',cell(1,1),'isi_fwhmx',NaN,'isi_half_width',NaN(1,4));
                sdatap = addToTable(sdatap,'burst_index',NaN,'burst_length_median',NaN,'burst_length_mean',NaN);
                sdatap = addToTable(sdatap,'intrinsic_theta_index',NaN,'intrinsic_theta_frequency',NaN,'intrinsic_theta_fit',NaN,'t500_spike_autocorr',cell(1,1),'t500_spike_autofit',cell(1,1));
                sdatap = addToTable(sdatap,'rpv_total',NaN,'rpv_proportion',NaN,'rpv_false_positive1',NaN,'rpv_false_positive2',NaN,'rpv_censored',NaN,'t25_spike_autocorr',cell(1,1));

                if config.sisi_asis && numel(pspx)>minspikes           
%% ########## %% inter-spike interval
                    % I don't remember exactly where I took this analysis from, but extracting the full width at half maximum is a fairly straightforward
                    % concept so I don't think it requires citation.
                    [~,idata,isis] = getISIhalfwidth(pspt);

                    % accumulate data
                    if isnan(idata.fwhmx) 
                        sdatap = addToTable(sdatap,'isi_dist',{NaN(1,101,'single')},'isi_fdist',{NaN(1,101,'single')},'isi_fwhmx',idata.fwhmx,'isi_half_width',NaN(1,4));
                    else
                        sdatap = addToTable(sdatap,'isi_dist',{single(interp1(idata.adist(:,1),idata.adist(:,2),0:0.5:50))},'isi_fdist',{single(interp1(idata.fdist(:,1),idata.fdist(:,2),0:0.5:50))},'isi_fwhmx',idata.fwhmx,'isi_half_width',[idata.half_max idata.fdist(idata.hwidth_ps,1)']);  
                    end

%% ########## %% burst index
                    % from Mizuseki et al. (2012) Activity Dynamics and Behavioral Correlates of CA3 and CA1 Hippocampal Pyramidal Neurons 
                    % "The spike-burst index was defined as the fraction of spikes with <6 ms ISIs (Harris et al., 2001)"
                    % https://doi.org/10.1002/hipo.22002
                    bindx = zeros(size(pspt));
                    bindx([isis; NaN] < 6 | [NaN; isis] < 6) = 1; % bindx is an index of all spikes sharing an isi less than 6ms
                    bindx = logical(bindx);
                    burst_index = sum(bindx) / numel(pspt);
                    sts = regionprops(bindx,'Area');

                    % accumulate data
                    sdatap = addToTable(sdatap,'burst_index',burst_index,'burst_length_median',nanmedian([sts.Area].'),'burst_length_mean',nanmean([sts.Area].'));

%% ########## %% Intrinsic theta    
                    % from van der Meer and Redish (2011) Theta Phase Precession in Rat Ventral Striatum Links Place and Reward Information
                    % "To quantify the degree and frequency of theta modulation in single cells, we used the method used by Royer et al. (2010). 
                    % First we computed the autocorrelogram of the cell, in 10 ms bins from -500 to +500 ms, normalized it to the maximum value
                    % between 100 and 150 ms (corresponding to theta modulation), and clipped all values above 1. Then we fit the following function [function in paper]
                    % where t is the autocorrelogram time lag, and a-c, w, and t1–2 were fit using the fminsearch optimization function in MATLAB. A measure
                    % of theta modulation strength, the “theta index,” was defined as a/b, which intuitively corresponds to the ratio of the sine fit relative to
                    % the baseline in the autocorrelogram."
                    % https://doi.org/10.1523/JNEUROSCI.4869-10.2011
                    [thi,thf,thr,~,c500,f500] = getTHETAintrinsic(pspt); % we also extract the frequency of theta and the goodness of the theta fit

                    % accumulate
                    sdatap = addToTable(sdatap,'intrinsic_theta_index',thi,'intrinsic_theta_frequency',thf,'intrinsic_theta_fit',thr,'t500_spike_autocorr',{single(c500(:))'},'t500_spike_autofit',{single(f500(:))'});
                    
%% ########## %% Refractory period analyses           
                    % refractory period violation analyses
                    % from Navratilova et al. (2016) Grids from bands, or bands from grids? An examination of the effects of single unit contamination on grid cell firing fields.
                    % although they took these analyses from pre-existing papers
                    % "There are few reliable methods for estimating the contamination of a unit isolated from tetrode recordings (see DISCUSSION). The only method that does not rely 
                    % on the same measures that are used for spike sorting is to check the number of spikes that occur within the refractory period of another spike..."
                    % https://doi.org/10.1152/jn.00699.2015
                    [nrpv,prpv,fp1,fp2,censored,~,c25] = getRPVcount(pspt,isis,part_duration); % we extract the number of spikes in the refractory period, what proportion of the total this is etc                    

                    % accumulate data
                    sdatap = addToTable(sdatap,'rpv_total',nrpv,'rpv_proportion',prpv,'rpv_false_positive1',fp1,'rpv_false_positive2',fp2,'rpv_censored',censored,'t25_spike_autocorr',{single(c25(:))'});
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Speed modulation analysis   
                % This analysis is taken from Kropff et al (mEC speed cell paper)
                % It is unlikely to be of great use to HPC people, although many place cells are speed modulated
                % Related to this, you can also compute the relationship between theta power and running speed
                
                % preallocate - in order to concatenate tables after each loop they have to have the same columns, unfortunately this means all variables have to be preallocated in case they are not filled during the loop, the syntax here should be straightforward though if you want to add more
                sdatap = addToTable(sdatap,'speed_score',NaN,'speed_slope',NaN,'speed_y_intercept',NaN,'speed_frate_curve',cell(1,1));
                
                if config.sped_asis && numel(pspx)>minspikes
                    [svals,stime,sscore,sslope,sintcpt,crve,~] = getSPEEDmod(ppot,ppov,pspt);

                    % accumulate data
                    pdata.(part_now).speed_dwell_time = single(interp1(svals,stime,0:1:50,'linear'));     
                    sdatap = addToTable(sdatap,'speed_score',sscore,'speed_slope',sslope,'speed_y_intercept',sintcpt,'speed_frate_curve',{single(interp1(svals,crve,0:1:50,'linear'))});
                end                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Part figure and finish loop
                % figPART is the main figure, it creates a summary for each cluster in each part, using the info in mtint, pdata and sdatap
                % I have tried to optimise it as far as possible
                if numel(pspx)>minspikes
                    if run_figPART
                        fig_parts_hills(mtint,pdata,sdatap);
                        %putvar(mtint,pdata,sdatap,fdatap);                        
                        %keyboard
                    end
                end

                % concatenate sdatap (this part data) into sdatac (this cluster's data)
                if ~isempty(sdatap)
                    sdatac = [sdatac;sdatap];
                end
                if ~isempty(fdatap)
                    fdatac = [fdatac;fdatap];
                end
                if ~isempty(tdatap)
                    tdata = [tdata;tdatap];
                end                
            end % this ends the parts loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Cluster figures and finish loop            
            if run_figCLUS
                if nparts > 1 % if there is more than one part (otherwise this figure is not useful)
                    fig_clus_hills(mtint,pdata,sdatac); % overall figure with ratemaps etc
                    %putvar(mtint,pdata,sdatac);                        
                    %keyboard
                end 
            end

            % concatenate sdatac (this cluster's data) into sdata (data for all clusters)
%             sdatac2 = sdatac;
%             fdatac2 = fdatac;   
            if ~isempty(sdatac)            
                sdata = [sdata;sdatac];
            end          
            if ~isempty(fdatac)            
                fdata = [fdata;fdatac];
            end

            if step_save>0
                save([pwd '\klustest\' cname '\' cname '_fdata.mat'],'fdata'); % save field data
                save([pwd '\klustest\' cname '\' cname '_sdata.mat'],'sdata'); % save cluster data
                save([pwd '\klustest\' cname '\' cname '_tdata.mat'],'tdata'); % save cluster data                
                save([pwd '\klustest\' cname '\' cname '_pdata.mat'],'pdata'); % save session data                
            end
%             loopout = looper(loopout);
        end % this ends the cluster loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Electrode figures and finish loop  

    %     %% Cluster cross-correlation figure
    %     if run_figCROSS
    %         if sum(clus ~= 0) > 1
    %             [~,~,~] = mkdir([pwd '\klustest\' sdata.combined_name '\figures\']);
    %             figfile = [pwd '\klustest\' sdata.combined_name '\figures\' sdata.combined_name '_E' num2str(tet) '_cross_correlogram.png'];            
    %             figCROSS(sdata,tet,figfile,fig_vis); % overall figure with ratemaps etc
    %         end
    %     end
    % 
    %     %% Cluster space figure
    %     if run_figCSPACE
    %         [~,~,~] = mkdir([pwd '\klustest\' sdata.combined_name '\figures\']);
    %         figfile = [pwd '\klustest\' sdata.combined_name '\figures\' sdata.combined_name '_E' num2str(tet) '_cluster_space.png'];
    %         figCSPACE(sdata.(tetstr).fetdata,figfile,fig_vis); % overall figure with ratemaps etc
    %     end
    end % this ends the electrode loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the data and finish up
%     save([pwd '\klustest\' cname '\' cname '_fdata.mat'],'fdata'); % save field data
    save([pwd '\klustest\' cname '\' cname '_sdata.mat'],'sdata'); % save cluster data
    save([pwd '\klustest\' cname '\' cname '_pdata.mat'],'pdata'); % save session data
    save([pwd '\klustest\' cname '\' cname '_tdata.mat'],'tdata'); % save cluster data

    % finish up
    toc1 = toc/60;
    disp(sprintf('klustest has finished. It took %0.3f seconds or %0.3f minutes',toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',[pwd '\klustest\' cname '\figures'],' &'');">','figures folder','</a>'])
    disp('-------------------------------------------------------------------------------------');

end % ends the main klustest function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A stupid subfunction for preallocating a table
function tin = addToTable(tin,varargin)
    for i=1:2:length(varargin)
        tin.(varargin{i}) = varargin{i+1}; % table in (variable name) = variable value to assign
    end
end % ends addToTable function




































