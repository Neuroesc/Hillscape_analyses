%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
%GIT_audit  main analysis script for:
% Grieves, Duvelle and Taube (202X) 
% How the brain represents irregular terrain: hippocampal place cells make mountains out of (rat)hills
%
% USAGE:
%           % process with default settings
%           GIT_audit 
%
% EXAMPLES:
%
%           % run function using default values:
%           GIT_audit
%
% See also: 

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
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Overwrite settings
    overwrite_sdata                                     = 0; % generate the sdata file from scratch
    maintain_sdata                                      = 1; % if the sdata table is loaded into memory don't reload or regenerate it
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Settings for different analyses
    config = struct;
    % You will need to edit config.main_dir to reflect the directory where you want data/figures to be saved
    config.main_dir                                     = 'C:\Users\roddyg\OneDrive - University of Glasgow\Projects in prep\2021 Place cells irregular terrain\';
    % config.main_dir                                     = 'C:\Users\roddy\OneDrive - Dartmouth College\Projects in prep\2021 Place cells irregular terrain';
    cd(config.main_dir);
    
    config.data_dir                                     = 'X:\2021 HPC irregular terrain project'; % the raw (rat) data directories can be found here
    config.data_out_dir                                 = [config.main_dir '\associated data\']; % output data can be saved here
    config.sname                                        = [config.data_dir '\cluma.mat']; % location of the accumulated sdata mat file  
    config.media_out_dir                                = [config.main_dir '\associated media\outputs']; % figure directory 1, for figure parts and general media
    config.fig_dir                                      = [config.main_dir '\main figures']; % figure directory 2, for main, final figures

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Functions to run
    funconfig = struct;

    % Func. to run some code on every data directory
    if 0
        run_fun_for_audit(config); 
        return
    end
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Figure settings
    print_now = 1; % should the function print the figures? If 1 the figures will be set to invisible
    print_qual = '-r200'; % quality for printing in DPI   
    print_form = '-tiff'; % format for printing figures
    close all;
    fast_figs = 0;
    res = 350;
    
    % figure colors
    plot_set = {rgb('DarkSlateGray'),rgb('Orange'),rgb('DarkSlateGray'),rgb('Cyan');'<','o','>','s'};
    maze_names = {'Arena 1','Hillscape','Arena 2', 'A1','Hs','A2','Arena'};
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Main figures
    overwrite_clumaa = 0;
    if ~exist('clumaa','var') || overwrite_clumaa
        load('C:\Users\roddyg\OneDrive - University of Glasgow\Projects in prep\2021 Place cells irregular terrain\associated data\clumaa.mat','clumaa')
        posdata = clumaa.Properties.CustomProperties.posdata;    
    end

    % cell indices
    pidx = clumaa.cell_type(:,2)==3 & clumaa.planar_spatial_info_shuffles(:,2)>2 & clumaa.f_rate_hz(:,1)>0.1;
    
%% Functions to run
    funconfig = struct;
    funconfig.PIT_remapping                         = 0; % stability between mazes, remapping between mazes
    funconfig.PIT_localisation                      = 0; % field size and number
    funconfig.PIT_elongation                        = 0; % alignment to geometry/hills
    funconfig.PIT_repetition                        = 0; % field repetition
    funconfig.PIT_field_anisotropy                  = 0; % field anisotropy and the relation to behaviour    
    funconfig.PIT_testing                           = 0;

    % paper figures v1
    funconfig.PIT_fig_1_v1                          = 0; % Fig 1: Schematic of hypotheses, mazes, example cells
    funconfig.PIT_fig_2_v1                          = 0; % Fig 2: Correlations, place fields per cell, field radius
    funconfig.PIT_fig_3_v4                          = 0; % Fig 3: Field elongation, field orientation, anisotropy analysis
    funconfig.PIT_fig_5_v1                          = 0; % Fig 4: Repetition example cells, autocorrelations, repetition score    
    funconfig.PIT_fig_8_v1                          = 0; % Fig 5: BVC results, Field to wall, angle plots

    funconfig.PIT_sup_fig_16_v1                     = 0; % Fig SX: Maze schematics and photos    
    funconfig.PIT_sup_fig_5_v1                      = 0; % Fig S1: Rat and cell stats, histology  
    funconfig.PIT_sup_fig_1_v3                      = 1; % Fig S2: Plot example place cells        
    funconfig.PIT_sup_fig_3_v1                      = 0; % Fig S3: Correlations between arena and hills restricted to edges, corners and alignment 
    funconfig.PIT_sup_fig_4_v1                      = 0; % Fig S4: Plot place cells, ranked by similarity between the two mazes 
    funconfig.PIT_sup_fig_13_v1                     = 0; % Fig S5: Field to wall, distance plots
    funconfig.PIT_sup_fig_14_v1                     = 0; % Fig S6: Field to wall, angle plots
    funconfig.PIT_sup_fig_8_v1                      = 0; % Fig S7: Example trajectory, movement direction and head direction sampling   
    funconfig.PIT_fig_7_v1                          = 0; % Fig S8 & S9: Directionality and Pitch coding vs repetition score
    funconfig.PIT_sup_fig_15_v1                     = 0; % Fig S13: Plot example BVC model place cells        
    funconfig.PIT_sup_fig_11_v1                     = 0; % Fig SX: Scatter graphs  
    funconfig.PIT_over_glm                          = 0; % Fig SX: Reviewer requested, GLM, residuals and stats
    funconfig.PIT_sup_1sess_v1                      = 0; % Fig SX: Reviewer requested, single session per rat

    % unused code
    funconfig.PIT_sup_fig_2_v1                      = 0; % UNUSED: Non repetitive example cells, field distribution all place cells
    funconfig.PIT_sup_fig_7_v1                      = 0; % UNUSED: behavioural anisotropy vs field anisotropy 
    funconfig.PIT_sup_fig_6_v1                      = 0; % UNUSED: behaviour coverage
    funconfig.PIT_fig_6_v1                          = 0; % UNUSED: behaviour figure, movement bias, rearing    
    funconfig.PIT_sup_fig_1_v1                      = 0; % UNUSED: Plot all place cells
    funconfig.PIT_sup_fig_9_v1                      = 0; % UNUSED: schematic of simulation
    funconfig.PIT_sup_fig_10_v1                     = 0; % UNUSED: height vs slope
    funconfig.PIT_sup_fig_12_v1                     = 0; % UNUSED: ?  

    % additional analyses


%% run the requested analyses
    fnames = fieldnames(funconfig);
    for ff = 1:length(fnames)
        val = funconfig.(fnames{ff});
        if val
            disp(sprintf('\trunning: %s',fnames{ff}))
            eval([fnames{ff} ';']);
        end
    end

    % save concatenated data
    if 0
        clumaa.Properties.CustomProperties.posdata = posdata;
        cname = [config.data_out_dir '\clumaa.mat'];
        save(cname,'clumaa'); % save session data
    end







































