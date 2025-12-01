function [clumaa,posdata] = PIT_3D_heading(config,clumaa,posdata,overwrite)
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
    if ~exist('overwrite','var') || isempty(overwrite) || isnan(overwrite)
        overwrite = 0;
    end
    min_speed = 5*10;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> 3d HEAD direction dwell maps
    if ~any(ismember(posdata.Properties.VariableNames,'hd_3d_dwell')) || overwrite % if the column(s) do not exist yet
        disp(sprintf('\t3D HD maps...'))
        
        posdata.hd_3d_dwell = cell(size(posdata,1),3); % preallocate

        for ss = 1:size(posdata,1) % for each recording session
            disp(sprintf('\tSession %d of %d (%.f%%)',ss,size(posdata,1),ss/size(posdata,1)*100))
            session_times = posdata.session_times{ss};
            pos = posdata.pos{ss};
    
            for pp = 1:size(session_times,1) % for each part
                idx = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part
                ppox = pos.pox(idx);
                ppoy = pos.poy(idx);
                ppoz = pos.poz(idx);       
                ppit = pos.pitch(idx);
                pyaw = pos.yaw(idx);   

                npoints = numel(pyaw); % number of points
                [X,Y,Z] = sph2cart(deg2rad(pyaw),deg2rad(ppit),ones(npoints,1)); % convert to XYZ
                [Xd,Yd,Zd,F1] = projectEIGS([X,Y,Z],64,10); % project these onto sphere     
                F1 = (F1 ./ sum(F1(:))) .* npoints .* (1/50);
                posdata.hd_3d_dwell(ss,pp) = { single(F1) };
                % negative pitch will be in lower y-axis/row values
            end
        end

        clumaa.Properties.CustomProperties.posdata = posdata;
        cname = [config.data_out_dir '\clumaa.mat'];
        save(cname,'clumaa'); % save session data        
    end

%% >>>>>>>>>> 3d MOVEMENT direction dwell maps
    if ~any(ismember(posdata.Properties.VariableNames,'mov_3d_dwell')) || overwrite
        disp(sprintf('\t3D movement maps...'))
        
        posdata.mov_3d_dwell = cell(size(posdata,1),3); % preallocate

        for ss = 1:size(posdata,1) % for every session
            disp(sprintf('\tSession %d of %d (%.f%%)',ss,size(posdata,1),ss/size(posdata,1)*100))
            pos = posdata.pos{ss};
            pos.pod = NaN(size(pos,1),1);
            session_times = posdata.session_times{ss};   

            for pp = 1:3 % for every session           
                pidxn = pos.pot > session_times(pp,1) & pos.pot < session_times(pp,2); % index for position data in this part

                % first half map
                pox = pos.pox_planar(pidxn);
                poy = pos.poy_planar(pidxn);
                poz = pos.poz(pidxn);

                % calculate change in position along each axis
                [speed,azimuth,tilt,dxyz] = insta_speed(pox,poy,'poz',poz,'srate',50,'wsize',0.25);
                pos.pod(pidxn) = azimuth;
                V = normr( dxyz );

%                 pxd = padarray(pox,25,NaN,'both');
%                 pxd = [circshift(pxd(:),25) circshift(pxd(:),-25)];
%                 pxd = pxd(25:end-26,1)-pxd(25:end-26,2);
%                 pxd = pxd./50;
%                 pyd = padarray(poy,25,NaN,'both');
%                 pyd = [circshift(pyd(:),25) circshift(pyd(:),-25)];
%                 pyd = pyd(25:end-26,1)-pyd(25:end-26,2);
%                 pyd = pyd./50;
%                 pzd = padarray(poz,25,NaN,'both');
%                 pzd = [circshift(pzd(:),25) circshift(pzd(:),-25)];
%                 pzd = pzd(25:end-26,1)-pzd(25:end-26,2);
%                 pzd = pzd./50;
% 
%                 % normalise these vectors
%                 V = [pxd(:) pyd(:) pzd(:)];
%                 V = normr(V);
%                 [podn,~] = cart2pol(V(:,1),V(:,2)); 
% keyboard
%                 pos.pod(pidxn) = rad2deg(podn); % movement direction   
% 
%                 % Determine average speed in every heading 
%                 % Calculate instantaneous speed
%                 pod = sqrt(sum([pxd pyd pzd].^2,2)) / 1000; % distance between each position point
%                 pov = pod' .* 50; % calculate speed in metres/second
            
                % determine dwell time at every elevation and azimuth using spherical coordinates
                [X,Y,Z,dwellmap] = projectEIGS(V(speed>min_speed,:),64,10);
                posdata.mov_3d_dwell(ss,pp) = { single(dwellmap) };

                if 0
                    figure
                    imagesc([-180 180],[-90 90],dwellmap,'alphadata',~isnan(dwellmap))
                    xlabel(sprintf('Azimuth (%c)',176))
                    ylabel(sprintf('Pitch (%c)',176))
                    axis xy tight
                    daspect([1 1 1])
                    view(0,90);
                    set(gca,'FontSize',12)
                    ax = gca;
                    ax.CLim = [0,max(dwellmap(:))];
                    colormap(gca,turbo)
                    ax.XTick = -180:90:180;
                    ax.YTick = -90:45:90;
                    keyboard
                end
            end
        end

        clumaa.Properties.CustomProperties.posdata = posdata;
        cname = [config.data_out_dir '\clumaa.mat'];
        save(cname,'clumaa'); % save session data
    end
   






















