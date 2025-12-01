function [Xs,Ys,Zs,spikemap,dwellmap,ratemap] = map_hd_3d(pos,spk,varargin)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short desc.
% long description
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
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
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
    addRequired(p,'pos',@(x) ~isempty(x) && ~all(isnan(x(:))) && size(x,2)==3); 
    addRequired(p,'spk',@(x) size(x,2)==3); 
    addParameter(p,'bins',64,@(x) isnumeric(x) && isscalar(x));       
    addParameter(p,'sigma',10,@(x) isnumeric(x) && isscalar(x)); 
    addParameter(p,'dwellmap',{},@(x) isnumeric(x));       
    parse(p,pos,spk,varargin{:});
    config = p.Results;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % generate dwell time map, if one was not provided to speed up calculation
    if isempty(config.dwellmap)
        [pox,poy,poz] = sph2cart(deg2rad(config.pos(:,1)),deg2rad(-config.pos(:,2)),ones(size(config.pos(:,1)))); % elevation is inverted so pitch up = positive angle   
        [~,~,~,dwellmap] = projectEIGS([pox,poy,poz],config.bins,config.sigma);
    else
        [pox,poy,poz] = sph2cart(deg2rad(config.pos(:,1)),deg2rad(-config.pos(:,2)),ones(size(config.pos(:,1)))); % elevation is inverted so pitch up = positive angle           
        dwellmap = config.dwellmap;
    end
    
    % generate spike map, this will need to be generated every time
    if isempty(config.spk(:,1))
        spikemap = zeros(size(dwellmap));
        ratemap = NaN(size(dwellmap));
        Xs = NaN;
        Ys = NaN;
        Zs = NaN;
    else
        [spx,spy,spz] = sph2cart(deg2rad(config.spk(:,1)),deg2rad(-config.spk(:,2)),ones(size(config.spk(:,1)))); % elevation is inverted so pitch up = positive angle    
        [Xs,Ys,Zs,spikemap] = projectEIGS([spx,spy,spz],config.bins,config.sigma);
        
        % calculate ratemap from the spikemap and dwellmap (both of which are KSDEs not histograms)
        spikemap = spikemap ./ sum(spikemap(:),'omitnan') .* numel(spx);
        dwellmap = dwellmap ./ sum(dwellmap(:),'omitnan') .* numel(pox) .* (1/50);
        ratemap = spikemap ./ dwellmap;        
    end

% figure
% plot3(pox,poy,poz,'k'); hold on;
% plot3(spx,spy,spz,'r.','MarkerSize',30);
% daspect([1 1 1])
% keyboard


    
    if 0
        figure
        subplot(2,2,1)
            s = surf(Xs,Ys,Zs,ratemap);
            s.EdgeColor = 'k';
            s.EdgeAlpha = 0.0;
            axis on vis3d; 
            rotate3d on;
            view(3);    
            daspect([1 1 1]); 
            grid on
            axis tight on       
            colormap(gca,inferno);
            caxis([0 nanmax(ratemap(:))])

        subplot(2,2,2)  
            s = surf(Xt,Yt,Zt,ratemap);
            s.EdgeColor = 'k';
            s.EdgeAlpha = 0.0;
            axis on vis3d; 
            rotate3d on;
            view(3);    
            daspect([1 1 1]); 
            grid on
            axis tight on       
            colormap(gca,inferno);
            caxis([0 nanmax(ratemap(:))])  

        subplot(2,2,3)
            [A,E,~] = cart2sph(Xs,Ys,Zs);
            [Xt,Yt,Zt] = sph2cart(A,E,ratemap./nanmax(ratemap(:)));       
            s = surf(A,E,zeros(size(A)),ratemap);
            s.EdgeColor = 'k';
            s.EdgeAlpha = 0.0;
            axis on vis3d; 
            rotate3d on;
            view(0,90);    
            daspect([1 1 1]); 
            grid on
            axis tight on       
            colormap(gca,inferno);
            caxis([0 nanmax(ratemap(:))])         
    end


































