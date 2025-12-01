function s = PIT_spatial_selectivity(r,d)
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
minTime = 0;
d(isnan(d)) = eps;
r(isnan(r)) = eps;

% Compute the spatial specificity of the map, based on the formula proposed by Skaggs et al. (1993):
%  specificity = SUM { p(i) . lambda(i)/lambda . log2(lambda(i)/lambda) }
T = sum(d(:),'omitnan');
p_i = d/(T+eps); % Probability of the animal occupying bin 'i'
lambda_i = r;
lambda = lambda_i(:)'*p_i(:);

if T == 0 || lambda == 0
    s = 0;
else
    ok = d>minTime; 
    s = sum(sum(p_i(ok).*lambda_i(ok)/lambda.*log2(lambda_i(ok)/lambda),'omitnan'),'omitnan');
end




























