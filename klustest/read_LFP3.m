 function [lfp,t,Fs] = getLFP3(fname,ds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION  short descr.
% long descr.
%
% USAGE:
%         [out] = template(in,in2)
%
% INPUT:
%         fname - the filename of the LFP to read, can be a full file path
%         ds - the desired sampling rate of the lfp data
%
% OUTPUT:
%    lfp - the lfp values (not in Volts)
%    t - a vector of time values
%    Fs - the sampling rate of the lfp
%
% EXAMPLES:
%
% See also: KLUSTEST FUNCTION3

% HISTORY:
% version 1.0.0, Release 26/02/16 initial release of getEEG and getEGF
% version 2.0.0, Release 11/01/17 created as a faster and more compact alternative to getEEG and getEGF
% version 2.1.0, Release 01/08/17 took advantage of some shortcuts in a function Jim sent, added downsampling
% version 2.2.0, Release 07/08/18 updated and cleaned
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2018 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
if ~exist('ds','var') || isempty(ds) || all(isnan(ds))
   ds = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
    % determine the file type
    [~,~,ext] = fileparts(fname);

    % if the file is .egf use int16 (2 bytes er sample), otherwise use int8 (1 byte per sample)
    if strcmp(ext,'.eeg')
        rform = 'int8';
    else
        rform = 'int16';
    end

    % get data headers and find the first line of the actual data
    [h,d] = get_dacq_headers(fname);
    Fs = h.sample_rate; % the sampling rate
    nsamp = h.num_EEG_samples;

    % try and open the file for reading
    fid = fopen(fname,'r');
    if fid < 0
       error('ERROR: Could not open file: %s... exiting',fname);
    end

    % move through the file to the correct line (actualy read the line immediately before the 'data_start' line
    for t = 1:d-1
       tmp = fgetl(fid);
    end 

    % read the actual data
    fseek(fid,10,0); % move forward 10 characters to remove 'data_start' text
    lfp = fread(fid,nsamp,rform);
    t = (0:length(lfp)-1)'/Fs;
    fclose(fid);

    % downsample data if required
    if ds
        [lfp,t] = resample(lfp,t,ds,'pchip'); % use interpolation and an anti-aliasing filter to resample the signal at a uniform sample rate
        Fs = ds; % the new sampling rate
    end
    lfp = double(lfp(:));
    t = double(t(:));
    Fs = double(Fs);














