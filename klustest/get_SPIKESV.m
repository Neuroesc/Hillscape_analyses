function [ts,ch1,ch2,ch3,ch4,hdata] = get_SPIKESV(fname,bfast)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getSPIKESV  function to import waveform data, quite quickly, from tetrode files
% This function acts as a wrapper for getspikes, whilst also converting the waveform data 
% to microvolts using the method outlined at the end of this file (supplied by Jim)
% This function can also just load the spike times for a marginal speed increase 
% (the waveforms still have to loaded but not processed or converted)
%
% USAGE:
%         [ts,ch1,ch2,ch3,ch4,hdata] = getSPIKESV(fname,bfast)
%
% INPUT:
%         fname - filename of .eeg or .egf file to convert, i.e. '180119a_square.egf', full file paths should also work
%         bfast - (optional) if set to 1 the function will only output the spike times, ts, the other outputs will be set to NaN
%
% OUTPUT:
%    ts - the spike times (in seconds)
%    ch1-ch4 - the spike waveforms in microvolts, each row is a spike, 50 columns, saved as single
%    hdata - the set file header data
%
% EXAMPLES:
% load in spike waveforms for a tetrode in microvolts:
% [ts,ch1,ch2,ch3,ch4,hdata] = getSPIKESV('180417a_square.9',0);
%
% load in spike times for a tetrode without the waveform information:
% ts = getSPIKESV('180417a_square.9',1);
%
% See also: getspikes getLFPV getDACQDATAHEADERS

% HISTORY:
% version 1.0.0, Release 02/04/19 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT ARGUMENTS CHECK
if ~exist('bfast','var') || isempty(bfast) || all(isnan(bfast))
   bfast = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION BODY
% load the spike data from the file
% [ts,ch1,ch2,ch3,ch4] = getspikes(fname,bfast);
[a,b,c] = fileparts(fname);
if isempty(a)
    a = pwd;
end
tet = str2double(regexp(c,'\d+','match'));

if bfast % if bfast, we don't want waveforms, ch1-ch4 will be NaN and conversion to uV is not necessary
    ts = Nlx2MatSpike([a '\' b '\TT' num2str(tet) '.ntt'],[1 0 0 0 0],0,1,1);    
    return
else
    if exist([a '\' b '\TT' num2str(tet) '.ntt'],'file')
        [ts,chs,hdata] = Nlx2MatSpike([a '\' b '\TT' num2str(tet) '.ntt'],[1 0 0 0 1],1,1,1);    
    else
        [ts,chs,hdata] = Nlx2MatSpike([a '\TT' num2str(tet) '.ntt'],[1 0 0 0 1],1,1,1);            
    end
end

% get the data for this session
idx = find(contains(hdata,'ADBitVolts'));
adcline = hdata{idx}(13:end);
adbv = str2double(split(adcline,' ')); % ADBitVolt value

% extract and convert waveforms to microvolts
chs = single( chs .* adbv' .* 10^6 );

ch1 = squeeze(chs(:,1,:))';
ch2 = squeeze(chs(:,2,:))';
ch3 = squeeze(chs(:,3,:))';
ch4 = squeeze(chs(:,4,:))';













