function hill_field_analysis(data_dir,rnow,dnow)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short description
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
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 21/02/23 Initial release

%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK  
    outname = 'klustest';
    skipfigs = 1; % 1 = skip making a figure if it exists already    
    save_figs = 1; % 1 = make and save figures (VERY time consuming, takes about 90% of klustest's time)
    fast_figures = 1; % 1 = 4-5x faster figure saving, but at a slightly lower quality
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'); tic;
    stk = dbstack;
    tnow = datestr(now,'yyyy-mm-dd-HH-MM-SS');
    disp(sprintf('Running %s at %s...',stk.name,tnow))
    
%% >>>>>>>>>> Prepare the data
    sname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    disp(sprintf('Loading cluma: %s',sname))    
    load(sname,'cluma'); % load saved session data
    pdata = cluma.Properties.CustomProperties.pdata;
    part_names = {'arena1','hills','arena2'};
    nparts = size(pdata.session_times,1);
    disp(sprintf('\t...%d sessions',nparts));   
    disp(sprintf('\t...recording date: %s',dnow));            
    disp(sprintf('\t...done'))    
    
    %% Map settings    
    mapset = pdata.mapset;
    pos = pdata.pos;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Run through clusters
    disp(sprintf('Analysing clusters...'))
    ucis = unique(cluma.uci); % list of unique cells in sdata    
    
    for pp = 1:nparts % for every part
        part_now = part_names{pp};
        disp(sprintf('\t%s',part_now))      
            
        loopout = looper(length(ucis));
        for uu = 1:length(ucis)
            uci = ucis{uu};
            idx = find( ismember(cluma.uci,uci) & cluma.partn==pp );
            % session_times = pdata.session_times;
            % pidx = pos.pot > session_times(pp,1) && pos.pot < session_times(pp,2); % index for position data in this part

%% >>>>>>>>>> Get (planar) firing rate map etc
            mapset.zcut = 1.2; % z-score for field detection in z-scored firing rate maps   
            mapset.frcut = 1; % place fields must have a peak firing rate at least greater than this value
            mapset.arcut = 400; % cm2, place fields must have a total area greater than this value   
            pdata.mapset = mapset;
            
            % planar map
            map = cluma.ratemap_planar{idx};
            zmap = ( map - mean(map(:),'omitnan') ) / std(map(:),'omitnan'); % zscore ratemap
            thresh_ratemap = imbinarize(zmap,mapset.zcut); % 2 s.d. threshold

            datout = regionprops('table',thresh_ratemap,map,'Area','Centroid','WeightedCentroid','MajorAxisLength','MinorAxisLength','Orientation','MaxIntensity','BoundingBox','PixelIdxList'); % detect contiguous regions
            if ~isempty(datout) % filter out small or low firing regions
                datout.Area(:) = datout.Area(:) .* ((mapset.binsize/10)^2); % convert field area to cm2
                nindx = datout.Area(:) < mapset.arcut | datout.MaxIntensity(:) < mapset.frcut;                        
                datout(nindx,:) = [];      
            end
            if ~any(ismember(cluma.Properties.VariableNames,'planar_fields')) % if the column(s) do not exist yet
                cluma.planar_fields = NaN(size(cluma,1),1); % preallocate
                cluma.planar_field_data = cell(size(cluma,1),1); % preallocate 
                cluma.planar_field_area = NaN(size(cluma,1),2); % preallocate                
            end               
            cluma.planar_fields(idx) = size(datout.Area,1);
            cluma.planar_field_data(idx) = { datout };             
            cluma.planar_field_area(idx,:) = [sum(datout.Area(:),'all','omitmissing') sum(datout.Area(:),'all','omitmissing')./numel(map(:))];

            % surficial map
            if ~any(ismember(cluma.Properties.VariableNames,'surficial_fields')) % if the column(s) do not exist yet
                cluma.surficial_fields = NaN(size(cluma,1),1); % preallocate
                cluma.surficial_field_data = cell(size(cluma,1),1); % preallocate   
                cluma.surficial_field_area(idx,:) = NaN(size(cluma,1),2); % preallocate               
            end             
            if pp==2
                map = cluma.ratemap_surficial{idx};
                zmap = ( map - mean(map(:),'omitnan') ) / std(map(:),'omitnan'); % zscore ratemap
                thresh_ratemap = imbinarize(zmap,mapset.zcut); % 2 s.d. threshold
    
                datout = regionprops('table',thresh_ratemap,map,'Area','Centroid','WeightedCentroid','MajorAxisLength','MinorAxisLength','Orientation','MaxIntensity','BoundingBox','PixelIdxList'); % detect contiguous regions
                if ~isempty(datout) % filter out small or low firing regions
                    datout.Area(:) = datout.Area(:) .* ((mapset.binsize/10)^2); % convert field area to cm2
                    nindx = datout.Area(:) < mapset.arcut | datout.MaxIntensity(:) < mapset.frcut;                        
                    datout(nindx,:) = [];      
                end              
                cluma.surficial_fields(idx) = size(datout.Area,1);
                cluma.surficial_field_data(idx) = { datout }; 
                cluma.surficial_field_area(idx,:) = [sum(datout.Area(:),'all','omitmissing') sum(datout.Area(:),'all','omitmissing')./numel(map(:))];                
            end

            loopout = looper(loopout);            
        end                     
    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% Save the data and finish up
    cluma = rmprop(cluma,'pdata');
    cluma = addprop(cluma,{'pdata'},{'table'});
    cluma.Properties.CustomProperties.pdata = pdata;
    cname = [data_dir '\' rnow '\' dnow '\cluma\cluma.mat'];
    save(cname,'cluma'); % save session data
    analysis_log({'hill_field_analysis'},1,'version',{'v1.0.0'});
    
    % finish up
    toc1 = toc/60;
    disp(sprintf('%s has finished. It took %0.3f seconds or %0.3f minutes',stk.name,toc,toc1)) % Stop counting time and display results
    disp(['Go to ','<a href = "matlab: [s,r] = system(''explorer ',pwd,' &'');">','current folder','</a>'])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');























