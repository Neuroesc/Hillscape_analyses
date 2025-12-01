% function out = templateNV(in,varargin)
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
%     p = inputParser;
%     addRequired(p,'in',@(x) ~isempty(x) && ~all(isnan(x(:))));  
%     addOptional(p,'param2',def_param1,@(x) isnumeric(x) && isscalar(x)); 
%     addParameter(p,'param3',def_param2,@(x) isnumeric(x) && isscalar(x));   
%     parse(p,in,varargin{:});
%     config = p.Results;
% close all;

% if planar:
% fields per m2 should be the same in arena and hills
% grid score should be higher in planar projection
% ellipticity should be higher in surficial projection
% surficial ellipse orientation should be horizontal (90/270 when using get_grid_score)

% if surficial:
% fields per m2 should be lower in arena than hills
% grid score should be higher in surficial projection
% ellipticity should be high in planar projection
% planar projection ellipse orientation should be vertical (180/0 when using get_grid_score)

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
close all;

    if 1
        ncells = 100;
        spacing_specs = [16 8 4 24]; % smaller = larger grid fields/scale

        dat = struct;    
        simdata = table;
        warning('off','MATLAB:table:RowsAddedExistingVars');
        for cc = 1:ncells
            %% Generate the basic hexagonal grid map
            simdata.cell(cc,1) = cc;

            k = round( normrnd(spacing_specs(1),spacing_specs(2),1,1) ); % grid spacing
            if k < spacing_specs(3)
                k = spacing_specs(3);
            end
            if k > spacing_specs(4)
                k = spacing_specs(4);
            end        
% k = 10            
            simdata.grid_spacing(cc,1) = k;
            disp(sprintf('\tcell %d (%.1f%%) spacing: %.1f',cc,cc/ncells*100,k))

            n = 1000; % full size grid edge length (we will cut to the middle to avoid edge effects)
            [Xp,Yp] = meshgrid(1:n,1:n); 
            I0 = 1/2*(1+sin(2*pi/n*k*Xp)); 
            I60 = imrotate(I0,60,'nearest','crop');
            I120 = imrotate(I0,120,'nearest','crop');
            large_grid = I0 .* I60 .* I120;
            large_grid(large_grid<0.15) = 0;
            large_grid = imrotate(large_grid,14,'nearest','crop');
            [Xr,Yr] = meshgrid((1:size(large_grid,2))-size(large_grid,2)/2,(1:size(large_grid,1))-size(large_grid,1)/2);
            orientation = rand(1)*360;    
            Fr = imrotate(large_grid,orientation,'nearest','crop');
            simdata.grid_orientation(cc,1) = orientation;    

            %% Load position data
            load(['C:\Users\F004KS7\OneDrive - Dartmouth College\Matlab\Functions\Grid cells irregular terrain\Simulations\path_data_1.mat'],'pos','pos_flat_curve','pos_planar','pos_surficial','fitresult','-mat');

% figure
% imagesc(large_grid)
% daspect([1 1 1])
% axis xy off
% keyboard
% exportgraphics(gcf,['C:\Users\F004KS7\Downloads\fig_grid.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  

% figure('visible','on','Units','pixels','Position',[10, 100, 800, 800]); % open/create figure
% plot3(pos(:,1),pos(:,2),pos(:,3),'k'); hold on;
% daspect([1 1 1])
% axis xy off
% view(3)   

figure('visible','on','Units','pixels','Position',[10, 100, 800, 800]); % open/create figure
subplot(2,3,2)
c = cline(pos(:,1),pos(:,2),pos(:,3),pos(:,3)); hold on;
set(c,'LineWidth',1.5)
daspect([1 1 1])
colormap('turbo')
caxis([50 430])
axis xy off
view(3)   

subplot(2,3,4)






keyboard
% exportgraphics(gcf,['C:\Users\F004KS7\Downloads\fig_pos.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  
%         
% 
figure('visible','on','Units','pixels','Position',[10, 100, 800, 800]); % open/create figure
plot3(pos_surficial(:,1),pos_surficial(:,2),pos_surficial(:,3),'k'); hold on;
daspect([1 1 1])
axis xy off
view(3)   

keyboard
% exportgraphics(gcf,['C:\Users\F004KS7\Downloads\fig_pos_surf.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  
%     
        
            
            
            ppos = pos_planar ./ 10;
            ppos = ppos - min(ppos,[],1,'omitnan');
            spos = pos_surficial ./ 10;
            spos = spos - min(spos,[],1,'omitnan');

        %     Fr = ( Fr./max(Fr(:)) ) .* 10;                    
        %     Vp = interp2(Xr,Yr,Fr,ppos(:,1),ppos(:,2),'nearest',0);
        %     Vs = interp2(Xr,Yr,Fr,spos(:,1),spos(:,2),'nearest',0); 

            % arena spikes
            coeff1 = 0.5; % weighting for pmap probability, bigger = more spikes
            denom1 = 0; % denominator for spike lag, bigger = smaller lag, zero disables
            Fr = ( Fr./max(Fr(:)) ) .* coeff1;                    
            Va = interp2(Xr,Yr,Fr,ppos(:,1)+normrnd(0,50,1,1),ppos(:,2)+normrnd(0,20,1,1),'nearest',0);        
            ind = 1:length(ppos(:,1));
            noiselam = 0; % lambda for poisson noise, higher = more noise          
            spk_prob = poissrnd(Va) + poissrnd(noiselam,size(Va));
            spk_prob(spk_prob<0) = 0;
            ind2 = repelem(ind,spk_prob);    
            aspk = pos_flat_curve(ind2,:); % planar (grid when seen from above)

            % planar spikes                  
            Vp = interp2(Xr,Yr,Fr,ppos(:,1),ppos(:,2),'nearest',0);
            spk_prob = poissrnd(Vp) + poissrnd(noiselam,size(Vp));
            spk_prob(spk_prob<0) = 0;
            ind2 = repelem(ind,spk_prob);    
            pspk = pos_flat_curve(ind2,:); % planar (grid when seen from above)
            pspk2 = pos_surficial(ind2,:); %   

            % surficial spikes
            Vs = interp2(Xr,Yr,Fr,spos(:,1),spos(:,2),'nearest',0);    
            spk_prob = poissrnd(Vs) + poissrnd(noiselam,size(Vs));
            spk_prob(spk_prob<0) = 0;
            ind2 = repelem(ind,spk_prob);    
            sspk = pos_flat_curve(ind2,:); % surficial (grid when curve is stretched and flattened) 
            sspk2 = pos_surficial(ind2,:); %      

            %% Map data
            % planar (vertical projection) ratemaps     
            % Map settings
            mapset.ppm          = 1000;
            mapset.method       = 'histogram';
            mapset.binsize      = 32; % (mm) firing rate map bin size
            mapset.ssigma       = 64; % (mm) firing rate map smoothing sigma
            mapset.padding      = 10; % (mm) how much to pad the edges of the firing rate map
            mapset.mindwell     = 0.01; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
            mapset.mindist      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
            mapset.smethod      = 1; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing              
            mapset.frcut        = 0.2; % cutoff % for regions to be considered a place field
            mapset.arcut        = 36; % cutoff area for regions to be considered a place field
    
            pos_all = {pos_planar, pos_flat_curve, pos_surficial, pos_flat_curve, pos_surficial};
            spk_all = {aspk, pspk, pspk2, sspk, sspk2};
            nmes = {'arena','planar_planar','planar_surficial','surficial_planar','surficial_surficial'};
            plotting_coords = cell(size(nmes));
            
            for ii = 1:5     
                % firing rate map
                pos = pos_all{ii}(:,1:2);
                spk = spk_all{ii}(:,1:2);     
                rmset = mapset;
                rmset.maplims = [min(pos(:,1),[],'omitnan') min(pos(:,2),[],'omitnan') max(pos(:,1),[],'omitnan') max(pos(:,2),[],'omitnan')];
                if isfield(dat,[nmes{ii} '_speedlift']) % if a previously computed dwellmap exists, use this to save time [data_projection]
                    speedlift = dat.([nmes{ii} '_speedlift']);
                    [rmap0,dmap0, ~, r, ~] = rate_mapper(pos,spk,rmset,speedlift);              
                else
                    [rmap0,dmap0, ~, r, speedlift] = rate_mapper(pos,spk,rmset);              
                    dat.([nmes{ii} '_speedlift']) = speedlift;
                end  
                simdata.([nmes{ii} '_rmap'])(cc,1) = { single(rmap0) };   
                
                % place fields
                thresh_ratemap = imbinarize(rmap0,max([0, mapset.frcut*max(rmap0(:),[],'omitnan')]) );
                datout = regionprops('table',thresh_ratemap,rmap0,'Area','Centroid','WeightedCentroid','MajorAxisLength','MinorAxisLength','Orientation','ConvexHull');
                datout.Area(:) = datout.Area(:) .* ((mapset.binsize/10)^2); % convert field area to cm2
                nindx = datout.Area(:) < mapset.arcut;
                datout(nindx,:) = [];
                simdata.fields(cc,ii) = size(datout.Area,1);
                
                % grid score etc
                amap0 = ndautoCORR(rmap0,rmap0);
                [g0,gdata0] = get_grid_score(amap0,mapset.binsize/10,'method','wills');         
                simdata.([nmes{ii} '_amap'])(cc,1) = { single(amap0) };        
                simdata.gscore(cc,ii) = g0;
                simdata.eccentricity(cc,ii) = gdata0.ellipticity;
                simdata.field_elongation(cc,ii) = gdata0.elongation;
                
                % coordinates for plotting
                x = movmean( r.xgrid+mean(rmset.maplims([1 3])), 2, 'Endpoints','discard');
                y = movmean( r.ygrid+mean(rmset.maplims([2 4])), 2, 'Endpoints','discard');
                [xx,yy] = meshgrid(x,y);
                if ismember(ii,[3 5])
                    zz = zeros(size(xx));
                else
                    zz = -fitresult( xx );  
                    zz = reshape(zz,size(xx));
                end

                plotting_coords(cc,ii) = { cat(3,xx,yy,zz) };        
            end

            
% figure
% t1 = 2;
% t2 = 4;
% t3 = 3;
% t4 = 5;
% 
% subplot(3,2,1)
% pox = pos_all{t1}(:,1);
% poy = pos_all{t1}(:,2);
% poz = pos_all{t1}(:,3);
% spx = spk_all{t1}(:,1);
% spy = spk_all{t1}(:,2);
% spz = spk_all{t1}(:,3);
% plot3(pox,poy,poz,'k'); hold on;
% plot3(spx,spy,spz,'r.');
% daspect([1 1 1])
% axis xy off
% 
% subplot(3,2,2)
% pox = pos_all{t2}(:,1);
% poy = pos_all{t2}(:,2);
% poz = pos_all{t2}(:,3);
% spx = spk_all{t2}(:,1);
% spy = spk_all{t2}(:,2);
% spz = spk_all{t2}(:,3);
% plot3(pox,poy,poz,'k'); hold on;
% plot3(spx,spy,spz,'r.');
% daspect([1 1 1])
% axis xy off
% 
% subplot(3,2,3)
% xx = plotting_coords{cc,t1}(:,:,1);
% yy = plotting_coords{cc,t1}(:,:,2);
% zz = plotting_coords{cc,t1}(:,:,3);
% map = simdata.([nmes{t1} '_rmap']){cc,1};
% surf(xx,yy,zz,map,'EdgeColor','none');
% daspect([1 1 1])
% axis xy off
%             
% subplot(3,2,4)
% xx = plotting_coords{cc,t2}(:,:,1);
% yy = plotting_coords{cc,t2}(:,:,2);
% zz = plotting_coords{cc,t2}(:,:,3);
% map = simdata.([nmes{t2} '_rmap']){cc,1};
% surf(xx,yy,zz,map,'EdgeColor','none');
% daspect([1 1 1])
% axis xy off          
%         
% subplot(3,2,5)
% xx = plotting_coords{cc,t3}(:,:,1);
% yy = plotting_coords{cc,t3}(:,:,2);
% zz = plotting_coords{cc,t3}(:,:,3);
% map = simdata.([nmes{t3} '_rmap']){cc,1};
% surf(xx,yy,zz,map,'EdgeColor','none');
% daspect([1 1 1])
% axis xy off
%             
% subplot(3,2,6)
% xx = plotting_coords{cc,t4}(:,:,1);
% yy = plotting_coords{cc,t4}(:,:,2);
% zz = plotting_coords{cc,t4}(:,:,3);
% map = simdata.([nmes{t4} '_rmap']){cc,1};
% surf(xx,yy,zz,map,'EdgeColor','none');
% daspect([1 1 1])
% axis xy off   
% 
% exportgraphics(gcf,['C:\Users\F004KS7\Downloads\fig_sim1.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  
% 
% 
% 
% 
%            
% figure
% t1 = 2;
% t2 = 4;
% t3 = 3;
% t4 = 5;
% 
% subplot(2,2,1)
% map = simdata.([nmes{t1} '_amap']){cc,1};
% imagesc(map);
% daspect([1 1 1])
% axis xy off
% caxis([-0.2 1])
%             
% subplot(2,2,2)
% map = simdata.([nmes{t2} '_amap']){cc,1};
% imagesc(map);
% daspect([1 1 1])
% axis xy off
% caxis([-0.2 1])
% 
% subplot(2,2,3)
% map = simdata.([nmes{t3} '_amap']){cc,1};
% imagesc(map);
% daspect([1 1 1])
% axis xy off
% caxis([-0.2 1])
%             
% subplot(2,2,4)
% map = simdata.([nmes{t4} '_amap']){cc,1};
% imagesc(map);
% daspect([1 1 1])
% axis xy off
% caxis([-0.2 1])
% 
% exportgraphics(gcf,['C:\Users\F004KS7\Downloads\fig_sim2.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  
% 
% 
% 
% 
% 
% 
% 
% 
% return
% 
% subplot(3,2,5)
% map = simdata.([nmes{t1} '_amap']){cc,1};
% imagesc(map);
% daspect([1 1 1])
% axis xy off
%             
% subplot(3,2,6)
% map = simdata.([nmes{t2} '_amap']){cc,1};
% imagesc(map);
% daspect([1 1 1])
% axis xy off
            
            
            
            
            
            
% t1 = 2;
% t2 = 4;
% figure
% subplot(3,2,1)
% pox = pos_all{t1}(:,1);
% poy = pos_all{t1}(:,2);
% poz = pos_all{t1}(:,3);
% spx = spk_all{t1}(:,1);
% spy = spk_all{t1}(:,2);
% spz = spk_all{t1}(:,3);
% plot3(pox,poy,poz,'k'); hold on;
% plot3(spx,spy,spz,'r.');
% daspect([1 1 1])
% axis xy off
% 
% subplot(3,2,2)
% pox = pos_all{t2}(:,1);
% poy = pos_all{t2}(:,2);
% poz = pos_all{t2}(:,3);
% spx = spk_all{t2}(:,1);
% spy = spk_all{t2}(:,2);
% spz = spk_all{t2}(:,3);
% plot3(pox,poy,poz,'k'); hold on;
% plot3(spx,spy,spz,'r.');
% daspect([1 1 1])
% axis xy off
% 
% subplot(3,2,3)
% xx = plotting_coords{cc,t1}(:,:,1);
% yy = plotting_coords{cc,t1}(:,:,2);
% zz = plotting_coords{cc,t1}(:,:,3);
% map = simdata.([nmes{t1} '_rmap']){cc,1};
% surf(xx,yy,zz,map,'EdgeColor','none');
% daspect([1 1 1])
% axis xy off
%             
% subplot(3,2,4)
% xx = plotting_coords{cc,t2}(:,:,1);
% yy = plotting_coords{cc,t2}(:,:,2);
% zz = plotting_coords{cc,t2}(:,:,3);
% map = simdata.([nmes{t2} '_rmap']){cc,1};
% surf(xx,yy,zz,map,'EdgeColor','none');
% daspect([1 1 1])
% axis xy off          
%         
% subplot(3,2,5)
% map = simdata.([nmes{t1} '_amap']){cc,1};
% imagesc(map);
% daspect([1 1 1])
% axis xy off
%             
% subplot(3,2,6)
% map = simdata.([nmes{t2} '_amap']){cc,1};
% imagesc(map);
% daspect([1 1 1])
% axis xy off
% 
% exportgraphics(gcf,['C:\Users\F004KS7\Downloads\fig.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  
%          keyboard     
        end
    end

    
    
    fig_clust = figure('visible','on','Units','pixels','Position',[50, 50, 1850, 920]); % open/create figure
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    dot_size = [22 30 22 30]*2;
    dot_sigma = 0.12;    
    
    ax1 = axes('Units','pixels','Position',[100,550,120,250]);  
        gdata = simdata.gscore;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,1),gdata(:,2));     
        cols = [rgb('Black');rgb('Orange');rgb('Black');rgb('Orange')];
        marks = {'o','s','o','s'};
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XLim = [0.5 2.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-1 1]; 
        ax.YTick = -1:0.25:1;
        ax.FontSize = 14;        
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular grid score'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6,'plot_omnibus',1,'omnibus_text_y_gap_coeff',1);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))    
    
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+230,ax1.Position(2),ax1.Position(3),ax1.Position(4)]);  
        gdata = simdata.eccentricity;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,1),gdata(:,2));     
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XLim = [0.5 2.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-1 1]; 
        ax.YTick = -1:0.25:1;
        ax.FontSize = 14;        
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular eccentricity'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6,'plot_omnibus',1,'omnibus_text_y_gap_coeff',1);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))    
    
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+230,ax1.Position(2),ax1.Position(3),ax1.Position(4)]);  
        gdata = simdata.fields;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,1),gdata(:,2));     
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XLim = [0.5 2.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-10 10]; 
        ax.YTick = -10:2:10;
        ax.FontSize = 14;
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular total fields'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6,'plot_omnibus',1,'omnibus_text_y_gap_coeff',1);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))    
    

    ax1 = axes('Units','pixels','Position',[100,100,120,250]);  
        gdata = simdata.gscore;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,3),gdata(:,4));     
        cols = [rgb('Black');rgb('Orange');rgb('Black');rgb('Orange')];
        marks = {'o','s','o','s'};
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XLim = [0.5 2.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-1 1]; 
        ax.YTick = -1:0.25:1;
        ax.FontSize = 14;        
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular grid score'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6,'plot_omnibus',1,'omnibus_text_y_gap_coeff',1);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))    
    
    ax2 = axes('Units','pixels','Position',[ax1.Position(1)+230,ax1.Position(2),ax1.Position(3),ax1.Position(4)]);  
        gdata = simdata.eccentricity;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,3),gdata(:,4));     
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XLim = [0.5 2.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-1 1]; 
        ax.YTick = -1:0.25:1;
        ax.FontSize = 14;        
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular eccentricity'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6,'plot_omnibus',1,'omnibus_text_y_gap_coeff',1);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))    
    
    ax3 = axes('Units','pixels','Position',[ax2.Position(1)+230,ax1.Position(2),ax1.Position(3),ax1.Position(4)]);  
        gdata = simdata.fields;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,3),gdata(:,4));     
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XLim = [0.5 2.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-10 10]; 
        ax.YTick = -10:2:10;
        ax.FontSize = 14;
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular total fields'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6,'plot_omnibus',1,'omnibus_text_y_gap_coeff',1);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))    
    
            
exportgraphics(gcf,['C:\Users\F004KS7\Downloads\fig2.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  
                    
        
    
    
    
    
    
    
    
    

    fig_clust = figure('visible','on','Units','pixels','Position',[50, 50, 1850, 920]); % open/create figure
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    dot_size = [22 30 22 30]*2;
    dot_sigma = 0.12;
    
    ax = axes('Units','pixels','Position',[100,450,250,250]);  
        gdata = simdata.gscore;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,1),gdata(:,2),gdata(:,3),gdata(:,4));     
        cols = [rgb('Black');rgb('Black');rgb('Blue');rgb('Blue')];
        marks = {'o','s','o','s'};
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-1 1]; 
        ax.YTick = -1:0.25:1;
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular grid score'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))

    ax = axes('Units','pixels','Position',[500,450,250,250]);    
        gdata = simdata.eccentricity;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,1),gdata(:,2),gdata(:,3),gdata(:,4));     
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-1 1]; 
        ax.YTick = -1:0.25:1;
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular eccentricity'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])    
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
    
%     ax = axes('Units','pixels','Position',[550,650,200,180]);    
%         [ds,gs] = vectorDATAGROUP([],simdata.orientation(:,1),simdata.orientation(:,2),simdata.orientation(:,3),simdata.orientation(:,4),simdata.orientation(:,5));     
%         cols = winter(5);
%         cmat = mat2cell(cols,ones(size(cols,1),1));
%         meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',18,'dotalpha',0.5,'dotsigma',0.1,'dotcolor',cmat); hold on;
% 
%         ax = gca;
%         ax.XLim = [.5 max(gs(:))+.5];
%         ax.XTick = 1:max(gs(:));
% %         ax.YLim(1) = 0; 
% %         ax.YTick = 0:0.25:1.5;
%         pnames = {'B','P-P','P-S','S-P','S-S'};
%         ax.XTickLabel = pnames;
%         ax.XTickLabelRotation = 25;
%         box off
%         ylabel(sprintf('Orientation'))  
%         line(ax.XLim,[0 0],'Color',[.5 .5 .5])    
        
    ax = axes('Units','pixels','Position',[900,450,250,250]);
        gdata = simdata.fields;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,1),gdata(:,2),gdata(:,3),gdata(:,4));     
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat,'dotsmarker',marks); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XTick = 1:max(gs(:));
        ax.YLim = [-10 10]; 
        ax.YTick = -10:2:10;
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Horizontal - Irregular total fields'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])    
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-2.6);        
        set(ax,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
        
    ax = axes('Units','pixels','Position',[1050,650,200,180]);
        gdata = simdata.field_elongation;
        gdata = gdata(:,2:end) - gdata(:,1);
        [ds,gs] = vectorDATAGROUP([],gdata(:,1),gdata(:,2),gdata(:,3),gdata(:,4));     
        meanplot(ds,gs,'linecolor',{'k'},'meansize',6,'dots',1,'dotsize',dot_size,'dotalpha',0.5,'dotsigma',dot_sigma,'dotcolor',cmat); hold on;

        ax = gca;
        ax.XLim = [.5 max(gs(:))+.5];
        ax.XTick = 1:max(gs(:));
%         ax.YLim = [-10 10]; 
%         ax.YTick = -10:2:10;
        pnames = {'P-P','P-S','S-P','S-S'};
        ax.XTickLabel = pnames;
        ax.XTickLabelRotation = 25;
        box off
        ylabel(sprintf('Field elongation'))  
        line(ax.XLim,[0 0],'Color',[.5 .5 .5])    
        [result1,result2] = plotsigbrackets(ds,gs,'bracket_text_y_gap_coeff',-1.8);        
        
exportgraphics(gcf,['C:\Users\F004KS7\Downloads\fig2.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
keyboard    
    
   




    
    keyboard
    
figure
    % planar data spikes
    subplot(2,4,1)
        plot(pos_flat_curve(:,1),pos_flat_curve(:,2),'k'); hold on;
        plot(pspk(:,1),pspk(:,2),'r.','MarkerSize',10);     
        daspect([1 1 1])    
        title(sprintf('Planar data (%d spikes)',numel(pspk(:,1))))    
    
    % planar data map (planar projection)
    subplot(2,4,2)
        imagesc(rmap1);
        axis tight on vis3d; 
        rotate3d on;
        view(0,90);    
        daspect([1 1 1]); 
        grid on
        colormap(gca,inferno);
        title('Planar data planar projection')    
    
    % planar data map (surficial projection)
    subplot(2,4,3)
        imagesc(rmap2);
        axis tight on vis3d; 
        rotate3d on;
        view(0,90);    
        daspect([1 1 1]); 
        grid on
        colormap(gca,inferno);
        title('Planar data surficial projection')  
        
        
    % surficial data spikes
    subplot(2,4,5)
        plot(pos_flat_curve(:,1),pos_flat_curve(:,2),'k'); hold on;
        plot(sspk(:,1),sspk(:,2),'r.','MarkerSize',10);     
        daspect([1 1 1])  
        title('Surficial data')    
    
    % surficial data map (planar projection)        
    subplot(2,4,6)
        imagesc(rmap3);    
        axis tight on vis3d; 
        rotate3d on;
        view(0,90);    
        daspect([1 1 1]); 
        grid on
        colormap(gca,inferno);
        title('Surficial data planar projection')

    % surficial data map (planar projection)
    subplot(2,4,7)
        imagesc(rmap4);
        axis tight on vis3d; 
        rotate3d on;
        view(0,90);    
        daspect([1 1 1]); 
        grid on
        colormap(gca,inferno);
        title('Surficial data surficial projection')

return

















