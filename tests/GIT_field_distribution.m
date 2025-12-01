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
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY

% pnow = 1;
% % rmaps1 = cat(3,sdata.ratemap_surficial{sdata.frate(:,1)>active_cutoff & sdata.partn==pnow & sdata.repetition_score(:,1)>0.5});
% rmaps1 = cat(3,sdata.ratemap_planar{pidx & sdata.partn==pnow});
% % rmaps1 = cat(3,sdata.ratemap_planar{pidx & sdata.partn==pnow & sdata.repetition_score(:,1)>0.5});
% 
% % rmaps1 = {sdata.ratemap_surficial{pidx & sdata.partn==pnow}};
% % rmaps1 = cellfun(@(x) imresize(x,[48 116],'bilinear'),rmaps1,'UniformOutput',0);
% % rmaps1 = cat(3,rmaps1{:});
% 
% offs = NaN(size(rmaps1,3),1);
% cents = NaN(size(rmaps1,3),size(rmaps1,2));
% for ii = 1:size(rmaps1,3)
%     m = rmaps1(:,:,ii);
%     if isempty(m)
%         continue
%     end
%     s = sum(m,1,'omitnan');
%     [mval,maxidx] = max(s);
%     
%     xpoints = linspace(1,size(m,2),7);
% %     xpoints = xpoints(2:2:end);
%     xpoints = xpoints(4);
% 
%     kindx = knnsearch(xpoints(:),maxidx);
% 
%     offs(ii) = maxidx - xpoints(kindx);
%     cents(ii,:) = s;
% 
%     figure
%     subplot(2,2,1)
%     imagesc(m); hold on;
%     xpoints = linspace(1,size(m,2),7);
%     ax = gca;
%     plot([xpoints(1:2:end)',xpoints(1:2:end)'],ax.YLim,'w-',[xpoints(2:2:end)',xpoints(2:2:end)'],ax.YLim,'w--')
% 
%     subplot(2,2,2)
%     imagesc(s)  
%     
%     subplot(2,2,3)
%     plot(s); hold on;
%     plot(maxidx,s(maxidx),'ko');
%     ax = gca;
%     line([xpoints(kindx) xpoints(kindx)],ax.YLim,'Color','b')
%     title(sprintf('%.1f',offs(ii)))
%     keyboard
%     
%     
%     
% end
% [~,sidx] = sort(offs,'descend');
% 
% 
% 
% figure
% cmat = cents(sidx,:);
% cmat = cmat ./ max(cmat,[],2,'omitnan');
% % cmat = cmat ./ sum(cmat,2,'omitnan');
% imagesc(cmat); hold on;
% axis xy tight
% colormap('turbo')
% caxis([0 1])
% 
% xpoints = linspace(1,size(m,2),7);
% ax = gca;
% plot([xpoints(1:2:end)',xpoints(1:2:end)'],ax.YLim,'w-',[xpoints(2:2:end)',xpoints(2:2:end)'],ax.YLim,'w--')
% 
% xlabel('x-position')
% ylabel('cell')
% 
% % exportgraphics(fig_clust,fname,'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',150);  
% 
% 
% return
% 
% 
% 
% 
% 
% figure
% tiledlayout(8,8,'TileSpacing','compact','Padding','compact');
% for ii = 1:size(rmaps1,3)
%     nexttile
%     imagesc(rmaps1(:,:,sidx(ii))); hold on;
%     axis xy off tight
%     daspect([1 1 1])
%     colormap('turbo')
% 
%     xspac = size(med_map,2)/6;
%     xpoints = 0:xspac:size(med_map,2);
%     ax = gca;
%     plot([xpoints(1:2:end)',xpoints(1:2:end)'],ax.YLim,'k-',[xpoints(2:2:end)',xpoints(2:2:end)'],ax.YLim,'k--')
% end
% 
% 
% 
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
% 





pnow = 2;
rmaps1 = cat(3,sdata.ratemap_planar{pidx & sdata.partn==pnow});
med_map = median(rmaps1,3,'omitnan');

fdat = cat(1,sdata.planar_field_data{pidx & sdata.partn==pnow}); % [area, centroid x, centroid y, majaxis, minaxis, orientation, weight centroid x, weight centroid y, max intensity]

figure
subplot(2,2,1)
    imagesc(med_map); hold on;
    axis xy
    daspect([1 1 1])

    xspac = size(med_map,2)/6;
    xpoints = 0:xspac:size(med_map,2);
    ax = gca;
    plot([xpoints(1:2:end)',xpoints(1:2:end)'],ax.YLim,'k-',[xpoints(2:2:end)',xpoints(2:2:end)'],ax.YLim,'k--')

subplot(2,2,2)
    plot(fdat(:,7),fdat(:,8),'ko'); hold on;
    axis ij
    daspect([1 1 1])

    xspac = size(med_map,2)/6;
    xpoints = 0:xspac:size(med_map,2);
    ax = gca;
    plot([xpoints(1:2:end)',xpoints(1:2:end)'],ax.YLim,'k-',[xpoints(2:2:end)',xpoints(2:2:end)'],ax.YLim,'k--')

subplot(2,2,4)
    h = histcounts2(fdat(:,8),fdat(:,7),0.5:1:size(med_map,1)+0.5,0.5:1:size(med_map,2)+0.5);
    h = imgaussfilt(h,2);
    imagesc(h); hold on;
    axis ij
    daspect([1 1 1])

    xspac = size(med_map,2)/6;
    xpoints = 0:xspac:size(med_map,2);
    ax = gca;
    plot([xpoints(1:2:end)',xpoints(1:2:end)'],ax.YLim,'k-',[xpoints(2:2:end)',xpoints(2:2:end)'],ax.YLim,'k--')
    
    return


figure
tiledlayout(10,10,'TileSpacing','compact','Padding','compact');
for ii = 1:100
    nexttile
    imagesc(rmaps1(:,:,ii)); hold on;
    axis xy off tight
    daspect([1 1 1])
    colormap('turbo')

    xspac = size(med_map,2)/6;
    xpoints = 0:xspac:size(med_map,2);
    ax = gca;
    plot([xpoints(1:2:end)',xpoints(1:2:end)'],ax.YLim,'k-',[xpoints(2:2:end)',xpoints(2:2:end)'],ax.YLim,'k--')
end




figure
plot(sum(med_map,1))




m = squeeze(sum(rmaps1,1,'omitnan'))';
% z = m ./ sum(m,2,'omitnan');
z = zscore(m,[],2);

figure
subplot(2,1,1)
plot(z')

subplot(2,1,2)
plot(mean(z,1))


figure
idx = 1;
subplot(3,1,1)
imagesc(rmaps1(:,:,idx))

subplot(3,1,2)
plot(sum(rmaps1(:,:,idx),1,'omitnan'))

subplot(3,1,3)
plot(m(:,idx))
























































