





idx = prctile(pdata.pos.pitch,[33 66]);
part_now = 'hills';


figure
subplot(1,2,1)
lw = 1.5;
c = cline(pdata.(part_now).pox,pdata.(part_now).poy,pdata.(part_now).poz,pdata.(part_now).poz); hold on;
set(c,'LineWidth',lw)
daspect([1 1 1])
colormap('turbo')
caxis([5 43])
% caxis([-pi/2 pi/2])
axis xy off
view(3) 

subplot(1,2,2)
lw = 1.5;
c = cline(pdata.(part_now).pox,pdata.(part_now).poy,pdata.(part_now).poz,pdata.(part_now).pit); hold on;
set(c,'LineWidth',lw)
daspect([1 1 1])
colormap('turbo')
% caxis([50 430])
caxis([-pi/2 pi/2])
axis xy off
view(3)  

exportgraphics(gcf,['C:\Users\F004KS7\Downloads\pos1.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  




figure
idx = prctile(pdata.arena1.pit,[50]);
% idx = 0;
% part_now = 'arena1';
% [x,y,z] = sph2cart(pdata.(part_now).yaw,pdata.(part_now).pit-idx,ones(size(pdata.(part_now).pit)));
% [Xs,Ys,Zs,F1] = projectEIGS([x(:) y(:) z(:)]);
% surf(Xs,Ys,Zs,F1,'EdgeColor','none');
% daspect([1 1 1])
% axis xy tight
% 
% figure
% imagesc(F1); hold on;
% axis xy
% xlabel('Yaw')
% ylabel('Pitch')
% ax = gca;
% line(ax.XLim,repmat(size(F1,1)/2,1,2),'Color','w','LineStyle','--')
% xp = linspace(1,size(F1,2),5);
% ax.XTick = xp;
% ax.XTickLabel = {sprintf('-%c',960),sprintf('-%c/2',960),'0',sprintf('%c/2',960),sprintf('%c',960)};
% yp = linspace(1,size(F1,1),5);
% ax.YTick = yp;
% ax.YTickLabel = {sprintf('-%c/2',960),sprintf('-%c/4',960),'0',sprintf('%c/4',960),sprintf('%c/2',960)};

figure
subplot(1,3,1)
idx = prctile(pdata.arena1.pit,[33 66]);
idxn = pit<idx(1);
pit = pdata.(part_now).pit;
poxn = pdata.(part_now).pox;
poyn = pdata.(part_now).poy;
pozn = pdata.(part_now).poz;
pitn = pit;
poxn(~idxn) = NaN;
poyn(~idxn) = NaN;
pozn(~idxn) = NaN;
pitn(~idxn) = NaN;
c = cline(poxn,poyn,pozn,pitn); hold on;
set(c,'LineWidth',lw)
caxis([-pi/2 pi/2])
daspect([1 1 1])
axis xy off
view(3)
rotate3d on

subplot(1,3,2)
idxn = pit>idx(1) & pit<idx(2);
pit = pdata.(part_now).pit;
poxn = pdata.(part_now).pox;
poyn = pdata.(part_now).poy;
pozn = pdata.(part_now).poz;
pitn = pit;
poxn(~idxn) = NaN;
poyn(~idxn) = NaN;
pozn(~idxn) = NaN;
pitn(~idxn) = NaN;
c = cline(poxn,poyn,pozn,pitn); hold on;
set(c,'LineWidth',lw)
caxis([-pi/2 pi/2])
daspect([1 1 1])
axis xy off
view(3)
rotate3d on

subplot(1,3,3)
idxn = pit>idx(2);
pit = pdata.(part_now).pit;
poxn = pdata.(part_now).pox;
poyn = pdata.(part_now).poy;
pozn = pdata.(part_now).poz;
pitn = pit;
poxn(~idxn) = NaN;
poyn(~idxn) = NaN;
pozn(~idxn) = NaN;
pitn(~idxn) = NaN;
c = cline(poxn,poyn,pozn,pitn); hold on;
set(c,'LineWidth',lw)
caxis([-pi/2 pi/2])
daspect([1 1 1])
axis xy off
view(3)
rotate3d on



exportgraphics(gcf,['C:\Users\F004KS7\Downloads\pos2.png'],'BackgroundColor',[1 1 1],'Colorspace','rgb','Resolution',250);  


































