    


close all

    % Create figure
    fig_now = figure('Units','pixels','Position',[50 50 1200 1000],'visible','on');
    set(gcf,'InvertHardCopy','off'); % gives the figure a grey background but means it will save white lines as white    
    set(gcf,'color','w'); % makes the background colour white

    xnow = 20;
    ynow = 580;
    sizx = 200;
    sizy = 180;
    xbuff = sizx+30;

    ax = axes('Units','pixels','Position',[xnow ynow sizx sizy]);
        PIT_plot_mazes(ax,1) % arena
            % camlight(ax,'right')

    xnow = xnow+xbuff;
    ax = axes('Units','pixels','Position',[xnow ynow sizx sizy]);
        PIT_plot_mazes(ax,2) % hills
            % camlight(ax,'left')

    xnow = xnow+xbuff;
    ax = axes('Units','pixels','Position',[xnow ynow sizx sizy]);
        PIT_plot_mazes(ax,3) % hemisphere
            % camlight(ax,'right')

    xnow = xnow+xbuff;
    ax = axes('Units','pixels','Position',[xnow ynow sizx sizy]);
        PIT_plot_mazes(ax,4) % saddle
            % camlight(ax,'right')






















