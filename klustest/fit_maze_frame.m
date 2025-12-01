function [lx,ly,lz] = fit_maze_frame(pox,poy,poz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes some position data and allows the user to generate a lattice which can be used in later figures etc
%   It can also be used to generate a plain cuboid for open field environments (with additional inputs)
%   [lx,ly,lz] = makeLATTICE(pox,poy,poz,lsize,lcomp,latrot)
%
%%%%%%%% Inputs
%   pox = position x data
%   poy = position y data
%   poz = position z data
%   pratio = pixel ratio (pixels per metre), should be 1000 for 3D reconstructed data
%   lsize = environment size in mm [x y z]
%   lcomp = the number of compartments in the maze, lengthways (6 for lattice maze, 1 for open field)
%   lbase = the height of the base of the maze (should be 0 if the cameras are calibrated with a checkerboard placed on the same stools etc).
%
%%%%%%%% Outputs
%   lx = lattice x data
%   ly = lattice y data
%   lz = lattice z data
%
%%%%%%%% Comments
%   26/03/17 created from makeLATTICE
%   26/03/17 added more parameters so it can be used for other environments
%   06/04/18 adapted for diagonal lattice
%   06/04/18 added figure saving
%   06/04/18 added switch for different maze orientations
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
    print_now = 1; % save a figure at the end or not

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate a lattice
    % vectors for lattice generation
    x = linspace(0,3000,2);
    y = linspace(0,1500,2);
    z = linspace(0,500,2);

    [X1,Y1,Z1] = meshgrid(x([1 end]),y,z);
    X1 = permute(X1,[2 1 3]); 
    Y1 = permute(Y1,[2 1 3]); 
    Z1 = permute(Z1,[2 1 3]);
    X1(end+1,:,:) = NaN; 
    Y1(end+1,:,:) = NaN; 
    Z1(end+1,:,:) = NaN;

    [X2,Y2,Z2] = meshgrid(x,y([1 end]),z);
    X2(end+1,:,:) = NaN; 
    Y2(end+1,:,:) = NaN; 
    Z2(end+1,:,:) = NaN;

    [X3,Y3,Z3] = meshgrid(x,y,z([1 end]));
    X3 = permute(X3,[3 1 2]); 
    Y3 = permute(Y3,[3 1 2]); 
    Z3 = permute(Z3,[3 1 2]);
    X3(end+1,:,:) = NaN; 
    Y3(end+1,:,:) = NaN; 
    Z3(end+1,:,:) = NaN;

    lx = [X1(:); X2(:); X3(:)];
    ly = [Y1(:); Y2(:); Y3(:)];
    lz = [Z1(:); Z2(:); Z3(:)];

    % centre on the position data
    lx = lx-nanmean(lx);
    lx = lx + nanmean([min(pox) max(pox)]);
    ly = ly-nanmean(ly);
    ly = ly + nanmean([min(poy) max(poy)]);    
    lz = lz-nanmean(lz);
    lz = lz - min(lz(:)); % place the maze frame so its base it at 0mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overlay lattice on data so user can position it
    fig_check = figure('visible','on','Position',[100, 100, 1800, 800]);
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    sval = 120;
    axlims = [nanmin(lx)-sval nanmax(lx)+sval nanmin(ly)-sval nanmax(ly)+sval nanmin(lz)-sval nanmax(lz)+sval];
    pos_lwidth = 1;
    lat_lwidth = 1;
    pos_color = [0 0 1];
    lat_color = [1 0 0];

%% XY view of data
    ax_xy = axes('Units','pixels','Position',[50 200 500 500]);
    plot3(pox,poy,poz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
    hold on
    l1 = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
    axis(axlims);
    daspect([1 1 1])
    view(0,90)
    xlabel('X mm')
    ylabel('Y mm')
    zlabel('Z mm')
    grid on
    ax_xy.ZDir = 'reverse';    
    title('X vs Y','FontSize',12,'FontWeight','normal')
    
%% XZ view of data
    ax_xz = axes('Units','pixels','Position',[590 200 500 500]);
    plot3(pox,poy,poz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
    hold on
    l2 = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
    axis(axlims);
    daspect([1 1 1])
    view(0,0)
    xlabel('X mm')
    ylabel('Y mm')
    zlabel('Z mm')
    grid on   
    title('X vs Z','FontSize',12,'FontWeight','normal')
    
%% YZ view of data
    ax_yz = axes('Units','pixels','Position',[1130 200 500 500]);
    plot3(pox,poy,poz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
    hold on
    l3 = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
    axis(axlims);
    daspect([1 1 1])
    view(90,0)
    xlabel('X mm')
    ylabel('Y mm')
    zlabel('Z mm')
    grid on  
    title('Y vs Z','FontSize',12,'FontWeight','normal')

%% Sliders to control lattice position    
    % Slider to control lattice X
    sld1 = uicontrol(fig_check,'Style', 'slider','Min',ax_xz.XLim(1),'Max',ax_xz.XLim(2),'Value',nanmean(lx),'Position',[590 90 500 20],'String','X','Callback',@lat_shift); % @lat_shift is the callback function which can be seen below
    set(sld1,'SliderStep',[1/360,10/360]);    
   
    % Slider to control lattice Y
    sld2 = uicontrol(fig_check,'Style', 'slider','Min',ax_yz.XLim(1),'Max',ax_yz.XLim(2),'Value',nanmean(ly),'Position',[1130 90 500 20],'String','Y','Callback',@lat_shift); % @lat_shift is the callback function which can be seen below
    set(sld2,'SliderStep',[1/360,10/360]);      
    
    % Slider to control lattice Z
    sld3 = uicontrol(fig_check,'Style', 'slider','Min',ax_yz.XLim(1),'Max',ax_yz.XLim(2),'Value',nanmean(lz),'Position',[1680 200 20 500],'String','Z','Callback',@lat_shift); % @lat_shift is the callback function which can be seen below
    set(sld3,'SliderStep',[1/360,10/360]);      
    
%% Complete process      
    % Create a button which will allow the user to continue once they have positioned the lattice
    uicontrol(fig_check,'Style','pushbutton','String','Finished','Position',[700 30 150 30],'Callback','uiresume(gcbf)');  
    uiwait(gcf); % at this point start waiting for the button to be pressed    
    
    % when the button is pressed get the amount of X-shift
    xshift = get(sld1,'value');
    yshift = get(sld2,'value');
    zshift = get(sld3,'value');

    % add this to lattice  shift
    lx = lx-nanmean(lx);
    lx = lx + xshift;
    ly = ly-nanmean(ly);
    ly = ly + yshift;        
    lz = lz-nanmean(lz);
    lz = lz + zshift;      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to control lattice XY position
    function lat_shift(source,~)
        % get the amount of X-shift
        xshift = get(sld1,'value');
        yshift = get(sld2,'value');
        zshift = get(sld3,'value');
               
        % add any shift
        lx = lx-nanmean(lx);
        lx = lx + xshift;
        ly = ly-nanmean(ly);
        ly = ly + yshift;        
        lz = lz-nanmean(lz);
        lz = lz + zshift;    
        
        % update XY plot
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_xy);    
        delete(l1);
        hold on
        l1 = line(lx,ly,lz,'Color',[.5 .5 .5 .5],'LineWidth',lat_lwidth);
        view(0,90)
    
        % update XZ plot
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_xz);         
        delete(l2);
        hold on
        l2 = line(lx,ly,lz,'Color',[.5 .5 .5 .5],'LineWidth',lat_lwidth);    
        view(0,0)
    
        % update YZ plot
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_yz);          
        delete(l3);
        hold on
        l3 = line(lx,ly,lz,'Color',[.5 .5 .5 .5],'LineWidth',lat_lwidth);   
        view(90,0)
    end

    % save figure
    if print_now
        [~,~,~] = mkdir([pwd '\klustest\']);
        print(fig_check,'-dpng','-r250',[pwd '\klustest\Lattice_frame_' datestr(now,30) '.png'])
    end    
    close(fig_check);

end



