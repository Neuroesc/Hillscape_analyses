function [lx,ly,lz] = fit_hill_frame_v2(pdata,skipman,mname,fname)
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
    pos = pdata.pos;
    pox = pos.pox(:,1); % pos x for all data, cm
    poy = pos.poy(:,1); % pos y for all data, cm
    session_times = pdata.session_times;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate a lattice
    % vectors for lattice generation
    x = linspace(0,3000,2);
    y = linspace(0,1500,2);
    z = linspace(0,750,2);

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
    fig_check = figure('visible','on','Position',[100, 100, 1400, 800]);
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    sval = 120;
    axlims = [nanmin(lx)-sval nanmax(lx)+sval nanmin(ly)-sval nanmax(ly)+sval nanmin(lz)-sval nanmax(lz)+sval];
    pos_lwidth = 1;
    lat_lwidth = 1;
    pos_color = [0 0 0];
    lat_color = [1 0 0];
    xbuff = 40;
    ybuff = 200;

%% >>>>>>>>>> arena1 data
    %% XY view of arena1 data
    ax_xy1 = axes('Units','pixels','Position',[50 500 400 250],'Clipping','off');
        tindx = pos.pot > session_times(1,1) & pos.pot < session_times(1,2); % index for position data in this part
        ppox = double( pos.pox( tindx ) ); % pos x for this part, cm
        ppoy = double( pos.poy( tindx ) ); % pos y for this part, cm
        ppoz = -double( pos.poz( tindx ) ); % pos z for this part, cm 
        part_now = 'arena 1';
        
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
        hold on
        l1_a = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(0,90)
        grid on
        axis off tight
        ax_xy1.ZDir = 'reverse';    
        fpos = 1.2;
        text(0.5,fpos,sprintf('%s - X vs Y',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
    
    %% XZ view of arena1 data
    ax_xz1 = axes('Units','pixels','Position',[ax_xy1.Position(1) ax_xy1.Position(2)-ybuff ax_xy1.Position(3) ax_xy1.Position(4)],'Clipping','off');
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
        hold on
        l2_a = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(0,0)
        grid on   
        axis off tight    
        text(0.5,fpos,sprintf('%s - X vs Z',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
        
    %% YZ view of arena1 data
    ax_yz1 = axes('Units','pixels','Position',[ax_xz1.Position(1) ax_xz1.Position(2)-ybuff ax_xz1.Position(3) ax_xz1.Position(4)],'Clipping','off');
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
        hold on
        l3_a = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(90,0)
        grid on  
        axis off tight    
        text(0.5,fpos,sprintf('%s - Y vs Z',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
    
%% >>>>>>>>>> hills data
    %% XY view of hills data
    ax_xy2 = axes('Units','pixels','Position',[ax_xy1.Position(1)+ax_xy1.Position(3)+xbuff ax_xy1.Position(2) ax_xy1.Position(3) ax_xy1.Position(4)],'Clipping','off');
        tindx = pos.pot > session_times(2,1) & pos.pot < session_times(2,2); % index for position data in this part
        ppox = double( pos.pox( tindx ) ); % pos x for this part, cm
        ppoy = double( pos.poy( tindx ) ); % pos y for this part, cm
        ppoz = -double( pos.poz( tindx ) ); % pos z for this part, cm 
        part_now = 'hills';   
        
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
        hold on
        l1_b = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(0,90)
        grid on
        axis off tight 
        ax_xy2.ZDir = 'reverse';    
        text(1,1.1,'X vs Y','Units','normalized','HorizontalAl','right','VerticalAl','top','FontSize',8)
        text(0.5,fpos,sprintf('%s - X vs Y',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
    
    %% XZ view of arena1 data
    ax_xz2 = axes('Units','pixels','Position',[ax_xy2.Position(1) ax_xz1.Position(2) ax_xy2.Position(3) ax_xy2.Position(4)],'Clipping','off');
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
        hold on
        l2_b = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(0,0)
        grid on   
        axis off tight     
        text(0.5,fpos,sprintf('%s - X vs Z',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
        
    %% YZ view of arena1 data
    ax_yz2 = axes('Units','pixels','Position',[ax_xz2.Position(1) ax_yz1.Position(2) ax_xz2.Position(3) ax_xz2.Position(4)],'Clipping','off');
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
        hold on
        l3_b = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(90,0)
        grid on  
        axis off tight     
        text(0.5,fpos,sprintf('%s - Y vs Z',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
        
%% >>>>>>>>>> arena2 data
    %% XY view of arena2 data
    ax_xy3 = axes('Units','pixels','Position',[ax_xy2.Position(1)+ax_xy2.Position(3)+xbuff ax_xy2.Position(2) ax_xy2.Position(3) ax_xy2.Position(4)]);
        tindx = pos.pot > session_times(3,1) & pos.pot < session_times(3,2); % index for position data in this part
        ppox = double( pos.pox( tindx ) ); % pos x for this part, cm
        ppoy = double( pos.poy( tindx ) ); % pos y for this part, cm
        ppoz = -double( pos.poz( tindx ) ); % pos z for this part, cm 
        part_now = 'arena 2';  
        
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth,'Clipping','off'); % plot position data for user to see
        hold on
        l1_c = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(0,90)
        grid on
        axis off tight 
        ax_xy3.ZDir = 'reverse';    
        text(0.5,fpos,sprintf('%s - X vs Y',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
    
    %% XZ view of arena1 data
    ax_xz3 = axes('Units','pixels','Position',[ax_xy3.Position(1) ax_xz1.Position(2) ax_xy3.Position(3) ax_xy3.Position(4)],'Clipping','off');
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
        hold on
        l2_c = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(0,0)
        grid on   
        axis off tight     
        text(0.5,fpos,sprintf('%s - X vs Z',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
        
    %% YZ view of arena1 data
    ax_yz3 = axes('Units','pixels','Position',[ax_xz3.Position(1) ax_yz1.Position(2) ax_xz3.Position(3) ax_xz3.Position(4)],'Clipping','off');
        plot3(ppox,ppoy,ppoz,'Color',pos_color,'LineWidth',pos_lwidth); % plot position data for user to see
        hold on
        l3_c = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        axis(axlims);
        daspect([1 1 1])
        view(90,0)
        grid on  
        axis off tight     
        text(0.5,fpos,sprintf('%s - Y vs Z',part_now),'Units','normalized','HorizontalAl','center','VerticalAl','top','FontSize',12)
    
%% Sliders to control lattice position    
    % Slider to control lattice X
    sld1 = uicontrol(fig_check,'Style', 'slider','Min',ax_xz1.XLim(1),'Max',ax_xz1.XLim(2),'Value',nanmean(lx),'Position',[50 110 500 20],'String','X','Callback',@lat_shift); % @lat_shift is the callback function which can be seen below
    set(sld1,'SliderStep',[1/360,10/360]);    
   
    % Slider to control lattice Y
    sld2 = uicontrol(fig_check,'Style', 'slider','Min',ax_yz1.XLim(1),'Max',ax_yz1.XLim(2),'Value',nanmean(ly),'Position',[50 85 500 20],'String','Y','Callback',@lat_shift); % @lat_shift is the callback function which can be seen below
    set(sld2,'SliderStep',[1/360,10/360]);      
    
    % Slider to control lattice Z
    sld3 = uicontrol(fig_check,'Style', 'slider','Min',ax_yz1.XLim(1),'Max',ax_yz1.XLim(2),'Value',nanmean(lz),'Position',[50 60 500 20],'String','Z','Callback',@lat_shift); % @lat_shift is the callback function which can be seen below
    set(sld3,'SliderStep',[1/360,10/360]);      
    
%% Complete process      
    % Create a button which will allow the user to continue once they have positioned the lattice
    uicontrol(fig_check,'Style','pushbutton','String','Finished','Position',[50 20 150 30],'Callback','uiresume(gcbf)'); 
    if skipman
        uiresume(gcf); % just keep going
    else
        uiwait(gcf); % at this point start waiting for the button to be pressed    
    end
    
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
        
        % update XY plot (arena1)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_xy1);    
        delete(l1_a);
        hold on
        l1_a = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        view(0,90)
    
        % update XZ plot (arena1)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_xz1);         
        delete(l2_a);
        hold on
        l2_a = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);    
        view(0,0)
    
        % update YZ plot (arena1)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_yz1);          
        delete(l3_a);
        hold on
        l3_a = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);   
        view(90,0)
        
        % update XY plot (hills)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_xy2);    
        delete(l1_b);
        hold on
        l1_b = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        view(0,90)
    
        % update XZ plot (hills)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_xz2);         
        delete(l2_b);
        hold on
        l2_b = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);    
        view(0,0)
    
        % update YZ plot (hills)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_yz2);          
        delete(l3_b);
        hold on
        l3_b = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);   
        view(90,0)        
        
        % update XY plot (arena2)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_xy3);    
        delete(l1_c);
        hold on
        l1_c = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);
        view(0,90)
    
        % update XZ plot (arena2)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_xz3);         
        delete(l2_c);
        hold on
        l2_c = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);    
        view(0,0)
    
        % update YZ plot (arena2)
        set(0,'currentfigure',fig_check);
        set(fig_check,'currentaxes',ax_yz3);          
        delete(l3_c);
        hold on
        l3_c = line(lx,ly,lz,'Color',lat_color,'LineWidth',lat_lwidth);   
        view(90,0)         
        
    end

    % save data
    save(mname,'lx','ly','lz','-mat'); 

    % save figure
    if print_now
        print(fig_check,'-dpng','-r250',fname)
    end    
    close(fig_check);

    analysis_log({'fit_hill_frame_v2'},1,'version',{'v2.0.0'});


end



