function [led_dat,rat_dat,ax_dat] = GIT_rat_calib(rat_name)



switch rat_name
    case {'generic_hpc'}
        % estimate the position of the rat's head nodes (ears and nose) based on what we know about the implant
        node = {'nose';'r_ear';'l_ear'};
        rat_dat = table(node); % prepare table  
        L = [0,0,0; -42,-21,0; -42,21,0]; % r,g,b 'idealized' coordinates   
        L = L - mean(L,1,'omitnan'); % center on their COM        
        rat_dat.coord(:,:) = L;        
        
        % get the positions of the LEDs (idealized coordinates based on the LED tracking frame)
        node = {'r';'g';'b'};
        yr = L + [-30 0 30]; % move backward [X _ _] mms, this is the distance of the implant behind the nose and up [_ _ Z] mms, this is the height of the LEDs above the head plane
        led_dat = table(node); % prepare table  
        led_dat.coord(:,:) = yr;

        % specify the axes of the rat's head
        node = {'x';'y';'z'};
        ax_dat = table(node); % prepare table          
        ax_dat.coord(:,:) = [1 0 0; 0 -1 0; 0 0 -1];      
    
    case {'generic_mec'}
        % estimate the position of the rat's head nodes (ears and nose) based on what we know about the implant
        node = {'nose';'r_ear';'l_ear'};
        rat_dat = table(node); % prepare table  
        L = [0,0,0; -42,-21,0; -42,21,0]; % r,g,b 'idealized' coordinates   
        L = L - mean(L,1,'omitnan'); % center on their COM        
        rat_dat.coord(:,:) = L;        
        
        % get the positions of the LEDs (idealized coordinates based on the LED tracking frame)
        node = {'r';'g';'b'};
        yr = L + [-30 0 30]; % move backward [X _ _] mms, this is the distance of the implant behind the nose and up [_ _ Z] mms, this is the height of the LEDs above the head plane
        yr = ( roty(8)*yr' )'; % rotate 8 degrees around the y-axis (pitch up), this is because the mEC implants are angled forward 8-10 degrees        
        led_dat = table(node); % prepare table  
        led_dat.coord(:,:) = yr;

        % specify the axes of the rat's head
        node = {'x';'y';'z'};
        ax_dat = table(node); % prepare table          
        ax_dat.coord(:,:) = [1 0 0; 0 -1 0; 0 0 -1]; 

    case {'RG7_MEC_TD'}
        led_color = {'r';'g';'b'};
        mm_to_left_ear = [55; 65; 60];
        mm_to_right_ear = [55; 58; 65];
        mm_to_nose = [70; 100; 100];
        tracker_calib = table(led_color,mm_to_left_ear,mm_to_right_ear,mm_to_nose);
        
        % get the positions of the LEDs (idealized coordinates based on the LED tracking frame)
        node = {'r';'g';'b'};
        led_coords = table(node); % prepare table  
        L = [0,0,0; -42,-21,0; -42,21,0]; % r,g,b 'idealized' coordinates
        L = L - mean(L,1,'omitnan'); % center on their COM
        led_coords.coord(:,:) = L;

        % estimate the position of the rat's head nodes (ears and nose) based on what we know about the implant
        node = {'l_ear';'r_ear';'nose'};
        rat_calib = table(node); % prepare table  
        yr = L([3 2 1],:) + [30 0 -30]; % move forward [X _ _] mms, this is the distance of the implant behind the nose and down [_ _ Z] mms, this is the height of the LEDs above the head plane
        yr = ( roty(8)*yr' )'; % rotate 8 degrees around the y-axis, this is because the mEC implants are angled forward 8-10 degrees
        rat_calib.coord(:,:) = yr;

        % specify the axes of the rat's head
        node = {'x';'y';'z'};
        rat_axes = table(node); % prepare table          
        rat_axes.coord(:,:) = [1 0 0; 0 -1 0; 0 0 -1];
        
    case {'RG8_MEC_AX'}
        led_color = {'r';'g';'b'};
        mm_to_left_ear = [55; 65; 55];
        mm_to_right_ear = [52; 55; 58];
        mm_to_nose = [75; 95; 102];
        tracker_calib = table(led_color,mm_to_left_ear,mm_to_right_ear,mm_to_nose);
        
        % get the positions of the LEDs (idealized coordinates based on the LED tracking frame)
        node = {'r';'g';'b'};
        led_coords = table(node); % prepare table  
        L = [0,0,0; -42,-21,0; -42,21,0]; % r,g,b 'idealized' coordinates
        L = L - mean(L,1,'omitnan'); % center on their COM
        led_coords.coord(:,:) = L;

        % estimate the position of the rat's head nodes (ears and nose) based on what we know about the implant
        node = {'l_ear';'r_ear';'nose'};
        rat_calib = table(node); % prepare table  
        yr = L([3 2 1],:) + [30 0 -30]; % move forward [X _ _] mms, this is the distance of the implant behind the nose and down [_ _ Z] mms, this is the height of the LEDs above the head plane
        yr = ( roty(8)*yr' )'; % rotate 8 degrees around the y-axis, this is because the mEC implants are angled forward 8-10 degrees
        rat_calib.coord(:,:) = yr;      
        
        % specify the axes of the rat's head
        node = {'x';'y';'z'};
        rat_axes = table(node); % prepare table          
        rat_axes.coord(:,:) = [1 0 0; 0 -1 0; 0 0 -1];
        
    case {'RG6_HPC_TD'}
        led_color = {'r';'g';'b'};
        mm_to_left_ear = [64; 74; 60];
        mm_to_right_ear = [64; 64; 68];
        mm_to_nose = [70; 100; 95];
        tracker_calib = table(led_color,mm_to_left_ear,mm_to_right_ear,mm_to_nose);
        
        % get the positions of the LEDs (idealized coordinates based on the LED tracking frame)
        node = {'r';'g';'b'};
        led_coords = table(node); % prepare table  
        L = [0,0,0; -42,-21,0; -42,21,0]; % r,g,b 'idealized' coordinates
        L = L - mean(L,1,'omitnan'); % center on their COM
        led_coords.coord(:,:) = L;

        % estimate the position of the rat's head nodes (ears and nose) based on what we know about the implant
        node = {'l_ear';'r_ear';'nose'};
        rat_calib = table(node); % prepare table  
        yr = L([3 2 1],:) + [30 0 -30]; % move forward [X _ _] mms, this is the distance of the implant behind the nose and down [_ _ Z] mms, this is the height of the LEDs above the head plane
        rat_calib.coord(:,:) = yr;
        
        % specify the axes of the rat's head
        node = {'x';'y';'z'};
        rat_axes = table(node); % prepare table          
        rat_axes.coord(:,:) = [1 0 0; 0 -1 0; 0 0 -1];
        
    case {'RG5_MEC_AX'}
        led_color = {'r';'g';'b'};
        mm_to_left_ear = [55; 72; 58];
        mm_to_right_ear = [55; 60; 63];
        mm_to_nose = [66; 99; 93];
        tracker_calib = table(led_color,mm_to_left_ear,mm_to_right_ear,mm_to_nose);
        
        % get the positions of the LEDs (idealized coordinates based on the LED tracking frame)
        node = {'r';'g';'b'};
        led_coords = table(node); % prepare table  
        L = [0,0,0; -42,-21,0; -42,21,0]; % r,g,b 'idealized' coordinates
        L = L - mean(L,1,'omitnan'); % center on their COM
        led_coords.coord(:,:) = L;

        % estimate the position of the rat's head nodes (ears and nose) based on what we know about the implant
        node = {'l_ear';'r_ear';'nose'};
        rat_calib = table(node); % prepare table  
        yr = L([3 2 1],:) + [30 0 -30]; % move forward [X _ _] mms, this is the distance of the implant behind the nose and down [_ _ Z] mms, this is the height of the LEDs above the head plane
        yr = ( roty(8)*yr' )'; % rotate 8 degrees around the y-axis, this is because the mEC implants are angled forward 8-10 degrees
        rat_calib.coord(:,:) = yr;
        
        % specify the axes of the rat's head
        node = {'x';'y';'z'};
        rat_axes = table(node); % prepare table          
        rat_axes.coord(:,:) = [1 0 0; 0 -1 0; 0 0 -1];
        
    otherwise
        error('Unknown rat... exiting')
        
end

% x = led_coords; % r, g, b
% d = table2array(tracker_calib(:,2:4));
% node = {'l_ear';'r_ear';'nose'};
% rat_calib = table(node);
% for k = 1:3
%     f = @(i)[d(1,k)^2 - (i(1)-x(1,1))^2 - (i(2)-x(1,2))^2 - (i(3)-x(1,3))^2;...
%     d(2,k)^2 - (i(1)-x(2,1))^2 - (i(2)-x(2,2))^2 - (i(3)-x(2,3))^2;...
%     d(3,k)^2 - (i(1)-x(3,1))^2 - (i(2)-x(3,2))^2 - (i(3)-x(3,3))^2;];
%     options = optimset('Algorithm',{'levenberg-marquardt',.1},'Display','off','TolX',1e-6);        
%     i0 = [0,0,-30];
%     p = lsqnonlin(f,i0,[-inf -inf -inf],[inf inf 0],options);
%     rat_calib.coord(k,:) = p;
% end   



















