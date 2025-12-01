



fs = 1; % Sampling frequency (samples per second) 
 dt = 1/fs; % seconds per sample 
 StopTime = 3000; % seconds 
 t = (0:dt:3000)'; % seconds 
 F = 3/3000; % Sine wave frequency (hertz) 
 data = sin(2*pi*F*t);
 data = (data+1)*225;
 plot(t,data)

f = gradient(data(:));
s = rad2deg( atan( f ) );

cline(t,data,s);
colorbar;



























