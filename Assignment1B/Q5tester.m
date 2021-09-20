%% Question 5
% Generate Time Vector 
t3 = linspace(-2, 2, 999 + 1); t3(end) = [];
% Generate signal vector s3(t) in time domain
s3t = (rectangularPulse(t3./6) + 6.*exp(-9*t3).*heaviside(t3).*cos(2*pi*33*t3))...
    .* cos(2 .* pi .* 33 .* t3);
% Determine the sampling period
Ts = t3(2) - t3(1);
% Determine the Sampling Frequency
fs= 1/Ts;
% Generate frequency vector
k3 = linspace(-fs/2, fs/2, 999 + 1); k3(end) = [];
% Perform Discrete Fourier Transform of s3t
S3k= fft(s3t);

% store the magnitude plot as 'mag_S3k'
figure(1)
mag_S3k= plot(k3, fftshift(abs(S3k))/fs, 'b');
title('Magnitude of Signal in Frenquency Domain')
xlabel('Frequency')
ylabel('Magnitude')
% store the phase plot as 'phase_S3k'
figure(2)
phase_S3k= plot(k3, fftshift(angle(S3k))/fs, 'b');
title('Phase of Signal in Freqency Domain')
xlabel('Frequency')
ylabel('Phase')

%% Question 5 Don's Example
% Generate Time Vector
t3 = linspace(-2,2, 999 + 1);
t3(end)=[];
% Generate signal vector s3(t) in time domain
s3t = (( 2.*exp(-8*t3) .* heaviside(t3) + 8 .* rectangularPulse(-5/2 , 5/2 ,t3) .* cos(60*pi*t3)) .* cos(60*pi*t3)) ; % Check this function
% Determine the sampling period
Ts = t3(2) - t3(1) ;
% Determine the Sampling Frequency
fs= 1/Ts;
% Generate frequency vector
k3 =linspace(-fs/2 , fs/2, 999 + 1);
k3(end)=[];
% Perform Discrete Fourier Transform of s3t
S3k= fft(s3t) ;

magnitude = abs(S3k);
phase = angle(S3k);

% storing the magnitude as 'mag_S3k'
figure(3)
mag_S3k= plot(k3, fftshift(magnitude)/fs , 'blue');
xlabel('frequency')
ylabel('magnitude')
title('Plotting Magnitude');
% storing the phase as 'phase_S3k'
figure(4)
phase_S3k= plot(k3, fftshift(phase),'blue') ;
xlabel('frequency')
ylabel('phase')
title('Plotting Phase');