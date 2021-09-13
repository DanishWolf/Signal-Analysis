% Generate Time Vector 
t3 = linspace(-2, 2, 999 + 1); t3(end) = [];
% Generate signal vector s3(t) in time domain
s3t = (rectangularPulse(t3./6) + 6.*exp(-9.*t3).*heaviside(t3).*cos(2.*pi.*33.*t3))...
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
mag_S3k= plot(t3, fftshift(abs(S3k))/fs);
title('Magnitude of Signal in Frenquency Domain')
xlabel('Time')
ylabel('Amplitude')
% store the phase plot as 'phase_S3k'
figure(2)
phase_S3k= plot(t3, fftshift(angle(S3k))/fs);
title('Phase of Signal in Freqency Domain')
xlabel('Time')
ylabel('Amplitude')