%% Sample solution for EGB242 - Comp Lab 7
% Confirmed 2015

clear all
close all
clc

%% Part 1 - Time Domain Analysis
% a) Read in noisy sound file
[x, fs] = audioread('NoisyLaCampanella.wav');

% b) Play sound
% sound(x,fs)

% c) Create time vector
Ls = length(x); % length of sound signal
Ts = 1/fs;
t = 0:Ts:Ls*Ts;
t = t(1:end-1);

% d) Plot time signal
figure
plot(t,x)
xlabel('Time');ylabel('Amplitude')
title('Noisy signal in time domain')

% e) Nothing much can be gained by looking at the time signal

%% Part 2 - Frequency Domain Analysis
% a) Fourier Transform
% Here we are not applying fftshift or scaling, since we will be further
% processing the signal. Only apply fftshift or scaling when plotting
X = fft(x); 

% b) Frequency vector
f = linspace(-fs/2,fs/2,Ls+1);
f = f(1:end-1);

% c) Plot magnitude spectrum
figure
plot(f,abs(fftshift(X))/fs) % here we apply fftshift and scale amplitude
xlabel('Frequency');ylabel('Magnitude')
title('Noisy signal in frequency domain')

% d) The clean signal has bandwidth of about 6kHz. We see 3 major noise
% peaks beginning at 9kHz, 15kHz, 20kHz respectively. Each noise peak has
% bandwidth 2kHz.

%% Part 3 - LTI Systems, Impulse Responses
% Load data file contains impulse responses. There are 4 filters, named
% hn1, hn2, hn3, hn4
load LTI_ImpRes.mat

% Investigate frequency response of each filter. Remember to create
% frequency vector for each filter. We are looking at frequency
% plots, so fftshift is applied.  Scaling is not applied to these.
% Pay attention to the magnitude response of these filters - when we
% pass the noisy signal through them, how would we expect the magnitudes
% of the output frequency to differ?

% Filter 1
fh1 = linspace(-fs/2,fs/2,length(hn1)+1);
fh1 = fh1(1:end-1);
Hf1 = fftshift(fft(hn1));

% Filter 2
fh2= linspace(-fs/2,fs/2,length(hn2)+1);
fh2 = fh2(1:end-1);
Hf2 = fftshift(fft(hn2));

% Filter 3
fh3 = linspace(-fs/2,fs/2,length(hn3)+1);
fh3 = fh3(1:end-1);
Hf3 = fftshift(fft(hn3));

% Filter 4
fh4 = linspace(-fs/2,fs/2,length(hn4)+1);
fh4 = fh4(1:end-1);
Hf4 = fftshift(fft(hn4));

% Plot frequency response of filters, on top of noisy signal spectrum

figure
subplot(2,2,1)
plot(fh1,abs(Hf1))
xlabel('Frequency');ylabel('Magnitude')
title('Spectrum of filter 1')

subplot(2,2,2)
plot(fh2,abs(Hf2))
xlabel('Frequency');ylabel('Magnitude')
title('Spectrum of filter 2')

subplot(2,2,3)
plot(fh3,abs(Hf3))
xlabel('Frequency');ylabel('Magnitude')
title('Spectrum of filter 3')

subplot(2,2,4)
plot(fh4,abs(Hf4))
xlabel('Frequency');ylabel('Magnitude')
title('Spectrum of filter 4')

% >> We can make the following observation:
% Filter 1: low pass filter, stop band frequency @ 15kHz
% Filter 2: band stop filter, stop band frequency @ 11kHz, 13kHz
% Filter 3: band pass filter, stop band frequency @ 6kHz, 20kHz
% Filter 4: high pass filter, stop band frequency @ 11kHz
% 
% To remove all noise peaks, we can apply Filter 1 to clean the 15kHz and
% 20kHz noise, then apply Filter 2 to clean the 9kHz noise.


%% Part 4 - Convolution
%  The output of a LTI system is the convolution between the input signal
%  and the impulse response of the system. Convolution in MATLAB is 
%  performed using "conv(h,x)" command.
%
%  NOTE: the number of points in the output of conv(h,x) is
%  length(h)+length(x)-1. Therefore we need to truncate the output to
%  length(x) points. You'll learn about this in ENB342.

figure
% plot original spectrum
subplot(3,1,1)
plot(f,abs(fftshift(X))/fs) 
xlabel('Frequency');ylabel('Magnitude')
title('Before filtering')

% Apply Filter 1
x1 = conv(x,hn1); 
x1 = x1(1:Ls); % truncate
% Compute spectrum
X1 = fftshift(fft(x1))/fs;

% plot spectrum after filter 1
subplot(3,1,2)
plot(f,abs(X1))
xlabel('Frequency');ylabel('Magnitude')
title('Spectrum after applying Filter 1')

% Apply Filter 2
x2 = conv(x1,hn2);
x2 = x2(1:Ls);
% Compute spectrum
X2 = fftshift(fft(x2))/fs;

% plot spectrum after filter 1
subplot(3,1,3)
plot(f,abs(X2))
xlabel('Frequency');ylabel('Magnitude')
title('Spectrum after applying Filter 2')

% Finally, play the cleaned up sound
% sound(x2,fs)