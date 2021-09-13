%% EGB242 Week 7 In Class Challenge
clear all, close all, clc

% NOTE: We have chosen to look at the signal first as this allows us to
% find the frequency range that the filters were made for.

%% 1) Load sound file and observe spectrum
% Load sound file
[audio,fs] = audioread('AudioNoisy.wav');
% sound(audio, fs)

% Find spectrum as it's easier to see what's happening in the frequency
% domain, especially when considering filters
Audio = ;
f = ; f(end) = [];

% Plot magnitude spectrum to look at what the filter is going to have to do
figure(1), clf, plot() % <-------------  Time in class to attempt
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('Magntiude spectrum of noisy sound file');


%% 2) Identify filters
% Apply impulse and observe spectrum for each filter, create impulse of 50 pts long
imp = ; % <------------- Time in class to attempt

% NOTE: We use an impulse of height 1 because the fourier transform of an 
% impulse of amplitude 1 is 1 for every frequency. Looking at the
% definition for the output of an LTI system, in the frequency domain it
% looks like this:
% Y(f) = H(f) * X(f)
% If we input an impulse then X(f) is just a multiply by 1 so the output
% becomes:
% Y(f) = H(f)
% and therefore, we can view the exact behaviour of the system, also called
% the impulse response of the system. This is the main reason the behaviour
% of an LTI system is called the impulse response, its literally the
% response of the system to an impulse and shows exactly what is happening
% at every frequency

% System 1
h1 = ; % <------- Time in class to attempt
H1 = ;
f1 = linspace(-fs/2, fs/2, length(H1)+1); f1(end) = [];

% System 2
h2 = ;
H2 = ;
f2 = linspace(-fs/2, fs/2, length(H2)+1); f2(end) = [];

% System 3
h3 = ;
H3 = ;
f3 = linspace(-fs/2, fs/2, length(H3)+1); f3(end) = [];

% System 4
h4 = ;
H4 = ;
f4 = linspace(-fs/2, fs/2, length(H4)+1); f4(end) = [];

% System 5
h5 = ;
H5 = ;
f5 = linspace(-fs/2, fs/2, length(H5)+1); f5(end) = [];

% System 6
h6 = ;
H6 = ;
f6 = linspace(-fs/2, fs/2, length(H6)+1); f6(end) = [];

% Plot spectrums on one figure for easy observation
% Note: dirac delta as the input is just a pulse so we do not have to
% divide by anything as there is technically no sampling frequency

figure(2), clf
subplot(3,2,1), plot()
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('System 1')

subplot(3,2,2), plot()
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('System 2')

subplot(3,2,3), plot()
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('System 3')

subplot(3,2,4), plot()
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('System 4')

subplot(3,2,5), plot()
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('System 5')

subplot(3,2,6), plot()
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('System 6')

%% 3) Apply noisy sound file through necessary Systems
% Filter through system
clean1 = ;
clean2 = ;
clean3 = ;

% observe spectrum
Clean = ;
f = linspace(-fs/2, fs/2, length(Clean)+1); f(end) = [];

figure(3), clf
plot(f, fftshift(abs(Clean))/fs)
xlabel('Frequency (Hz)'), ylabel('Magnitude')
title('Magntiude spectrum of filtered sound file')

% sound(clean3, fs)