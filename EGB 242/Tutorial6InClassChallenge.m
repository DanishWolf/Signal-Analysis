%% EGB242 In Class Challenge 6 - Skeleton code

clear all, close all, clc

%% Reading in audio file
[audioReversed,fs] = audioread('audio1.wav');
% sound(audi0Reversed, fs)
%% Q1) The sound file 'audioReversed' appears to be reversed

% Transform audio file to frequency domain
AudioReversed = fft(audioReversed);

% Generate frequency vector f for plotting the magnitude spectrum
% in the next step
f = linspace(-fs/2, fs/2, length(audioReversed) + 1); f(end) = []; % This is the frequency vector from -fs/2 to fs/2

% Plot the magnitude of 'AudioReversed'
figure(1);
subplot(2,1,1);
plot(f, fftshift(abs(AudioReversed))/fs)
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude spectrum of AudioReversed');
grid on;

% Plot the phase of 'AudioReversed'
figure(1)
subplot(2,1,2)
plot(f, fftshift(angle(AudioReversed)))
xlabel('Frequency (Hz)')
ylabel('Phase')
title('Phase spectrum of Audio Reversed')
grid on;

% Using the relation found in the question sheet attempt to reverse
% 'AudioReversed' back to normal.

Audio = conj(AudioReversed);

% Plot the magnitude of Audio
figure(2)
subplot(2,1,1)
plot(f, fftshift(abs(Audio))/fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of Audio')
grid on;

% Plot the phase of Audio
subplot(2,1,2)
plot(f, fftshift(angle(Audio)))
xlabel('Frequency (Hz)')
ylabel('Phase')
title('Phase spectrum of Audio')
grid on;

% Transform 'Audio' back to the time domain
audio = ifft(Audio);

% Listen to audio file using the SOUND function
sound(audio,fs)

% What can you say about the two plots ?

%% Question 2

% Reading in audio file 'audioTone'
[audioTone,fs] = audioread('audio2.wav');

% Transform audio file to frequency domain
AudioTone = fft(audioTone);

% Generate frequency vector f for plotting the magnitude spectrum
% in the next step
f = linspace(-fs/2, fs/2, length(audioReversed) + 1); f(end) = []; % This is the frequency vector from -fs/2 to fs/2


% Plot the magnitude of 'AudioTone'
figure(3)
subplot(2,1,1)
plot(f, fftshift(abs(AudioTone))/fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of AudioTone')
grid on;

% Plot the phase of 'AudioTone'
subplot(2,1,2)
plot(f, fftshift(angle(AudioTone)))
xlabel('Frequency (Hz)')
ylabel('Phase')
title('Phase spectrum of AudioTone')
grid on;

% Recreate the tone
amplitude = 5e-3; %mV
frequency = 1500;
phase = 0.295;

Ts = 1/fs;
t = 0:Ts:length(audioTone)*Ts;
t(end) = [];
tone = amplitude*cos(2*pi*frequency*t + phase);

% Denoise the actual signal
denoised_audioTone = AudioTone - tone.';
Denoised = fft(denoised_audioTone);

% play the sound
sound(denoised_audioTone,fs);

% Plot the magnitude of 'Denoised'
figure(3)
subplot(2,1,1)
plot(f, fftshift(abs(denoised_audioTone))/fs)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude spectrum of AudioTone')
grid on;

% Plot the phase of 'Denoised'
subplot(2,1,2)
plot(f, fftshift(angle(denoised_audioToneclc
)))
xlabel('Frequency (Hz)')
ylabel('Phase')
title('Phase spectrum of AudioTone')
grid on;