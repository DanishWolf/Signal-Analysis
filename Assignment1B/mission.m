%% Assignment 1 - Part B, Section B2 (Combining Signals)
%  Do not change before line 38
%  You will need to have generated Data1B.mat from 
%  GenerateDataAssignment1B.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from Data1B.mat
load('Data1B.mat', 'fs', 'muxSignal');% 

%  The variable loaded are:
%    fs         Sampling frequency for Section B2
%    muxSignal  Multiplexed signals, for Section B2

%==================================================================
%
% Names of variables you will need for marking.
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the list below.
% ---------------------------------------------
% ===== Part 2 =====
%   Ts              Sampling period
%   t               time vector
%   MUXF            Fourier transform of muxSignal
%   k               Frequency vector
%   freqshifting    Frequency shifts (vector)
%   MagSpec         amplitude of peaks (vector)
%   PhaseSpec       phases of peaks (vector)
%   xdm             output of FDMDemux (matrix)
%   XDM             Fourier transform of xdm (matrix)
%   B               Bandwidth
%   filteredoutput  Filtered Output Signal
%   decodedtext     The Decoded Text Output
% ---------------------------------------------
%====Enter your code below this line================================
%% B2.1

% Finding the Sampling Period
Ts = 1/fs;

% Creating the Time Vector
t = linspace(0, Ts*length(muxSignal), length(muxSignal) + 1); t(end) = [];

% Plotting muxSignal over Time
hold off
figure(1)
plot(t, muxSignal, 'b')
grid on
title('muxSignal in the Time Domain')
xlabel('Time')
ylabel('Magnitude')

%% B2.2

% Computing the Fourier Transformation of muxSignal
MUXF = fft(muxSignal);

% Creating the Frequency Vector of mux
k = linspace(-fs/2, fs/2, length(muxSignal) + 1); k(end) = [];

% Unshifting the Fourier Transformation of mux
MUXFshift = fftshift(MUXF);

% Plotting the Fourier Transformation over Frequency
hold off
figure(2)
plot(k, abs(MUXFshift)/fs, 'b')
grid on
title('muxSignal in the Frequency Domain')
xlabel('Frequency')
ylabel('Magnitude')

%% B2.3

% Find and Store the Frequency shift Peaks
% Using MATLAB's inbuilt function 'findpeaks'
% Locating the Frequency Shift Peaks
[Peaks, PeakLocations] = findpeaks(abs(MUXFshift)/fs, 'MinPeakHeight', 0.4);

% To get the positive values and sort in the correct order
index = sort(PeakLocations(7:12));

% Getting the required Frequency Shifted Values
freqshifting = k(index);

%% B2.4

% Find the Magnitude and Store in MagSpec and Find Phase and Store in
% PhaseSpec the Values in these vectors correspond to the peaks in
% freqshift
MagSpec = abs(MUXFshift(index))/fs;
PhaseSpec = angle(MUXFshift(index))/fs;

% Plot Locations and Magnitudes of these Frequency Shifts using Red Circles
% on the Magntidue Spectrum of MUXF
hold off
figure(3)
hold on
plot(k, abs(MUXFshift)/fs, 'b')
plot(freqshifting, MagSpec, 'ro')
grid on
title('Magnitude against Frequency')
xlabel('Frequency')
ylabel('Magnitude')

%% B2.5

% FDMDemux function implemented into MATLAB so that mission.m can run
% properly
    xdm = zeros(length(freqshifting),length(t));
    
    for i = 1:length(freqshifting)
        xdm(i,:) = muxSignal .* cos(2 * pi * freqshifting(i) * t + PhaseSpec(i)) * MagSpec(i); 
    end

% frequency shifting module called FDMDemux.m
[xdm] = FDMDemux(muxSignal, t, MagSpec, freqshifting, PhaseSpec);

%% B2.6

% Compute the Fourier transform for each data stream in xdm (row by row)
% Store the result in the matrix XDM.
XDM = zeros(size(xdm));

 for count = 1:size(xdm, 1)
     XDM(count, :) = fft(xdm(count, :));
     
 end   
     
% Plot the magnitude spectrum for each data stream. 
figure(4);
for count = 1: size(XDM)
   subplot(2, 3, count);
   plot(k, abs(fftshift(XDM(count, :)) / fs));
   title('Magnitudes from Data stream')
end

%% B2.7

% Storing the bandwidth value in the vector B (Bandwidth)
B = 3500;
B2 = 2000;
B3 = 5000;

% Creating the required filter
Filter = zeros(1,length(t));
Filter(-B<k &  k<B) = 1;

Filter2 = zeros(1,length(t));
Filter2(-B2<k &  k<B2) = 1;

Filter3 = zeros(1,length(t));
Filter3(-B3<k &  k<B3) = 1;


Filtered_Signal1 = fftshift(Filter).*(XDM(1,:));
Filtered_Signal2 = fftshift(Filter2).*(XDM(2,:));
Filtered_Signal3 = fftshift(Filter2).*(XDM(3,:));
Filtered_Signal4 = fftshift(Filter2).*(XDM(4,:));
Filtered_Signal5 = fftshift(Filter2).*(XDM(5,:));
Filtered_Signal6 = fftshift(Filter3).*(XDM(6,:));

figure(5)
subplot(6,1,1)
plot(k, fftshift(abs(Filtered_Signal1))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from Stream 1');

subplot(6,1,2)
plot(k, fftshift(abs(Filtered_Signal2))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from Stream 2');

subplot(6,1,3)
plot(k, fftshift(abs(Filtered_Signal3))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from Stream 3');

subplot(6,1,4)
plot(k, fftshift(abs(Filtered_Signal4))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from Stream 4');

subplot(6,1,5)
plot(k, fftshift(abs(Filtered_Signal5))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from Stream 5');

subplot(6,1,6)
plot(k, fftshift(abs(Filtered_Signal6))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from Stream 6');

figure(6)
subplot(6,1,1)
signal1 = ifft(Filtered_Signal1) - mean(ifft(Filtered_Signal1));
plot(t, signal1)
xlabel('Time (s)')
ylabel('Magnitude')
title('Decoded Signal Stream 1')

subplot(6,1,2)
signal2 = ifft(Filtered_Signal2) - mean(ifft(Filtered_Signal2));
plot(t, signal2)
xlabel('Time (s)')
ylabel('Magnitude')
title('Decoded Signal Stream 2')

subplot(6,1,3)
signal3 = ifft(Filtered_Signal3) - mean(ifft(Filtered_Signal3));
plot(t, signal3)
xlabel('Time (s)')
ylabel('Magnitude')
title('Decoded Signal Stream 3')

subplot(6,1,4)
signal4 = ifft(Filtered_Signal4) - mean(ifft(Filtered_Signal4));
plot(t, signal4)
xlabel('Time (s)')
ylabel('Magnitude')
title('Decoded Signal Stream 4')

subplot(6,1,5)
signal5 = ifft(Filtered_Signal5) - mean(ifft(Filtered_Signal5));
plot(t, signal5)
xlabel('Time (s)')
ylabel('Magnitude')
title('Decoded Signal Stream 5')

subplot(6,1,6)
signal6 = ifft(Filtered_Signal6) - mean(ifft(Filtered_Signal6));
plot(t, signal6)
xlabel('Time (s)')
ylabel('Magnitude')
title('Decoded Signal Stream 6')
%% B2.8

% Signal 1 
sound(signal1,fs)

% Signal 2
text2 = A1BTextdecode(signal2,fs)

% Signal 3
text3 = A1BTextdecode(signal3,fs)

% Signal 4
text4 = A1BTextdecode(signal4,fs)

% Signal 5
text5 = A1BTextdecode(signal5,fs)

% Signal 6
sound(signal6,fs)

decodedtext = [text2 text3 text4 text5];
