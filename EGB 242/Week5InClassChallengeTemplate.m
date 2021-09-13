% Week 5 Code challenge
clear;
close all;
clc;

% Note some common variable naming standards:
% Start with lower case letter for time domain signal
% i.e x, y
% Start with Upper case letter for frequency domain signal
% i.e. X, Y

%% Part A
% --- Q1 ---
% Define commonly used variables and magic numbers
T = 10;
samples = 2000;
% Time vector goes from 0 to 10 seconds
tA = linspace(0, T, samples+1); tA(end) = [];


% --- Q2 ---
% Define commonly used variables and magic numbers
f = [5 14 17 40 63];
n = 1:5;
% Calculate the composite waveform
phase = n*pi/4;
phase = repmat(phase.', [1 length(tA)]);
x = cos(2*pi*f.'*tA + phase);
x = sum(x);

figure
plot(tA, x)

% --- Q3 --- 
% Find sampling frequency to be able to determine the frequency vector
TsA = tA(2) - tA(1);
fsA = 1/TsA;
% Calculate frequency vector
fA = linspace(-fsA/2, fsA/2, samples + 1); fA(end) = [];

% --- Q4 --- 
% Transform to frequency domain
X = fft(x);
% NOTE: Careful to note that if 'x' is a 1 x samples we can fft this straight 
% up as shown above but if 'a' is a matrix of signals, say 
% num_signals x samples, we must fft each row individually and NOT do fft(a)

% --- Q5 --- 
% NOTE:
% Divide by scaling factor when viewing the magnitude spectrum as MATLAB does
% scaling when using fft(). 
% Divide by length of signal if the signal is simple and made up of constant
% pure sinosoids
% (i.e. cos(...) + sin(...)
% Divide by fs if otherwise
% (i.e. a speech signal)

% Plot the signal in the time and frequency domains
figure
subplot(3, 1, 1)
plot(tA, x)
title('Time Domain')

subplot(3, 1, 2)
plot(fA, fftshift(abs(X))/fsA)  % Only use fftshift when plotting 
title('Frequency Domain')

subplot(3, 1, 3)
plot(fA, fftshift(angle(x)))
title('Frequency Domain Phase')
% Lets have a closer look at the time domain as it's hard to see what's
% happening


%% Part B

% --- Q1 ---
% Define commonly used variables and magic numbers
tB = linspace(-50, 50, samples + 1); tB(end) = [];
% Time vector goes from -50 to +50 seconds
fc = 5;
y = sinc(tB);
yshift = y .* cos(2*pi*fc*tB);
yshift2 = yshift .* cos(2*pi*fc*tB);
% --- Q2 --- 
% Define commonly used variables and magic numbers

% Generate the signal and its shifted version


% --- Q3 --- 
% Find sampling frequency to be able to determine the frequency vector
TsB = tB(2) - tB(1);
fsB = 1/TsB;
% Calculate the frequency vector 
fB = linspace(-fsB/2, fsB/2, samples + 1); fB(end) = [];

% --- Q4 --- 
% Transform to frequency domain
Y = fft(y);
YShift = fft(yshift);
YShift2 = fft(yshift2);


% --- Q5 --- 
% Remember from before:
% Divide by scaling factor when viewing the magnitude spectrum as MATLAB does
% scaling when using fft(). 
% Divide by length of signal if the signal is simple and made up of constant pure sinosoids
% (i.e. cos(...) + sin(...)
% Divide by fs if otherwise
% (i.e. a speech signal)

figure
subplot(3, 1, 1)
plot(fB, fftshift(abs(Y))/fsB)
title('Y')
xlabel('Frequnecy')
ylabel('Magnitude')

subplot(3, 1, 2)
plot(fB, fftshift(abs(YShift))/fsB)
title('YShift')
xlabel('Frequnecy')
ylabel('Magnitude')

subplot(3, 1, 3)
plot(fB, fftshift(abs(YShift2))/fsB)
title('YShift2')
xlabel('Frequnecy')
ylabel('Magnitude')