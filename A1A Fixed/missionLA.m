%% Assignment 1 Part A - Section A2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.%
%  Do not change before line 30.
%  If you have not generated Data1A from GenerateDataAssignment1A.m,
%  do that now.
%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data.
load Data1A;  

% VARIABLES:
% k - Time vector
% kPeriod - Time vector for 1 period of the interfering noise waveform
% T - Period of signal
% addNoise1st - Your first noise waveform
% addNoise2nd - Your second noise waveform
% c0, cn - Complex Fourier series coefficient vectors 
% OR
% a0, an, bn - Trig Fourier series coefficient vectors 
% FS1st - Fourier series approximation vector of first noise waveform 
% FS2nd - Fourier series approximation vector of second noise waveform 
% FSTotal - Fourier series approximation vector of total noise waveform 
% dnMsg - De-noised resulting wave (with 5 harmonics)

%==================================================================
% Refer to the assignment sheet for details on variable naming.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
%====Enter your code below this line================================
%% Question 1: Plotting the Recieved
% Length of Noise Signal 
lgthNS = length(noiseSound);

% Time vector
t = linspace (0 ,20 ,lgthNS + 1);
t (end) = [ ];

% Plotting recieved noise
% Sound ( noiseSound ,44100);
figure(1)
plot(t, noiseSound) % Plot of the received noise over time
title('Received Speech Signal')
xlabel('Time (s)')
ylabel('Amplitude')


%% Question 2: Noise Signals Plotted over 1 Period
% Period
samples = lgthNS * 0.25; 
T = 5;
% Time Vector
kPeriod = linspace(0, T, samples + 1); kPeriod(end) = [];

% Time Vector in Halves
kPeriodHalf = linspace(0, T/2, samples/2 + 1); kPeriodHalf(end) = [];
kPeriodHalf2 = linspace(T/2, T, samples/2 + 1); kPeriodHalf2(end) = [];

% The Piecewise function s3(t)
PW1 = (A*kPeriodHalf -1.25*A);  % First half of the piecewise function
PW2 = (-A*kPeriodHalf2 + 3.75*A); % Second half of the piecewise function

% Addition of Noise Vectors
addNoise1st = [PW1 PW2]; % Piece Wise Function
addNoise2nd = kPeriod/A; % Upward Sloping Function

% Plot of the above
hold off
figure(2)
hold on

subplot(1, 2, 1)
plot(kPeriod, addNoise1st)
title('First Additive Noise')
xlabel('Time (s)')
ylabel('Amplitude')
legend('S3(t)')

subplot(1, 2, 2)
plot(kPeriod, addNoise2nd)

title('Second Additive Noise')
xlabel('Time (s)')
ylabel('Amplitude')
legend('S1(t)')
hold off


%% Question 3 Classifing Signals
% Results inside the report


%% Question 4 Fourier Series Approximations of AddNoise2nd
% Complete Complex Fourier Series

% Number of Harmonices
Harm = 5;
compN = -Harm:Harm;
compN(Harm + 1) = []; % Used to remove n = 0 as c0 is present

% Time step
Ts = kPeriod(2) - kPeriod(1);

% Fundamental Frequnecy
f0 = 1/T;

% Complex Fourier Series
c0 = 1/T * sum(addNoise2nd) * Ts;
cn = 1/T * addNoise2nd * exp(1j * -2 * pi * f0 * compN' * kPeriod).' * Ts;
FS2nd = c0 + cn * exp(1j * 2 * pi * f0 * compN' * kPeriod);

% Plotting the Complex Fourier Series
hold off
figure(3)
hold on

plot(kPeriod, real(FS2nd), 'b') % Only plots real numbers
plot(kPeriod, addNoise2nd, 'm')

title('Original Signal vs Complex Fourier Series Approximation')
xlabel('Time (s)')
ylabel('Magnitude')
legend('CFS Approximation', 'S1(t)')
hold off


%% Question 5
% Complex Fourier Series Calculated by Hand
CHc0 = 0;
CHcn = (7.5 .* (exp(1j .* -pi .* compN) - 1) ./ (pi.^2 .* compN.^2));

% Zeros Vector Required for For Loop
CHapprox = zeros(1, samples);

% For Loop
for n = compN
    CHapprox = CHapprox + (7.5 .* (exp(1j .* pi .* n) - 1) ./ (pi.^2 .* n.^2))...
        .* exp(1j .* 2 .* pi .* n .* f0 .* kPeriod);
end

% Add the value of n = 0
CHapprox = CHapprox + CHc0;

% Plotting the Complex Fourier Series Calculated by Hand
hold off
figure(4)
hold on

plot(kPeriod, real(CHapprox))
plot(kPeriod, real(addNoise1st), '--')

title('Hand Calculation vs Real Signal')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Hand Calculation', 'Real Signal')


%% Question 6
% Fourier Series Approximations of the Additional Noise
% Additional Noise Vectors
Noise1st = [PW1 PW2]; % Piece Wise Function S3(t)
Noise2nd = kPeriod/A; % Upward Sloping Function S1(t)

% Complex Fourier Series of the 1st Noise
FS1c0 = 1/T * sum(addNoise1st) * Ts;
FS1cn = 1/T * addNoise1st * exp(1j * -2 * pi * f0 * compN' * kPeriod).' * Ts;

FS1st = FS1c0 + FS1cn * exp(1j * 2 * pi * f0 * compN' * kPeriod);
%FS1st = repmat(FS1stapp, [1, 2]);

% Complex Fourier Series of the 2nd Noise
FS2c0  = 1/T * sum(Noise2nd) * Ts;
FS2cn = 1/T * Noise2nd * exp(1j * -2 * pi * f0 * compN' * kPeriod).' * Ts;

FS2nd = FS2c0 + FS2cn * exp(1j * 2 * pi * f0 * compN' * kPeriod);
%FS2nd = repmat(FS2ndapp, [1, 2]);

%kPeriod2

% Plotting the Additional Noise Approximation
hold off
figure(5)
hold on

subplot(1, 2, 1)
plot(kPeriod, real(FS1st))

title('First Additional Noise Approximation')
xlabel('Time (s)')
ylabel('Amplitude')
legend('-5 < n < 5')

subplot(1, 2, 2)
plot(kPeriod, real(FS2nd))

title('Second Additional Noise Approximation')
xlabel('Time (s)')
ylabel('Amplitude')
legend('-5 < n < 5')


%% Question 7
% The Fourier Series Total
FSTotal = [FS1st FS1st FS2nd FS2nd];

% Plot of Fourier Series Total
hold off
figure(6)
hold on

plot(t, real(FSTotal))

title('Fourier Series Total')
xlabel('Time (s)')
ylabel('Amplitude')


%% Question 8
% De-Noising the Original Signal
dnMsg = noiseSound - real(FSTotal).';


%% Question 9
% Results inside the report


%% Question 10
% Plotting the De-Noised Signal
hold off
figure(7)
hold on

plot(t, dnMsg)

title('The De-Noised Speech Signal')
xlabel('Time (s)')
ylabel('Amplitude')

load chirp.mat
sound(real(dnMsg), 44100)

