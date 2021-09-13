%% Queensland University of Technology
%  EGB242 - Signal Analysis
%  Computer Lab W4.
%  Solution written by: 2013/2 team
%  Confirmed 2015
clear all % clears all existing variables
clc % clears the command window
close all % close all opened figures

% Highlight parts of the code and then hit the F9 key to run the selected
% areas.
% F5 key runs the entire script.

%% Question 1

% (Setup)
% It is always good practice to define the required parameters neatly on
% the top.
N = 10;
n = -N:N;
samples = 100;
T = 1;
f0 = 1/T;
t = linspace(0,T,samples+1);
t = t(1:end-1);    %  (b)
an = zeros(1,N); % Defining how big the empty array should be
bn = zeros(1,N);
cn = zeros(1,length(n)); % cn is 2*N+1, the extra point being the DC component
f1 = zeros(1,length(t));
f2 = zeros(1,length(t));

% (a)
% All operations does not include f0 for simplicity because f0=1.
% Expression for cn from hand calcs.
% Unlike the trigo coefs, the complex coefs goes from -N:N, N being the
% maximum number of harmonics. 
% We want a loop of 2*N+1 long, consisting of 1 to N harmonics and 
% accounting for the DC component. FOR loop indices start from 1 so an 
% offset is done to go from -N to N.
% When p=0, cn=0.6. 
for m = 1:length(n)
    p = m-(N+1);
    if p == 0
        cn(m) = 0.6;
    else
        cn(m) = (exp(-1j*0.4*pi*p)-1 + exp(-1j*2*pi*p) - exp(-1j*1.2*pi*p))/(-1j*2*pi*p);
    end
end

% (c)
% f1 needs to be of length(t) because the complex coefs are multiplied with 
% the exponential with a time vector. 
for m = 1:length(n)
    p = m-(N+1);
    f1 = f1 + cn(m)*exp(1i*2*pi*p*t);
end

%% Question 2

% (a)
% For m=1, an(1)=cn(12)+cn(10), which is the first positive and negative
% harmonics for cn. MATLAB array indices start from 1. The lowest
% harmonic value for cn is the first index of the array, the next lowest 
% harmonic value for cn is the second index of the array, so forth.
for m = 1:N
    an(m) = cn(m+N+1) + cn(N+1-m);
    bn(m) = 1i*(cn(m+N+1) - cn(N+1-m));
end
a0 = cn(N+1); % Extracting the corresponding DC component

% (b)
% Make sure your fcoefs.p is in the same directory
[a01,an1,bn1] = fcoefs(f1,f0,samples,N); 

% (c)
for m = 1:N
    f2 = f2 + an(m)*cos(2*pi*m*t) + bn(m)*sin(2*pi*m*t);
end
f2 = f2 + a0; % adding the DC component

figure(1)
subplot(2,1,1)
plot(t,f1,'r')
title('complex coefs plot')
xlabel('time')
ylabel('amplitude')
subplot(2,1,2)
plot(t,f2,'b')
title('trigo coefs plot')
xlabel('time')
ylabel('amplitude')

% (d)
% To plot the waveforms, we need to build the corresponding time vectors
t1 = linspace(0,T,samples+1);
t1 = t1(1:end-1);
t2 = linspace(0+T,T+T,samples+1);
t2 = t2(1:end-1);
t3 = linspace(0+2*T,T+2*T,samples+1);
t3 = t3(1:end-1);
t4 = linspace(0+3*T,T+3*T,samples+1);
t4 = t4(1:end-1);

figure(2)
plot(t1,f1,'r')
hold on
plot(t2,f2,'b')
plot(t3,f1,'r')
plot(t4,f2,'b')
hold off
title('Approximated waveform from trig (blue) and complex (red) coefs')
xlabel('time')
ylabel('amplitude')