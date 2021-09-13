%% Queensland University of Technology
%  EGB242 - Signal Analysis
%  Computer Lab Week 3.
%  Solution written by: 2013/2 team
%  Confirmed 2015
clear all % clears all existing variables
clc % clears the command window
close all % close all opened figures

% Highlight parts of the code and then hit the F9 key to run the selected
% areas.
% F5 key runs the entire script.

%% Question 1

% A cosine/sine waveform is a continuous function, made up of an infinite
% amount of samples. MATLAB however, operates on a discrete
% representation with a predefined finite amount of samples. To represent
% a cosine/sine waveform in MATLAB, we need to declare the limits of the
% time vector, i.e. from 0 to 0.5 secs, from 0 to 5 secs, from 2 to 10
% secs, with an amount of sample points. Having more sample points would
% produce in a better resolution, but result in a longer computation.

% The fundamental frequency, f0, is lowest frequency of a periodic waveform.
% It is the first harmonic of the waveform. The second harmonic is 2*f0,
% the third is 3*f0, so forth. All calculations are done with respect to
% the fundamental frequency. The different frequencies are all integer
% multiples of the fundamental frequency.

f0 = 50;
T = 1/50;
samples = 100;
t = linspace(0,5*T,5*samples+1);
t = t(1:end-1);

% RMS values are usually given. Peak values are used for waveform
% calculations.
s0 = sqrt(2)*240*cos(2*pi*f0*t);   %                             (a)
s1 = sqrt(2)*60*cos(2*pi*(f0*4)*t); %f0=50, f0*4=200Hz           (b)
s2 = sqrt(2)*100*cos(2*pi*(f0*10)*t); %f0=50, f0*10=500Hz        (c)

st = s0 + s1 + s2; %                                                 (d)

figure(1)
subplot(2,2,1)
plot(t,s0);
title('Plot of 50 Hz 240 volts cosine waveform')
xlabel('time')
ylabel('amplitude')
subplot(2,2,2)
plot(t,s1);
title('Plot of 200 Hz 60 volts cosine waveform')
xlabel('time')
ylabel('amplitude')
subplot(2,2,3)
plot(t,s2);title('Plot of 500 Hz 100 volts cosine waveform')
xlabel('time')
ylabel('amplitude')
subplot(2,2,4)
plot(t,st)
title('Plot of st')
xlabel('time')
ylabel('amplitude')

% (e)  Our function is even.
% (f)  If sine waveform was used instead, it would be an odd function.

%% Question 2

% Function fcoefs.p will return output arguments of a0, an, bn.
% The function will require input values of st, f0, samples, n.
% Because we are calculating the Fourier series coefficients for ONE
% period, please note that when approximating the waveform in the for loop
% summation below, the time vector has to be for ONE period.

% a)
N = 6; % number of coefficients required
[a0,an,bn] = fcoefs(st,f0,samples,N);

% (b)
st1 = zeros(1,samples); % note single period
t = linspace(0,1/f0,samples+1); % creating a new time vector of ONE period
t = t(1:end-1);

for m = 1:N
    st1 = st1 + an(m)*cos(2*pi*m*f0*t) + bn(m)*sin(2*pi*m*f0*t);
end

st1 = st1 + a0; % adding the DC component
st1 = repmat(st1,1,5); % expand the waveform to 5 periods

t = linspace(0,5*T,5*samples+1);
t = t(1:end-1); %rebuilding the time vector for plotting

figure(2)
plot(t,st1)
title('Plot of st1')
xlabel('time')
ylabel('amplitude')

%% Question 3

% a)
N = 20; % number of coefficients required
[a0,an,bn] = fcoefs(st,f0,samples,N);

% (b)
st2 = zeros(1,samples); % creating an array of appropriate length
t = linspace(0,1/f0,samples+1); % creating a new time vector of ONE period
t = t(1:end-1);

% reducing the number of loops used for faster computation
for m = 1:20
    st2 = st2 + an(m)*cos(2*pi*m*f0*t) + bn(m)*sin(2*pi*m*f0*t);
end

st2 = st2 + a0; % adding the DC component
st2 = repmat(st2,1,5); % expand the waveform to 5 periods, WITH FUNCTION
t = linspace(0,5*T,5*samples+1);
t = t(1:end-1); % rebuilding the time vector for plotting

figure(3)
plot(t,st2)
title('Plot of st2')
xlabel('time')
ylabel('amplitude')

% (c)
figure(4)
plot(t,st,'k')
hold on
plot(t,st1,'r')
plot(t,st2,'c--')
hold off
title('Overlapped plots of st, st1 & st2')
xlabel('time')
ylabel('amplitude')
legend('st=black','st1=red','st2=cyan','location','NorthEast')

% (d)
% As visible from figure 4, the st with 5 harmonics  was not able to
% approximate st, while st2 with 20 harmonics was able to give a better
% approximation which is somewhat similar to the original waveform.

% Try out and get the minimum amount of harmonics needed to approximate the
% orignal waveform.


