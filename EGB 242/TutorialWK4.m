%% Week 4 In Class Challenge
clear all, close all, clc
%% Using fcoefs and cfseries
N = 3;
trigN = 1:N;
compN = -N:N;

an = zeros(1, N);

bn = -1./(pi*trigN);

a0 = 0.5;

cn = 1j./(2*pi*compN);
cn(N+1) = 0.5;

[y, t] = cfseries(cn, 1);

[s, t1] = fseries(a0, an, bn, 1);

%% 1
samples = 100;
T = 4;

t1 = linspace(0, T, samples+1) ; t1(end) = [];

y1 = zeros([1 samples]);
y1(t1 < 1) = 0.5;
y1(t1 < 2 & t1 >= 1) = 2;
y1(t1 < 3 & t1 >= 2) = 0;
y1(t1 >= 3) = 1;

figure
plot(t1, y1)

%% 2
R = 5;
s1 = repmat(y1, [1 R]);

%% 3
t = linspace(0, R*T, R*samples + 1); t(end) = [];



%% 4
f0 = 1/T;
Harm = 5;
n_trig = 1:Harm;
wnt = 2*pi*f0*t1.'*n_trig;
Ts = t1(2) - t1(1);

a0 = 1/T * sum(y1) * Ts;
an = 2/T * y1 * cos(wnt) * Ts;
bn = 2/T * y1 * sin(wnt) * Ts;


%% 5
% Method 1
c0 = a0;
cn_pos = 1/2 * (an - 1j*bn);
cn_neg = 1/2 * (an + 1j*bn);
cn = [fliplr(cn_neg) c0 cn_pos];

% Method 2
n_exp = -Harm:Harm;
c0 = 1/T * sum(y1) * Ts;
cn = 1/T * y1 * exp(-1j*2*pi*f0*t1.'*n_exp) *Ts;
cn(n_exp == 0) = c0;


%% 6
figure
subplot(2, 1, 1)
stem(n_exp, abs(cn))
subplot(2, 2, 2)
stem(n_exp, angle(cn))



%% 7-8
wnt = 2*pi*f0*t.'*n_trig;
f1  = cn * exp(1j*2*pi*f0*t.'*n_exp).';
f2  = a0 + an * cos(wnt).' + bn * sin(wnt).';


%% 9
figure
hold on
plot(t, real(f1), 'b')
plot(t, f2, 'g--')
plot(t, s1)
title('Fourier Series Approximation of Y1')
xlabel('Time (s)')
ylabel('Amplitude')



