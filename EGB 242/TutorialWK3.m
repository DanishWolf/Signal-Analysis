close all, clear, clc
% Week 3 in class challenge
%%Part A

T  = 10;
samples = 100;
t_start = -T/2;
t_end = T/2;

t1 = linspace(t_start, t_end, samples + 1); t1(end) = [];

y1 = sin(pi*t1) ./ (pi*t1);

y1(t1 == 0) = 1;

y2 = (t1/10).^3;

s1 = repmat(y1,1,10);
s2 = repmat(y2,1,10);

t2_end = 10*T - T/2;

t2 = linspace(t_start,t2_end,T*samples+1); t2(end) = [];

figure()
hold on
plot(t2,s1)
plot(t2,s2)
hold off
title('Plot of s1 and s2')
xlabel('Time (s)')
ylabel('Amplitude')
xlim([t_start,t2_end]);
ylim,([-0.2,1]);






