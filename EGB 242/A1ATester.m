%% Assignment 1A Tester
clear all, close all, clc
%% Using fcoefs and cfseries
N = 2;
trigN = 1:N;
compN = -N:N;

an = (3*sin((3*pi*trigN)/2)/(pi*trigN)) - (3*sin((pi*trigN)/2)/(pi*trigN));

bn = (-2*sin((pi*trigN)/2)/(pi^2 * trigN^2)) + ((3/(pi*trigN))*(cos((pi*trigN)/2) - cos((3*pi*trigN)/2)));

a0 = 3/2;

cn = (((exp(3) - exp(3/2 + 1j*pi*compN))/(3 - 1j*2*pi*compN)) + ((1/(1j*pi*compN)) * (exp(-1j*pi*compN) - 1)));
cn(N+1) = ((exp(3) - exp(-3/2 + 3))/3  - 1);

[y, t] = cfseries(cn, 1);

[s, t1] = fseries(a0, an, bn, 1);


