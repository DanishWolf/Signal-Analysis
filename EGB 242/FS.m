function [a0, an, bn] = FS(s, Harm, f0)
% s = vector of values from your signal ('y'-values)
% Harm  = number of harmonics we want to calculate
% f0 = frequency of function

%Useful variables
T = 1/f0;
samples = length(s);

%time vector
t = linspace(0, 1, samples+1); t(end) = [];

%harmonics
n = 1:Harm;

%time step for riemanns sum
Ts = t(2) - t(1);

a0 = 1/T * sum(s) * Ts;
an = 2/T * s *cos(2*pi*n*f0.*t) * Ts;
bn = 2/T * s *sin(2*pi*n*f0.*t) * Ts; 

end

