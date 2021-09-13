clc

syms t n

T = 1; %T is the period

f0 = 1/T; % f0 is the fundamental frequency

n = 2;

a0 = (1/T)*(int(-2*t,t,-1/4,1/4) + int(3,t,1/4,3/4))

an = (2/T)*(int(-2*t*cos(2*pi*n*f0*t),t,-1/4,1/4) + (int(3*cos(2*pi*n*f0*t),t,1/4,3/4)))
 
bn = (2/T)*(int(-2*t*sin(2*pi*n*f0*t),t,-1/4,1/4) + (int(3*sin(2*pi*n*f0*t),t,1/4,3/4)))

c0 = (1/T)*( (int((exp(3*t+3) ),t,-0.5,0)) + (int((-2),t,0,0.5)) )
 
cn = (1/T)*( (int((exp(3*t+3))*(exp(-2*pi*i*n*t)),t,-0.5,0)) + (int((-2)*(exp(-2*pi*i*n*t)),t,0,0.5)) )

%cn = (1/T)*int((exp(3*t) - 2 )*(exp(-1j*2*pi*n*t)),t,0,1)
 
%cn = (1/(3 - 2*i*pi*n*1)) - (  exp(1.5 + i*pi*n*1)  / (3 - 2*i*pi*n*1)) - (  exp(-1.5 - i*pi*n*1) /  (3 + 2*i*pi*n*1)  ) + ( 1/(3 + 2*i*pi*n*1) ) + ((exp(-i*pi*n*1))/(i*pi*n*1)) - (1/ (i*pi*n*1))
% s2t = (cn.*exp(1j*2*pi*n*t))

% sid = 10528873;
% [c0, cn, s2t_approx, sid] = compFS(sid);

% This line will output the two functions generated through generateDataAssignment1A.m
% [s1t, s2t] = generateFunction(sid)

% This line will display a typeset version of your input solution 
% displayResultQ2(c0, cn, s2t_approx)


%%
A = 5;
B = 3;
t = linspace(0, 5, 88200 + 1); t(end) = [];

s1 = t./A;
s2 = -exp((t-B)./A);
% s3 = A.*t - 1.25*A;
% s3(t > 5 & t >= 2.5) = -A.*t + 3.75*A;

figure
subplot(3, 1, 1)
plot(t,s1)
subplot(3,1,2)
plot(t,s2)
% subplot(3,1,3)
% plot(t,s3)





