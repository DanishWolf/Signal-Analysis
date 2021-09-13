%% Matlab Grader Question 1
function [a0, an, bn, s1t_approx, sid] = trigFS(sid)
    syms t, syms n, pi = sym('pi');
    % Enter your student ID as sid:
    sid = 10528873;
    % As a symbolic expression, save the expressions 
    % for the trigonometric coefficients below:
    a0 = 3/2
    an = (3*sin((3*pi*n)/2)/(pi*n)) - (3*sin((pi*n)/2)/(pi*n))
    bn = (-2*sin((pi*n)/2)/(pi^2 * n^2)) + ((3/(pi*n))*(cos((pi*n)/2) - cos((3*pi*n)/2)))
    % Find the trigonometric FS, evaluate the FS with
    % 2 harmonics and save the expression below:
    f = 1;
     nn = 2;
     a1 = (3*sin((3*pi*n)/2)/(pi*n)) - (3*sin((pi*n)/2)/(pi*n));
     b1 = (-2*sin((pi*n)/2)/(pi^2 * n^2)) + ((3/(pi*n))*(cos((pi*n)/2) - cos((3*pi*n)/2)));
     a2 = (3*sin((3*pi*nn)/2)/(pi*nn)) - (3*sin((pi*nn)/2)/(pi*nn));
     b2 = (-2*sin((pi*nn)/2)/(pi^2 * nn^2)) + ((3/(pi*nn))*(cos((pi*nn)/2) - cos((3*pi*nn)/2)));
    s1t_approx = a0 + a1 + b1 + a2 + b2
end


%% Matlab Grader Question 2
function [c0, cn, s2t_approx, sid] = compFS(sid)
    syms n t; e = sym(exp(sym(1))); pi = sym('pi');
    % Enter your student ID as sid:
    sid = 10528873;
    % As a symbolic expression, save the expressions 
    % for the complex coefficients below:
    c0 = ((exp(3) - exp(-3/2 + 3))/3  - 1);
    cn = (((exp(3) - exp(3/2 + 1j*pi*n))/(3 - 1j*2*pi*n)) + ((1/(1j*pi*n)) * (exp(-1j*pi*n) - 1)));
    % Find the complex FS, evaluate the FS up to 
    % second harmonic and save the expression below:
    s2t_approx = (((exp(3) - exp(3/2 + 1j*pi*-2))/(3 - 1j*2*pi*-2)) + ((1/(1j*pi*-2)) * (exp(-1j*pi*-2) - 1))) + (((exp(3) - exp(3/2 + 1j*pi*-1))/(3 - 1j*2*pi*-1)) + ((1/(1j*pi*-1)) * (exp(-1j*pi*-1) - 1))) + ((exp(3) - exp(-3/2 + 3))/3  - 1) + (((exp(3) - exp(3/2 + 1j*pi))/(3 - 1j*2*pi)) + ((1/(1j*pi)) * (exp(-1j*pi) - 1))) + (((exp(3) - exp(3/2 + 1j*pi*2))/(3 - 1j*2*pi*2)) + ((1/(1j*pi*2)) * (exp(-1j*pi*2) - 1)));
end

%% Matlab Grader Question 3
function [t, s2_hinf, sid] = noiseFunc(sid);
    % Enter you student ID below:
    sid = 10528873;
    % Save the appropriate outputs below as defined above:
    N = 5;
    T = 1;
    samples = 100;
    time = linspace(-T/2 , T/2, samples+1); time(end) = [];
    t = linspace(0,5,500+1); t(end) = [];
    s2 = exp(3.*time + 3);
    s2(time >= 0) = -2;
    
    s2_hinf = repmat(s2, [1,N]);
end

%% Matlab Grader Question 4
function [a0, an, bn, s_approx, T] = trigFS(s_hinf, t, N)
Ts = t(2) - t(1);
T = t(end) - t(1) + Ts;
f0 = 1/T;
n_trig = 1:N;

a0 = 1/T * sum(s_hinf) * Ts;
an = 2/T * s_hinf * cos(2*pi*f0*n_trig'*t).' * Ts;
bn = 2/T * s_hinf * sin(2*pi*f0*n_trig'*t).' * Ts;

s_approx = a0 + an * cos(2*pi*f0*n_trig'*t) + bn * sin(2*pi*f0*n_trig'*t);
    
end

%% Matlab Grader Question 5
function [c0, cn, s_approx, T] = complexFS(s_hinf, t, N)
Ts = t(2) - t(1);
T = t(end) - t(1) + Ts;
f0 = 1/T;
n_comp = -N:N;

c0 = 1/T * sum(s_hinf) * Ts;

cn = 1/T *s_hinf * exp(1j * -2 * pi * f0 *n_comp' * t).' * Ts;

s_approx = cn * exp(1j * 2 * pi * f0 *n_comp' * t);
end
