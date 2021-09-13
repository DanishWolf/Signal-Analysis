close all

L = 1;
samples = 100; % number of time sample
t = linspace(-L,L,samples+1); t(end) = []; % remove t = L due to periodicity
x = zeros(1,samples);
for k = 1:samples
    if t(k) > -L && t(k) < -L/2
        x(k) = 0;
    elseif t(k) >= -L/2 && t(k) <= L/2
        x(k) = 1;
    elseif t(k) > L/2 && t(k) < L 
        x(k) = 0;
    end
end

figure;
plot(t,x)

%% Trigonometric fourier Series
N = 5; % number of terms in simulation
a0 = 1/2;
xtf = a0*ones(1,samples);
for n = 1:N
    xtf = xtf + 2/(n*pi)*sin(n*pi/2)*cos(n*pi*t/L);
end
hold on
plot(t,xtf)

%% Complex Fourier Series
c0 = 1/2;
xcf = c0*ones(1,samples);
for n = [-N:-1,1:N]
%     if n == 0 
%         continue;
%     end
   xcf = xcf + 1/(n*pi)*sin(n*pi/2)*exp(1j*n*pi*t/L);
end
xcf = real(xcf);
plot(t,xcf,'--', 'k')