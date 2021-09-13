%% Tutorial 1 Matlab Questions
% Question 3
%a
x = [10:20];

%b
y = [30:40];

%c
z = x + y;

%d
z_elementwise = x.*y;

%e
a = x*y.';

%f
b = x.'*y;

%g
size_a = size(a); 
size_b = size(b);

%% Question 4
%a
x = [2 5 1 6];

x_sqrt = sqrt(x); % x.^(0.5)

%b
y = x;

y(1:2:end) = y(1:2:end) * 1j;
%c
z = x+y;
z_re = real(z);
z_im = imag(z);

%d
z_ab = abs(z);


%% Question 5

t = linspace(0, 1, 1000 + 1); t(end) = [];

y1 = log(2 + t + t.^2); 
y2 = exp(t).*(1 + cos(3*t));
y3 = cos(t).^2 + sin(t).^2;



