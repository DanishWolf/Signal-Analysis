%% Tutorial 2 Take Home Lab
function y = myfunc(T)


samples = 100;
f = linspace(-0.5,0.5,samples+1); f=f(1:end-1);


y = sin(pi*f*T)./(pi*f);

figure
plot(f,y,'m');
xlabel('Frequency [Hz]');
ylabel('Mangitude');
title('Plot of f vs y');
