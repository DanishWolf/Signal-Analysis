%% Tutorial Week 2 In Class Challenge
function [st, t] = gencos(f0, phi, R)
   % Check if R is an integer
   if ~isequal(round(R),R)
       error('R is not an integer');
   end
   T = 1/f0;
   
   samples = 100;
   
   t = linspace(0, R*T, R * samples + 1); t(end) = [];  
   
   st = cos((2*pi*f0*t) + phi);
   
   figure
   plot(t,st)
   xlabel('Time')
   ylabel('Amplitude')
   title('GENCOS Output st')
end