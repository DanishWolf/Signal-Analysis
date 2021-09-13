%% Queensland University of Technology
%  EGB242 - Signals Analysis
%  Computer Lab Week 2.
%  Solution written by: 2013/2 team
%  Confirmed 2015
clear all % clears all existing variables
clc % clears the command window
close all % close all opened figures

% Highlight parts of the code and then hit the F9 key to run the selected
% areas.
% F5 key runs the entire script.

%% Question 1
% My first Matlab script.
T = 1;
samples = 100;
% As shown previously in Computer Lab Wk 1, the correct representation for
% f = [-0.5 0.5) is:
f = linspace(-0.5,0.5,samples+1);
f = f(1:end-1);
y = sin(pi*f*T)./(pi*f);

figure
plot(f,y)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Plot of y vs f')
legend('y function')


%% Question 2
% All function and the respective file name MUST be the same if not the
% function WILL NOT work.
% function [1,2..]  =  FunctionName(A,B,C..), written IN SCRIPT
% [1,2..]  =  FunctionName(A,B,C..), written IN MATLAB COMMAND WINDOW
% function WILL NOT work with insufficient input arguments

% copy the below code into a separate file and save as "myfunc.m"
function [y] = myfunc(T)
    samples = 100;
    f = linspace(-0.5,0.5,samples+1);
    f = f(1:end-1);
    y = T*sin(pi*f*T)./(pi*f*T);
    figure(1)
    plot(f,y)
end

%% Question 3
% NaN is returned when the computation output is not defined.
% An error message is displayed if input T is equal to 0,

% copy the below code and update "myfunc.m"
function [y] = myfunc(T)
if T == 0
    error('Invalid input: T cannot be zero.');
else
    samples = 100;
    f = linspace(-0.5,0.5,samples+1);
    f = f(1:end-1);
    y = T*sin(pi*f*T)./(pi*f*T);
    figure(1)
    plot(f,y)
end

%% Question 4
% Compute the first 15 numbers of the Fibonacci sequence
% Arrays in MATLAB starts from the index 1
fseq = zeros(1,15);
fseq(1) = 0;
fseq(2) = 1;
% fseq(x) will call out the xth index in the array

% In this case, the for loops will run from ii = 3,4,5,6,7,8,9,10,11,12,13,
% 14,15. 13 loops in total. Components within will operate on the index it is on.
% i.e. if ii=5, fseq(5) = fseq(5-2)+fseq(5-1);
% i.e. if ii=5, fseq(5) = fseq(3)  +fseq(4);
for ii=3:15
    fseq(ii) = fseq(ii-2)+fseq(ii-1);
end
end