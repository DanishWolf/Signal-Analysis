%% Assignment 1 - Part B, Section B2 (Combining Signals)
%  Do not change before line 38
%  You will need to have generated Data1B.mat from 
%  GenerateDataAssignment1B.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from Data1B.mat
load('Data1B.mat', 'fs', 'muxSignal');% 

%  The variable loaded are:
%    fs         Sampling frequency for Section B2
%    muxSignal  Multiplexed signals, for Section B2

%==================================================================
%
% Names of variables you will need for marking.
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the list below.
% ---------------------------------------------
% ===== Part 2 =====
%   Ts              Sampling period
%   t               time vector
%   MUXF            Fourier transform of muxSignal
%   k               Frequency vector
%   freqshifting    Frequency shifts (vector)
%   MagSpec         amplitude of peaks (vector)
%   PhaseSpec       phases of peaks (vector)
%   xdm             output of FDMDemux (matrix)
%   XDM             Fourier transform of xdm (matrix)
%   B               Bandwidth
%   filteredoutput  Filtered Output Signal
%   decodedtext     The Decoded Text Output
% ---------------------------------------------
%====Enter your code below this line================================
%% Question 1


