clear all
clc
cd('C:\Users\Iannelli\Desktop\MATLAB\DD_adaptive_LTV\GitHub-repo-files')
%% master file to reproduce results of paper
% "A hybrid systems framework for data-based adaptive control of linear time-varying systems"
% seeding randomness for reproducibility and debug
rng('default');
s = rng;

%%% Case studies:

CASE = 1; % 1 (switching plant: Figure 1), 2 (sinusoidal plant: Figure 2), 3 (cubic plant from ODDAC paper: Figure 3)

%%% Controller options

CONTROLLER = 3; 

% For CASE={1},{2} the following options are available:
% 0 (time-invariant controller designed with offline data), 
% 1 (proposed algorithm: event-triggered adaptive controller), 
% 2 (time-triggered adaptive controller)->see inside HS_DAC the period and starting time  

% For CASE={3} the following options are available:
% 1 (proposed algorithm: event-triggered adaptive controller), 
% 3 (robust time-triggered adaptive controller ODDAC) -> see inside HS_DAC tuning parameters 
% Note: CONTROLLER={3} will only run on CASE={3}! Because it was only implemented for a comparison with the plant analyzed in the paper that proposed ODDAC 



HS_DAC(CASE,CONTROLLER);





