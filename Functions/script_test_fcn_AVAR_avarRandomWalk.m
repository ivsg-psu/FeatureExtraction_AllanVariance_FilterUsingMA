%% script_test_fcn_AVAR_avarRandomWalk.m
% This script tests 'fcn_AVAR_avarRandomWalk'
%
% This script was written on 2021_05_28 by Satya Prasad
% Questions or comments? szm888@psu.edu
% Updated: 2022/02/15

%% Prepare workspace
clear all %#ok<CLALL>
close all
clc

%% Intialization
rng('default') % set random seeds

number_of_time_steps = 16385;
p = floor(log2(number_of_time_steps));
number_of_time_steps = 2^p+1;
list_of_correlation_intervals = 2.^(0:(p-1))'; % list of correlation intervals

%% Example 1: Regularly Sampled Random Walk
random_walk_coefficient = 0.025; % [unit/sqrt(s)]
sampling_frequency      = 20; % [Hz]
list_of_correlation_time = list_of_correlation_intervals/sampling_frequency;

random_walk = fcn_AVAR_generateRandomWalk(random_walk_coefficient, ...
    sampling_frequency, number_of_time_steps); % generate random walk
% estimate AVAR
avar_estimated = fcn_AVAR_avar(random_walk, list_of_correlation_intervals,12345); % Normal
% actual AVAR
avar_calculated = fcn_AVAR_avarRandomWalk(random_walk_coefficient, ...
    list_of_correlation_time,12346);

fcn_AVAR_plotCompareAvar('Calculated',avar_calculated,'Estimated',avar_estimated,...
                         list_of_correlation_time,12347); % Plot for AVAR
