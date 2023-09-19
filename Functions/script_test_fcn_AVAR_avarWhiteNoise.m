%% script_test_fcn_AVAR_avarWhiteNoise.m
% This script tests 'fcn_AVAR_avarWhiteNoise'
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
list_of_correlation_intervals = 2.^(0:(p-1))'; % list of correlation intervals

%% Example 1: White Noise
power_spectral_density = 0.0025; % PSD of white noise [unit^2 s]
sampling_frequency     = 20; % [Hz]
list_of_correlation_time = list_of_correlation_intervals/sampling_frequency;

white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density, ...
    sampling_frequency, number_of_time_steps); % generate white noise
% estimate AVAR
avar_estimated  = fcn_AVAR_avar(white_noise, list_of_correlation_intervals, 12345);
% actual AVAR
avar_calculated = fcn_AVAR_avarWhiteNoise(power_spectral_density, list_of_correlation_time, 12346);

fcn_AVAR_plotCompareAvar('Calculated',avar_calculated,'Estimated',avar_estimated,...
                         list_of_correlation_time,12347); % Plot for AVAR
