%% script_test_fcn_AVAR_favar.m
% This script tests 'fcn_AVAR_favar' on different type of noise signals.
%
% This script was written on 2021_05_15 by Satya Prasad
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
sampling_frequency   = 20; % [Hz]

white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density, ...
    sampling_frequency, number_of_time_steps); % generate white noise
white_noise = white_noise(number_of_time_steps-2^p:number_of_time_steps);
fcn_AVAR_favar(white_noise, list_of_correlation_intervals, 12345);

%% Example 2: Random Walk
random_walk_coefficient = 0.025; % [unit/sqrt(s)]
sampling_frequency      = 20; % [Hz]

random_walk = fcn_AVAR_generateRandomWalk(random_walk_coefficient, ...
    sampling_frequency, number_of_time_steps); % generate random walk
random_walk = random_walk(number_of_time_steps-2^p:number_of_time_steps);
fcn_AVAR_favar(random_walk, list_of_correlation_intervals, 12346);

%% Example 3: Error in input 'data'
power_spectral_density = 0.0025; % PSD of white noise [unit^2 s]
sampling_frequency   = 20; % [Hz]

white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density, ...
    sampling_frequency, number_of_time_steps); % generate white noise
fcn_AVAR_favar([1; white_noise], list_of_correlation_intervals);

%% Example 4: Error in 'correlations_intervals'
power_spectral_density = 0.0025; % PSD of white noise [unit^2 s]
sampling_frequency   = 20; % [Hz]

white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density, ...
    sampling_frequency, number_of_time_steps); % generate white noise
white_noise = white_noise(number_of_time_steps-2^p:number_of_time_steps);
fcn_AVAR_favar(white_noise, [list_of_correlation_intervals; 100]);
