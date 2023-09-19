%% script_test_fcn_AVAR_avarEmpiricalConfidenceCurves.m
% This script tests 'fcn_AVAR_avarEmpiricalConfidenceCurves' on different type of noise signals
%
% This script was written on 2022_02_07 by Satya Prasad
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
list_of_correlation_intervals = 2.^(1:(p-1))'; % list of correlation intervals

%% Example 1: White Noise
power_spectral_density = 0.0025; % PSD of white noise [unit^2 s]
sampling_frequency     = 20; % [Hz]
confidence_coefficient = 0.95;
noise_type             = 'wn';

white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density, ...
    sampling_frequency, number_of_time_steps); % generate white noise
white_noise = white_noise(number_of_time_steps-2^p:number_of_time_steps);
fcn_AVAR_avarEmpiricalConfidenceCurves(white_noise, list_of_correlation_intervals, ...
                                       confidence_coefficient, noise_type, 12345);

%% Example 2: Random Walk
random_walk_coefficient = 0.025; % [unit/sqrt(s)]
sampling_frequency     = 20; % [Hz]
confidence_coefficient = 0.95;
noise_type             = 'rw';

random_walk = fcn_AVAR_generateRandomWalk(random_walk_coefficient, ...
    sampling_frequency, number_of_time_steps); % generate random walk
random_walk = random_walk(number_of_time_steps-2^p:number_of_time_steps);
fcn_AVAR_avarEmpiricalConfidenceCurves(random_walk, list_of_correlation_intervals, ...
                                       confidence_coefficient, noise_type, 12346);
