%%%%%%%%%%%%%%%%%%% script_AVAR_demonstrateAreaOfAVAR.m %%%%%%%%%%%%%%%%%%%
%% Purpose
%   The purpose of this script is to demonstrate Trapezoidal and Euler area
%   under AVAR curve.
%
% This script was written on 2023_09_08 by Satya Prasad
% Questions or comments? szm888@psu.edu

%% Prepare the workspace
clear all %#ok<CLALL>
close all
clc

%% Initialization
rng('default')

%% Define inputs and other parameters
sampling_frequency   = 50; % [Hz]
number_of_time_steps = 2^19;

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]

%% Synthesize the test signal
white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
              sampling_frequency,number_of_time_steps); % White noise
random_walk = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
              sampling_frequency,number_of_time_steps); % Random walk
test_signal = random_walk + white_noise;

%% Estimate AVAR of the input signal
p = floor(log2(number_of_time_steps));
list_of_correlation_intervals = 2.^(0:p-3)'; % List of correlation intervals

% AVAR estimated
est_avar_test_signal = fcn_AVAR_favar([test_signal; 0],list_of_correlation_intervals);

%% Plot the results
%% Euler area under AVAR curve in log2 scale
figure(01)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(log2(list_of_correlation_intervals),est_avar_test_signal,'k','Linewidth',1.2)
plot(log2(list_of_correlation_intervals),est_avar_test_signal,'k*','Markersize',8)
for i = 1:numel(est_avar_test_signal)-1
    area(log2(list_of_correlation_intervals(i:i+1)),est_avar_test_signal(i)*ones(2,1),...
         'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);
end % NOTE: END FOR loop 'numel(est_avar_test_signal)-1'
xt = 0:4:p-3;
xticks(xt)
xticklabels(cellstr(num2str(xt(:),'2^{%d}')))
set(gca,'FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([-0.01 0.21])
xlim([0 p-3])

%% Trapezoidal area under AVAR curve in log2 scale
figure(02)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(log2(list_of_correlation_intervals),est_avar_test_signal,'k','Linewidth',1.2)
plot(log2(list_of_correlation_intervals),est_avar_test_signal,'k*','Markersize',8)
area(log2(list_of_correlation_intervals),est_avar_test_signal,...
     'EdgeColor','none','FaceColor','c','FaceAlpha',0.5);
xt = 0:4:p-3;
xticks(xt)
xticklabels(cellstr(num2str(xt(:),'2^{%d}')))
set(gca,'FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([-0.01 0.21])
xlim([0 p-3])
