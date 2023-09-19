%%%%%%%%%%%%%%%%% script_AVAR_effectOfSamplingFrequency.m %%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to demonstrate the effect of sampling 
%   frequency on AVAR.
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
list_of_sampling_frequency   = [25, 50, 100];
list_of_number_of_time_steps = [2^18, 2^19, 2^20];

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]

%% Plot the results
%% AVAR of test signals with different sampling frequency
legend_cell = cell(numel(list_of_sampling_frequency),1);
figure(01)
clf
width = 1056.2; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
for i = 1:numel(list_of_sampling_frequency)
    sampling_frequency   = list_of_sampling_frequency(i); % [Hz]
    number_of_time_steps = list_of_number_of_time_steps(i);
    
    p = floor(log2(number_of_time_steps));
    list_of_correlation_intervals = 2.^(0:p-3)'; % List of correlation intervals
    list_of_correlation_time      = list_of_correlation_intervals/sampling_frequency;
    % Calculate AVAR of the signal
    avar_test_signal = ...
        fcn_AVAR_avarRandomWalk(random_walk_coefficient,list_of_correlation_time) + ...
        fcn_AVAR_avarWhiteNoise(power_spectral_density,list_of_correlation_time);
    
    legend_cell{i} = [num2str(sampling_frequency) 'Hz'];
    subplot(1,2,1)
    if i==1
        hold on
        grid on
        plot(list_of_correlation_time,avar_test_signal'+(i-2)*0.0001,'c--','Linewidth',1.2)
    elseif i==2
        plot(list_of_correlation_time,avar_test_signal'+(i-2)*0.0001,'k','Linewidth',1.2)
    elseif i==3
        plot(list_of_correlation_time,avar_test_signal'+(i-2)*0.0001,'m-.','Linewidth',1.2)
        axis_position = [75/width, 0.1567, 415.6/width, 0.7683];
        legend(legend_cell,'Location','best','Interpreter','latex','FontSize',13)
        set(gca,'Position',axis_position,'XScale','log','YScale','log','FontSize',13)
        ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
        xlabel('Correlation Time $[s]$','Interpreter','latex','FontSize',18)
        ylim([1e-4 1e0])
        title('(a)','Interpreter','latex','FontSize',18)
    end % NOTE: END IF statement
    
    subplot(1,2,2)
    if i==1
        hold on
        grid on
        plot(list_of_correlation_intervals,avar_test_signal','c--','Linewidth',1.2)
    elseif i==2
        plot(list_of_correlation_intervals,avar_test_signal','k','Linewidth',1.2)
    elseif i==3
        plot(list_of_correlation_intervals,avar_test_signal','m-.','Linewidth',1.2)
        axis_position = [(2*75+415.584)/width, 0.1567, 415.6/width, 0.7683];
        set(gca,'Position',axis_position,'XScale','log','YScale','log','FontSize',13)
        xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
        ylim([1e-4 1e0])
        title('(b)','Interpreter','latex','FontSize',18)
    end % NOTE: END IF statement
end % NOTE: END FOR loop 'numel(list_of_sampling_frequency)'
