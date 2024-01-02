%%%%%%%%%%%%%%%%%%%%% script_AVAR_effectOfNoiseRatio.m %%%%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to demonstrate the effect of noise ratio 
%   (power spectral density divided by square of random walk coefficient) 
%   on Allan Variance.
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
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]
list_of_ratios          = 10.^(-2:2:2)';
confidence_coefficient  = 0.95;
noise_type              = 'rw';

%% Synthesize the test signal
time_vector = sampling_interval*(0:number_of_time_steps-1);
white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
              sampling_frequency,number_of_time_steps); % White noise
random_walk = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
              sampling_frequency,number_of_time_steps); % Random walk
test_signal = random_walk + white_noise;

%% Estimate AVAR of the test signal
p = floor(log2(number_of_time_steps));
list_of_correlation_intervals = 2.^(0:p-3)'; % List of correlation intervals
list_of_correlation_time      = sampling_interval*list_of_correlation_intervals;

% AVAR estimated
[est_avar_test_signal,avar_lb,avar_ub] = ...
    fcn_AVAR_avarEmpiricalConfidenceCurves([test_signal; 0],...
    list_of_correlation_intervals,confidence_coefficient,noise_type);

%% Plot the results
%% AVAR of test signal with different noise ratios
legend_cell = cell(numel(list_of_ratios),1);
figure(01)
clf
width = 540; height = 448; right = 100; bottom = 100;
axis_position = [0.1354, 70.515/height, 0.7696, 307.32/height];
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
for i = 1:numel(list_of_ratios)
    calc_avar_test_signal = ...
        fcn_AVAR_avarRandomWalk(random_walk_coefficient,list_of_correlation_time) + ...
        fcn_AVAR_avarWhiteNoise(list_of_ratios(i)*power_spectral_density,list_of_correlation_time);
    if 1==i
        plot(list_of_correlation_intervals,calc_avar_test_signal','c--','Linewidth',1.2)
    elseif 2==i
        plot(list_of_correlation_intervals,calc_avar_test_signal','k','Linewidth',1.2)
    elseif 3==i
        plot(list_of_correlation_intervals,calc_avar_test_signal','m-.','Linewidth',1.2)
    end % NOTE: END IF statement
    legend_cell{i} = ['$\frac{C_{wn}}{C_{rw}^{2}} = $ ' num2str(list_of_ratios(i))];
end % NOTE: END FOR loop 'numel(list_of_ratios)'
legend(legend_cell,'Location','north','Interpreter','latex','FontSize',18)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([10^-4.5 10^0.5])
ax1 = gca;
axes('Position',axis_position,'XAxisLocation','top',...
     'xLim',ax1.XLim*sampling_interval,'XScale','log','ytick',[],...
     'xtick',[1e-1 1e1 1e3],'Color','none','Fontsize',13,'Box','off');
ax2 = gca;
xlabel(ax2,'Correlation Time $[s]$','Interpreter','latex','FontSize',18)

%% Plot test signal against time
min_yaxis   = min(test_signal(1:2000*sampling_frequency));
max_yaxis   = max(test_signal(1:2000*sampling_frequency));
range_yaxis = max_yaxis-min_yaxis;
figure(02)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(time_vector,test_signal,'k.','Markersize',1)
set(gca,'Fontsize',13)
ylabel('Amplitude $[Unit]$','Interpreter','latex','FontSize',18)
xlabel('Time $[s]$','Interpreter','latex','FontSize',18)
ylim([min_yaxis-0.1*(range_yaxis) max_yaxis+0.1*(range_yaxis)])
xlim([0 2000])

%% Plot AVAR of the test signal with confidence bounds
figure(03)
clf
width = 540; height = 448; right = 100; bottom = 100;
axis_position = [0.1354, 70.515/height, 0.7696, 307.32/height];
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(list_of_correlation_intervals,est_avar_test_signal,'k','Linewidth',1.2)
plot(list_of_correlation_intervals,avar_lb,'k--','Linewidth',1.2)
plot(list_of_correlation_intervals,avar_ub,'k-.','Linewidth',1.2)
legend('Estimated',[num2str(100*confidence_coefficient) '$\%$ LB'],...
       [num2str(100*confidence_coefficient) '$\%$ UB'],...
       'Location','best','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])
ax1 = gca;
axes('Position',axis_position,'XAxisLocation','top',...
     'xLim',ax1.XLim*sampling_interval,'XScale','log','ytick',[],...
     'xtick',[1e-1 1e1 1e3],'Color','none','Fontsize',13,'Box','off');
ax2 = gca;
xlabel(ax2,'Correlation Time $[s]$','Interpreter','latex','FontSize',18)
