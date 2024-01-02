%%%%%%%%%%%%%%%%%%%%%% script_AVAR_mseRWinMAfilter.m %%%%%%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to show that MSE of a MA filter with
%   random walk input corrupted by white noise is independent of time.
%
% This script was written on 2024/01/02 by Satya Prasad
% Questions or comments? szm888@psu.edu

%% Prepare the workspace
clear all %#ok<CLALL>
close all
clc

%% Initialization
rng('default')

%% Define inputs and other parameters
number_of_iterations = 100;

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals   = 2.^(0:p-3)'; % List of correlation intervals
number_of_correlation_intervals = numel(list_of_correlation_intervals);
target_points = [round(number_of_time_steps/3), round(2*number_of_time_steps/3), ...
                 number_of_time_steps];
number_of_target_points = numel(target_points);

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]

% Initialize variables to store MSE
estimated_MSE  = zeros(number_of_time_steps,number_of_correlation_intervals);
calculated_MSE = NaN(number_of_correlation_intervals,1);
for index_iter = 1:number_of_iterations
    %% Synthesize the input signal
    if index_iter==1
        time_vector = sampling_interval*(0:number_of_time_steps-1)';
    end % NOTE: END IF statement 'index_iter==1'
    white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,sampling_frequency,...
                   number_of_time_steps+list_of_correlation_intervals(end)-1); % White noise
    random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,sampling_frequency,...
                   number_of_time_steps+list_of_correlation_intervals(end)-1); % Random walk
    random_walk  = random_walk - random_walk(list_of_correlation_intervals(end));
    input_signal = random_walk + white_noise;
    
    %% Estimate Mean Squared Error
    for i = 1:number_of_correlation_intervals
        window_length  = list_of_correlation_intervals(i);
        input_data     = input_signal(end-number_of_time_steps-window_length+2:end);
        moving_average = filter(ones(1,window_length)/window_length,1,input_data);
        moving_average = moving_average(window_length:end);
        
        % MSE with Random Walk as reference input
        estimated_MSE(:,i) = estimated_MSE(:,i) + ...
            (random_walk(end-number_of_time_steps+1:end)-moving_average).^2;
        
        if index_iter==1
            % Calculated MSE with Random Walk as reference input
            calculated_MSE(i) = power_spectral_density/(window_length*sampling_interval) + ...
                                (random_walk_coefficient^2)*window_length*sampling_interval/3;
        end % NOTE: END IF statement 'index_iter==1'
    end % NOTE: END FOR loop 'number_of_correlation_intervals'
end % NOTE: FOR loop 'number_of_iterations'
estimated_MSE = estimated_MSE/number_of_iterations;

%% Plot the results
%%% MSE of MA filter with random walk input corrupted by white noise
% default_color_map = jet(256);
% custom_color_map = default_color_map(1:floor(256/number_of_target_points):256,:);
custom_color_map = [1 0 0; 0 0 1; 0 1 0];
legend_cell      = cell(number_of_target_points+1,1);
jf = java.text.DecimalFormat;   % To add comma to numbers
figure(01)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(list_of_correlation_intervals,calculated_MSE,'Color',[0.7 0.7 0.7],'Linewidth',2.4)
legend_cell{1} = 'Calculated';
for i = 1:number_of_target_points
    plot(list_of_correlation_intervals,estimated_MSE(target_points(i),:)',...
         '--','Color',custom_color_map(i,:),'Linewidth',1.2)
    legend_cell{i+1} = ['Est. $@k = $ ' char(jf.format(target_points(i)-1))];
end % NOTE: FOR loop 'number_of_target_points'
legend(legend_cell,'NumColumns',2,'Location','best','Interpreter','latex','FontSize',13)
set(gca,'xtick',[1e0 1e2 1e4],'XScale','log','YScale','log','FontSize',13)
ylabel('Mean Squared Error $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Window Length $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])
