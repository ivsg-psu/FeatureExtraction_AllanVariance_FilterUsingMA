%%%%%%%%%%%%%%%%%%%% script_AVAR_analyzeRWinMAfilter.m %%%%%%%%%%%%%%%%%%%%
%% Purpose:
%   The purpose of this script is to analyze Moving Average filter using
%   Mean Squared Error and AVAR of error for Random walk input corrupted 
%   by White noise.
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
number_of_iterations = 100;
ma_noise_model       = 1;

sampling_frequency   = 50; % [Hz]
sampling_interval    = 1/sampling_frequency; % [second]
number_of_time_steps = 2^19;

p = floor(log2(number_of_time_steps));
list_of_correlation_intervals   = 2.^(0:p-3)'; % List of correlation intervals
number_of_correlation_intervals = numel(list_of_correlation_intervals);

% Noise parameters
power_spectral_density  = 0.0004; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)]

%%% Initialize variables to store MSE and AVAR
estimated_MSE   = zeros(number_of_time_steps,number_of_correlation_intervals);
calculated_MSE  = NaN(number_of_time_steps,number_of_correlation_intervals);
estimated_AVAR  = NaN(number_of_correlation_intervals);
calculated_AVAR = NaN(number_of_correlation_intervals);
for index_iter = 1:number_of_iterations
    %% Synthesize the test signal
    if 1==index_iter
        time_vector  = sampling_interval*(0:number_of_time_steps-1)';
    end % NOTE: END IF statement '1==index_iter'
    white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,sampling_frequency,...
                   number_of_time_steps+list_of_correlation_intervals(end)-1); % White noise
    random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,sampling_frequency,...
                   number_of_time_steps+list_of_correlation_intervals(end)-1); % Random walk
    random_walk  = random_walk - random_walk(list_of_correlation_intervals(end));
    input_signal = random_walk + white_noise;
    
    %% Estimate Mean Squared Error of MA filter
    for i = 1:number_of_correlation_intervals
        window_length  = list_of_correlation_intervals(i);
        input_data     = input_signal(end-number_of_time_steps-window_length+2:end);
        moving_average = filter(ones(1,window_length)/window_length,1,input_data);
        moving_average = moving_average(window_length:end);
        actual_error   = random_walk(end-number_of_time_steps+1:end)-moving_average;
        
        %% Estimate MSE of MA filter with Random Walk input corrupted by white noise
        estimated_MSE(:,i) = estimated_MSE(:,i) + actual_error.^2;
        if index_iter==number_of_iterations && i==number_of_correlation_intervals
            estimated_MSE = estimated_MSE/number_of_iterations;
        end % NOTE: END IF statement 'index_iter==number_of_iterations && i==number_of_correlation_intervals'
        
        if 1==index_iter
            %% Calculate MSE of MA filter with Random Walk input corrupted by white noise
            calculated_MSE(:,i) = ...
                fcn_AVAR_mseMA(power_spectral_density,random_walk_coefficient,...
                window_length,sampling_interval,ma_noise_model);
            
            %% Estimate AVAR of MA filter error with Random Walk input corrupted by white noise
            estimated_AVAR(:,i) = fcn_AVAR_favar([actual_error; 0],list_of_correlation_intervals);
            
            %% Calculate AVAR of MA filter error with Random Walk input corrupted by white noise
            calculated_AVAR(:,i) = ...
                fcn_AVAR_avarMA(power_spectral_density,random_walk_coefficient,...
                list_of_correlation_intervals,window_length,sampling_interval,ma_noise_model);
        end % NOTE: END IF statement '1==index_iter'
    end % NOTE: END FOR loop 'number_of_correlation_intervals'
end % NOTE: FOR loop 'number_of_iterations'

%% Plot the results
%%% Plot to show AVAR of error in MA filter with Random Walk input
%%% corrupted by White noise
default_color_map = jet(256);
custom_color_map  = default_color_map(1:floor(256/number_of_correlation_intervals):256,:);
legend_cell       = cell(number_of_correlation_intervals,1);
jf = java.text.DecimalFormat;   % To add comma to numbers
figure(01)
clf
width = 1056.2+25; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
subplot(1,2,1)
axis_position = [75/width, 0.1567, 415.6/width, 0.7683];
hold on
grid on
for i = 1:number_of_correlation_intervals
    if i~=7
        plot(list_of_correlation_intervals,estimated_AVAR(:,i),...
             'Color',custom_color_map(i,:),'Linewidth',1.2)
    else
        plot(list_of_correlation_intervals,estimated_AVAR(:,i),...
             'k--','Linewidth',2.4)
    end % NOTE: END IF statement
    legend_cell{i} = ['$M =$ ' char(jf.format(list_of_correlation_intervals(i)))];
end % NOTE: END FOR loop 'number_of_correlation_intervals'
legend(legend_cell,'NumColumns',3,'Location','best','Interpreter','latex','FontSize',13)
annotation('textbox',[0.2,0.2,0.1,0.1],'String',"64",...
           'EdgeColor','w','Interpreter','latex','FontSize',30)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
title('$(a)$','Interpreter','latex','FontSize',18)
ylim([10^-6.5 1e2])

%%% Area of AVAR
axis_position = [(75+100+415.584)/width, 0.1567, 415.6/width, 0.7683];
subplot(1,2,2)
hold on
grid on
plot(list_of_correlation_intervals,...
     sum(0.5*(estimated_AVAR(1:end-1,:)+estimated_AVAR(2:end,:)),1)','k','Linewidth',1.2)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Area of AVAR $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Moving Window $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])

%%% Area of AVAR of error and MSE in MA filter with random walk input corrupted by white noise
figure(02)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(list_of_correlation_intervals,...
     estimated_MSE(number_of_time_steps/2,:)','k--','Linewidth',1.2)
plot(list_of_correlation_intervals,...
     sum(0.5*(estimated_AVAR(1:end-1,:)+estimated_AVAR(2:end,:)),1)','k','Linewidth',1.2)
legend('MSE: 100 iterations','Area of AVAR (error): 1 iteration',...
       'Location','best','Interpreter','latex','FontSize',13)
set(gca,'xtick',[1e0 1e2 1e4],'XScale','log','YScale','log','FontSize',13)
ylabel('Magnitude $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Moving Window $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])

%%% Area of AVAR of error in MA filter with random walk input corrupted by white noise
figure(03)
clf
width = 1056.2+25; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
subplot(1,2,1)
axis_position = [75/width, 0.1567, 415.6/width, 0.7683];
hold on
grid on
xline(2,'--','Color',[0.8500 0.3250 0.0980],'Linewidth',1.2)
xline(64,'b','Linewidth',1.2)
xline(2048,'m-.','Linewidth',1.2)
plot(list_of_correlation_intervals,...
     sum(0.5*(estimated_AVAR(1:end-2,:)+estimated_AVAR(2:end-1,:)),1)','k','Linewidth',1.2)
legend('$M=2$','$M=64$','$M=2$,$048$','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Area of AVAR $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Moving Window $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
title('$(a)$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])

%%% Time domain demo
axis_position = [(75+100+415.584)/width, 0.1567, 415.6/width, 0.7683];
subplot(1,2,2)
hold on
grid on
plot(-1,-1,'Color',[0.7 0.7 0.7],'Linewidth',3)
plot(-1,-1,'.','Color',[0.8500 0.3250 0.0980],'Markersize',10)
plot(-1,-1,'b','Linewidth',1.2)
plot(-1,-1,'m','Linewidth',1.2)

moving_average = filter(ones(1,2)/2,1,input_signal);
moving_average = moving_average(end-number_of_time_steps+1:end);
plot(time_vector,moving_average,'.','Color',[0.8500 0.3250 0.0980],'Markersize',4)

plot(time_vector,random_walk(end-number_of_time_steps+1:end),...
     'Color',[0.7 0.7 0.7],'Linewidth',3)

moving_average = filter(ones(1,64)/64,1,input_signal);
moving_average = moving_average(end-number_of_time_steps+1:end);
plot(time_vector,moving_average,'b','Linewidth',1.2)

moving_average = filter(ones(1,2048)/2048,1,input_signal);
moving_average = moving_average(end-number_of_time_steps+1:end);
plot(time_vector,moving_average,'m','Linewidth',1.2)

legend('Reference Input','$M=2$','$M=64$','$M=2$,$048$','NumColumns',2,...
       'Location','best','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'FontSize',13)
ylabel('Amplitude $[Unit]$','Interpreter','latex','FontSize',18)
xlabel('Time $[s]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
xlim([0 100])

%% Plots for ACC
%%% Plot to show AVAR of error in MA filter with Random Walk input
%%% corrupted by White noise
default_color_map = jet(256);
custom_color_map  = default_color_map(1:floor(256/number_of_correlation_intervals):256,:);
legend_cell       = cell(number_of_correlation_intervals,1);
jf = java.text.DecimalFormat;   % To add comma to numbers
figure(11)
clf
width = 540; height = 770.04+25; right = 100; bottom = 10;
set(gcf, 'position', [right, bottom, width, height])
subplot(2,1,1)
axis_position = [0.1354, (432.72+25)/height, 0.7696, 307.32/height];
hold on
grid on
for i = 1:number_of_correlation_intervals
    if i~=7
        plot(list_of_correlation_intervals,estimated_AVAR(:,i),...
             'Color',custom_color_map(i,:),'Linewidth',1.2)
    else
        plot(list_of_correlation_intervals,estimated_AVAR(:,i),...
             'k--','Linewidth',2.4)
    end % NOTE: END IF statement
    legend_cell{i} = ['$M =$ ' char(jf.format(list_of_correlation_intervals(i)))];
end % NOTE: END FOR loop 'number_of_correlation_intervals'
legend(legend_cell,'NumColumns',3,'Location','best','Interpreter','latex','FontSize',13)
annotation('textbox',[0.35,0.55,0.1,0.1],'String',"64",...
           'EdgeColor','w','Interpreter','latex','FontSize',30)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
title('$(a)$','Interpreter','latex','FontSize',18)
ylim([10^-6.5 1e2])

%%% Area of AVAR
axis_position = [0.1354, 62.7/height, 0.7696, 307.32/height];
subplot(2,1,2)
hold on
grid on
plot(list_of_correlation_intervals,...
     sum(0.5*(estimated_AVAR(1:end-1,:)+estimated_AVAR(2:end,:)),1)','k','Linewidth',1.2)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Area of AVAR $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Moving Window $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])

%%% Area of AVAR of error in MA filter with random walk input corrupted by white noise
figure(13)
clf
width = 540; height = 770.04+25; right = 100; bottom = 10;
set(gcf, 'position', [right, bottom, width, height])
subplot(2,1,1)
axis_position = [0.1354, (432.72+25)/height, 0.7696, 307.32/height];
hold on
grid on
xline(2,'c--','Linewidth',2.4)
xline(64,'b','Linewidth',1.2)
xline(2048,'m-.','Linewidth',2.4)
plot(list_of_correlation_intervals,...
     sum(0.5*(estimated_AVAR(1:end-2,:)+estimated_AVAR(2:end-1,:)),1)','k','Linewidth',1.2)
legend('$M=2$','$M=64$','$M=2$,$048$','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'xtick',[1e0 1e2 1e4],'XScale','log',...
    'YScale','log','FontSize',13)
ylabel('Area of AVAR $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Moving Window $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
title('$(a)$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])

%%% Time domain demo
axis_position = [0.1354, 62.7/height, 0.7696, 307.32/height];
subplot(2,1,2)
hold on
grid on
plot(-1,-1,'Color',[0.7 0.7 0.7],'Linewidth',4)
plot(-1,-1,'c.','Markersize',10)
plot(-1,-1,'b','Linewidth',1.2)
plot(-1,-1,'m','Linewidth',1.2)

plot(time_vector,random_walk(end-number_of_time_steps+1:end),...
     'Color',[0.7 0.7 0.7],'Linewidth',4)

moving_average = filter(ones(1,2)/2,1,input_signal);
moving_average = moving_average(end-number_of_time_steps+1:end);
plot(time_vector,moving_average,'c.','Markersize',4)

moving_average = filter(ones(1,64)/64,1,input_signal);
moving_average = moving_average(end-number_of_time_steps+1:end);
plot(time_vector,moving_average,'b','Linewidth',1.2)

moving_average = filter(ones(1,2048)/2048,1,input_signal);
moving_average = moving_average(end-number_of_time_steps+1:end);
plot(time_vector,moving_average,'m','Linewidth',1.2)

legend('Reference Input','$M=2$','$M=64$','$M=2$,$048$','NumColumns',2,...
       'Location','best','Interpreter','latex','FontSize',13)
set(gca,'Position',axis_position,'FontSize',13)
ylabel('Amplitude $[Unit]$','Interpreter','latex','FontSize',18)
xlabel('Time $[s]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
xlim([0 100])
