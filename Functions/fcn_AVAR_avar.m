function allan_variance = fcn_AVAR_avar(data,list_of_correlation_intervals,...
                          varargin)
%% fcn_AVAR_avar
%   This function computes allan variance of regularly sampled data over 
%   correlation interval in 'list_of_correlation_intervals'. It uses normal 
%   algorithm.
%
% FORMAT:
%   allan_variance = fcn_AVAR_avar(data,list_of_correlation_intervals)
%
% INPUTS:
%   data: A N x 1 vector of data points.
%   list_of_correlation_intervals: A M x 1 vector containing list of 
%   correlation intervals.
%   varargin: figure number for debugging.
%
% OUTPUTS:
%   allan_variance: A M x 1 vector containing allan variance corresponding 
%   to the correlation intervals.
%
% EXAMPLES:
%   See the script:
%       script_test_fcn_AVAR_avar.m for a full test suite.
%
% This function was written on 2021_05_14 by Satya Prasad
% Questions or comments? szm888@psu.edu
% Updated: 2022/02/15

flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plot  = 0; % Flag to plot the results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1, 'STARTING function: %s, in file: %s\n', st(1).name, st(1).file);
end

%% Check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_check_inputs
    % Are there the right number of inputs?
    if 2>nargin || 3<nargin
        error('Incorrect number of input arguments')
    end
    
    % Check input type and domain
    try
        fcn_AVAR_checkInputsToFunctions(data,'avar data');
    catch ME
        assert(strcmp(ME.message,...
            'The data input must be a N x 1 vector of real numbers'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(list_of_correlation_intervals,'avar interval');
    catch ME
        assert(strcmp(ME.message,...
            'The list_of_correlation_intervals input must be a M x 1 vector of increasing positive integers'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
end

% Does the user want to make a plot at the end?
if 3 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%% Calculate Allan Variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of data points in the INPUT data
number_of_datapoints = numel(data);
% number of interested correlation intervals
number_of_correlation_intervals = numel(list_of_correlation_intervals);

% initialize variable to store allan variance
allan_variance = nan(number_of_correlation_intervals,1);
for i = 1:number_of_correlation_intervals
    % loop over the list of correlation intervals
    
    correlation_interval = list_of_correlation_intervals(i);
    allan_variance_sum = 0; % initialize Allan Variance sum to zero
    for m = 1:(number_of_datapoints-2*correlation_interval)
        allan_variance_sum = allan_variance_sum+...
            (mean(data((m+correlation_interval):(m+2*correlation_interval-1)))...
            -mean(data(m:(m+correlation_interval-1)))).^2;
    end
    % write Allan Variance to the output
    allan_variance(i) = 0.5*allan_variance_sum/...
                        (number_of_datapoints-2*correlation_interval);
end

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plot
    figure_title = inputname(1);
    figure_title('_' == figure_title) = ' ';
    figure(fig_num)
    clf
    width = 540; height = 400; right = 100; bottom = 400;
    set(gcf, 'position', [right, bottom, width, height])
    plot(list_of_correlation_intervals,allan_variance,'b','Linewidth',1.2)
    grid on
    set(gca,'XScale','log','YScale','log','FontSize',13)
    ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
    xlabel('Correlation Interval [Number of Samples]','Interpreter','latex','FontSize',18)
    title(figure_title,'Interpreter','latex','FontSize',18)    
end

if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end