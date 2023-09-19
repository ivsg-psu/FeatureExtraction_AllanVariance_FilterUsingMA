function mean_squared_error = ...
    fcn_AVAR_mseMA(power_spectral_density,random_walk_coefficient,...
    window_length,sampling_interval,noise_model)
%% fcn_AVAR_mseMA
%   This function calculates MSE in MA filter.
%
% FORMAT:
%   mean_squared_error = ...
%       fcn_AVAR_mseMA(power_spectral_density,random_walk_coefficient,...
%       window_length,sampling_interval,noise_model)
%
% INPUTS:
%   power_spectral_density: Power spectral density of white noise [unit^2 s].
%   random_walk_coefficient: Noise coefficient for random walk [unit/sqrt(s)].
%   window_length: Window length of MA filter.
%   sampling_interval: Sampling interval [s].
%   noise_model:
%       1-> Random walk corrupted by white noise
%       2-> Constant corrupted by random walk and white noise
%
% OUTPUTS:
%   mean_squared_error: Mean squared error.
%
% EXAMPLES:
%   See the script:
%       script_test_fcn_AVAR_mseMA.m for a full test suite.
%
% This script was written on 2023_09_08 by Satya Prasad
% Questions or comments? szm888@psu.edu

flag_do_debug = 0; % Flag to plot the results for debugging
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
    if 5~=nargin
        error('Incorrect number of input arguments')
    end
    
    % Check input type and domain
    try
        fcn_AVAR_checkInputsToFunctions(power_spectral_density,'positive');
    catch ME
        assert(strcmp(ME.message,...
            'The power_spectral_density input must be a positive number'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(random_walk_coefficient,'positive');
    catch ME
        assert(strcmp(ME.message,...
            'The random_walk_coefficient input must be a positive number'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(window_length,'positive integer');
    catch ME
        assert(strcmp(ME.message,...
            'The window_length input must be a positive integer'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
    try
        fcn_AVAR_checkInputsToFunctions(sampling_interval,'positive');
    catch ME
        assert(strcmp(ME.message,...
            'The sampling_interval input must be a positive number'));
        fprintf(1, '%s\n\n', ME.message)
        return;
    end
end

%% Calculate MSE of MA filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1==noise_model
    mean_squared_error = ...
        power_spectral_density/(window_length*sampling_interval) + ...
        (random_walk_coefficient^2)*window_length*sampling_interval/3;
elseif 2==noise_model
    disp(1);
end % NOTE: END IF statement

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
if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end