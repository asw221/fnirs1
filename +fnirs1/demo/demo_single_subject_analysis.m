
%% Example single subject analysis 
% Illustrate two ways of doing analyses: typical route and alternative
% modular structure which may be useful in some situations
%

%% Fit subject level model: typical route
%
fit = fnirs1.dlm('', ...
    'DownSampleRate', 10, ...
    'SpecificChannels', 9:10, ...
    'McmcControl', fnirs1.mcmc_control(500, true));
                    % (McmcControl controls ^^ TempDeriv specifier)

%% Alternative route: modular structure
% This may be useful for organizing analyses and running them later, or for
% manually tuning a model hyperparameter (like the expected number of
% knots; set via 'McmcControl')
%
setupdir = fnirs1.specify_model('', ...
    'DownSampleRate', 10, ...
    'SpecificChannels', 9:10, ...
    'McmcControl', fnirs1.mcmc_control(500));

setups = fnirs1.list_setup_files(setupdir);

% Single subject model - modular construction of setup file
% This fits the exact same model as above, but gives more of an idea about
% how fnirs1 organizes execution step-by-step
fit = fnirs1.dlm(setups);

