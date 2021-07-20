
%% Example group analysis 
% Group analyses require multiple participant data files and some guiding
% information about group level covariates. You can specify this group
% level information using data tables and an associated model formula
%

%% Locate and organize participant data files
% Retrieve data files from multiple directories and read patient
% demographic data
data_files = fnirs1.get_data_files;                % Gather data files from task 1
data_files = [data_files; fnirs1.get_data_files];  % Concatenate data files from task 2
data_files = reshape(data_files, numel(data_files), 1);

demo = readtable('~/Box/BayesianDataAnalysis/EN_CH_MA_PA_Tasks_NEW/Demographicinfo.xlsx');

% Subset the data and add a categorical column 
% demo = demo(demo.ID > 3000, :);
demo = fnirs1.expand_table_conditions(demo, 2, 'Task');

% Check files are in the correct order:
% [ok, mismatch] = fnirs1.check_filenames_id(data_files, demo.ID);

% Center the non-categorical columns:
cdemo = fnirs1.utils.center(demo);


%% Fit the model
% Demographic information: names Cond and TempDeriv are special/reserved
% and do not need to be part of table. In a model formula, put Cond last of
% all categorical variables: e.g. Task + Cond + Task:Cond
%

fit = fnirs1.dlm(data_files, ...
    'GroupData', cdemo, ...
    'GroupFormula', 'ID ~ Task * Cond + LWIDraw + Age', ...
    'DownSampleRate', 10, ...
    'SpecificChannels', 1, ...
    'McmcControl', fnirs1.mcmc_control(1000, 2));


%% Add Contrasts:
% Group Level Parameters are:
% >> fit(1).Descriptions(1:11)
% ans =
%   11×1 cell array
%     {'Population Effect: Cond_1'       }
%     {'Population Effect: LWIDraw'      }
%     {'Population Effect: Age'          }
%     {'Population Effect: Task_2'       }
%     {'Population Effect: Cond_2'       }
%     {'Population Effect: Cond_3'       }
%     {'Population Effect: Task_2:Cond_2'}
%     {'Population Effect: Task_2:Cond_3'}
%     {'Population Effect: TempDeriv_1'  }
%     {'Population Effect: TempDeriv_2'  }
%     {'Population Effect: TempDeriv_3'  }
%
% Condition 2 > Condition 1 is the main effect of Condition 2 in Task_1. 
% Let's add C2 > C1 in Task 2
%

v = [0 0 0 1 1 0 1 0] - [0 0 0 1 0 0 0 0]; 
% v auto-padded with trailing 0's to correct length
fit = add_contrast(fit, v);

% Add additional contrasts for (i) C2 > C1 in Task 1 and (ii) average of 
% {C1, C2} > C3 in Task 1. To do this, setup contrast vectors as
% **columns** of a contrast matrix:
V = [ ([0 0 0 0 1 0] - [1 0 0 0 0 0])', ...
      ([0.5 0 0 0 0.5 0] - [0 0 0 0 0 1])' ];
fit = add_contrast(fit, V);

