
function B = sml_toolbox_available
% Query if the Statistics and Machine Learning Toolbox is available for the
% current Matlab license. Return true/false
B = any(strcmp(toolboxes, 'Statistics and Machine Learning Toolbox'));
end
