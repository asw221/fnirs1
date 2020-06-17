
function B = interactive
% Returns logical(1) if Matlab is being run in an interactive session with
% figure windows enabled

% Best available answer at time of writing. May not cover all cases:
% https://stackoverflow.com/questions/6754430/determine-if-matlab-has-a-display-available
%
B = ~(usejava('jvm') && ~feature('ShowFigureWindows'));
end
