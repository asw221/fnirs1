
function G = geweke_sequence(x, varargin)
% Compute Geweke statistic sequences for increasing burnin lengths based on
% the available samples
if ~isnumeric(x)
    error('geweke diagnostic designed for numeric inputs');
end
if isvector(x)
    x = reshape(x, numel(x), 1);
end

p = [0.1, 0.5];
if (nargin > 1)
    try
        p(1) = double(varargin{1}(1));
        if (p(1) <= 0 || p(1) >= 1)
            error('geweke_sequence:BadProbability', ...
                'geweke proportion arguments should be on (0, 1)')
        end
    catch ME
        error(ME.identifier, '%s', ME.message);
    end
end
if (nargin > 2)
    try
        p(2) = double(varargin{2}(1));
        if (p(2) <= 0 || p(2) >= 1)
            error('geweke_sequence:BadProbability', ...
                'geweke optional arguments should be on (0, 1)');
        end
        if (p(2) > (1 - p(1)))
            error('geweke_sequence:NonIndependentSubsets', ...
                'geweke proportion arguments should sum to < 1');
        end
    catch ME
        error(ME.identifier, '%s', ME.message);
    end
end

N = size(x, 1);
Imax = diff(fix(N * p)) - 1;
step = 10;
if (step > Imax)
    step = 1;
end
A = NaN(length(1:step:Imax), size(x, 2) + 1);
row = 1;
for i = 1:step:Imax
    A(row, :) = [i - 1, fnirs1.mcmc.geweke(x, p(1), p(2), i - 1)];
    row = row + 1;
end
G = array2table(A, 'VariableNames', ...
    cellstr([ "Offset", "Var" + string(1:size(x, 2)) ]));
end
