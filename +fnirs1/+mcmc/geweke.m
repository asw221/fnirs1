
function G = geweke(x, varargin)
% Compute Geweke statistic to test for convergence of an MCMC chain. Values
% close to zero indicate approximate convergence.
%
% Details:
% By default, conducts an equal variance t-test for equality of means from
% the first 10% and last 50% of the MCMC chain. If the test fails to reject
% (returns statistics < ~2) then we conclude the chain is approximately
% stationary and convergent.
%
% fnirs1.mcmc.geweke(x)
%   Computes the Geweke statistic with default 10/50 proportions
%
% fnirs1.mcmc.geweke(x, 0.2)
%   Uses the first 20% of the and last 50% of the chain to compute the
%   Geweke statistic
%
% fnirs1.mcmc.geweke(x, 0.2, 0.4)
%   Uses the first 20% and last 40% of the chain to compute the Geweke
%   statistic. These proportions can be any number on (0, 1) as long as
%   they do not add to >= 1
%
%
% Example usage:
%   rng(8675309);
%   x = rand(100, 1);           % A stationary MCMC chain
%   y = [x; rand(100, 1) + 5];  % A non-stationary MCMC chain
%   fnirs1.mcmc.geweke(x)       % Returns 0.2802 (do not reject)
%   fnirs1.mcmc.geweke(y)       % Returns 7.5721 (reject)
% 
% 
% See also:
% fnirs1.mcmc.ess, fnirs1.mcmc.geweke_plot, fnirs1.mcmc.geweke_sequence
%
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
            error('geweke:BadProbability', ...
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
            error('geweke:BadProbability', ...
                'geweke optional arguments should be on (0, 1)');
        end
        if (p(2) > (1 - p(1)))
            error('geweke:NonIndependentSubsets', ...
                'geweke proportion arguments should sum to < 1');
        end
    catch ME
        error(ME.identifier, '%s', ME.message);
    end
end

N = size(x, 1);
offset = 0;
if (nargin > 3)
    try
        offset = abs(fix(double(varargin{3}(1))));
        if (offset > diff(fix(N * p) - 1))
            error('geweke:OffsetTooLarge', ...
                'geweke start offset will result in non-independent samples');
        end
    catch ME
        error(ME.identifier, '%s', ME.message);
    end
end

n0 = [1 + offset, fix(N * p(2))];
n1 = [fix(N * p(1)) + offset, N];

x0 = x(n0(1):n1(1), :);
x1 = x(n0(2):n1(2), :);
if (size(x0, 1) <= 1 || size(x1, 1) <= 1)
    error('geweke:InsufficientSamples', ...
        'Insufficient samples to estimate Geweke diagnostic');
end

v0 = var(x0) * (size(x0, 1) - 1) ./ fnirs1.mcmc.ess(x0);
v1 = var(x1) * (size(x1, 1) - 1) ./ fnirs1.mcmc.ess(x1);
% z = NaN(1, size(x, 2));
G = abs(mean(x1) - mean(x0))  ./ sqrt(v1 + v0);
end
