
classdef mcmc_control
    % FNIRS1.MCMC_CONTROL Organizes parameters that control the MCMC
    %   Used to fine tune MCMC implementation. For example, users can
    %   control how long the MCMC chains are run (longer chains can
    %   improve accuracy, but will also lengthen run-time), or whether or
    %   not to include the HRF temporal derivatives in model output.
    %
    % Example usage:
    %   maxIterations = 1e4;        % default = 20,000
    %   includeDerivatives = true;  % default = false
    %   expectedKnots = 12;         % default = 15
    %   burnin = 4e3;               % default = (maxIterations / 2)
    %   ctrl = fnirs1.mcmc_control(maxIterations, includeDerivataves, expectedKnots, burnin);
    %   disp(ctrl.maxIterations)
    %
    %   Note that default values exist for all of the above parameters, so
    %   that the user need not specify every parameter. Parameters can only
    %   be set by the object constructor.
    %
properties (SetAccess = private)
    burnin;              % Number of burnin iterations for MCMC chain (default = maxIterations / 2)
    expectedKnots;       % Expected number of spline-basis knots to model HRF (default = 15)
    includeDerivatives;  % Logical flag to include temporal derivatives (default = false)
    maxIterations;       % Maximum number of MCMC iterations (default = 20,000)
end
methods
    function obj = mcmc_control(varargin)
        % Construct an mcmc_control object; ensure any inputs are valid
        
        % defaults
        obj.maxIterations = 2e4;
        obj.expectedKnots = 15;
        obj.includeDerivatives = false;
        if (nargin >= 1)
            if (isnumeric(varargin{1}) && isscalar(varargin{1}) && varargin{1} > 0)
                obj.maxIterations = fix(varargin{1});
            else
                error("MCMC max iterations must be a positive integer");
            end
        end
        % burnin default: maxIterations / 2
        obj.burnin = fix(obj.maxIterations / 2);
        if (nargin >= 2)
            if (isscalar(varargin{2}))
                obj.includeDerivatives = logical(varargin{2});
            else
                error("MCMC include temporal derivatives must convert to logical");
            end
        end
        if (nargin >= 3)
            if (isnumeric(varargin{3}) && isscalar(varargin{3}) && varargin{3} > 0)
                obj.expectedKnots = fix(varargin{3});
            else
                error("MCMC expected number of knots must be a positive integer");
            end
        end
        if (nargin >= 4)
            if (isnumeric(varargin{4}) && isscalar(varargin{4}) && varargin{4} > 0)
                obj.burnin = fix(varargin{4});
                if (obj.burnin >= obj.maxIterations)
                    error("MCMC burnin should be between [0, maxIterations)");
                end
            else
                error("MCMC burnin must be a positive integer");
            end
        end
        if (nargin > 4)
            warning("MCMC - unused inputs");
        end
    end
end
end
