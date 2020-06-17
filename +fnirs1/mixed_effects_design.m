
classdef mixed_effects_design
    % FNIRS1.MIXED_EFFECTS_DESIGN Extract design information from a model
    % formula. Requires the Statistics and Machine Learning Toolbox
    %
    % Example usage:
    %   ID = 1:4;
    %   Var = rand(4, 1);
    %   Grp = categorical([1 1 2 2]);
    %   tbl = table(ID, Var, Grp)
    %   med = fnirs1.mixed_effects_design(tbl, 'Var ~ Grp')
    %   designMatrix(med)
    %
   properties (SetAccess = private) 
       Categories;           % Struct containing categories from any categorical covariates
       CoefficientNames;     % Names of fixed effects coefficients
       FixedEffectsMatrix;   % Fixed effects design matrix
       Formula;              % Formula the object was constructed with
       HasIntercept;         % Logical true if model has a global intercept
       MainEffect;           % Logical (fixed) main effects index
       RandomEffectsMatrix;  % Random effects design matrix
       VariableInfo;         % Table of variable information
   end
   properties (GetAccess = private)
       RandomEffectsNameInfo;  % table of information for constructing random effects names
   end
   methods
       function obj = mixed_effects_design(tbl, formula, varargin)
           % Constructor for fnirs1.mixed_effects_design objects
           %
           % Arguments tbl should be an object of class 'table', and 
           % formula shouldbe a valid model formula character string
           
           % Plan: to create an unfitted LinearMixedModel object and use
           %   the result to extract and save design matrix components
           if (~fnirs1.sml_toolbox_available)
               error('Formula parsing requires the Statistics and Machine Learning Toolbox');
           end
           if (~istable(tbl))
               error('First argument must be of class table');
           end
           
           DontOptimize = statset('LinearMixedModel');
           DontOptimize.MaxIter = 0;
           
           % Turn warnings off before "fitting" with fitlme
           warning('off', 'all');
           dummyLme = fitlme(tbl, formula, 'OptimizerOptions', ...
               DontOptimize, varargin{:});
           warning('on', 'all');
           
           obj.Categories = struct();
           obj.CoefficientNames = dummyLme.CoefficientNames;
           obj.FixedEffectsMatrix = designMatrix(dummyLme, 'Fixed');
           obj.Formula = dummyLme.Formula;
           obj.HasIntercept = contains('(Intercept)', obj.CoefficientNames);
           obj.MainEffect = ~contains(obj.CoefficientNames, ':');
           [~, obj.RandomEffectsNameInfo] = randomEffects(dummyLme);
           obj.RandomEffectsMatrix = designMatrix(dummyLme, 'Random');
           obj.VariableInfo = dummyLme.VariableInfo;
           
           catVars = obj.VariableInfo.Row(obj.VariableInfo.IsCategorical);
           if (~isempty(catVars))
               catCats = cell(size(catVars));
               for i = 1:numel(catVars)
                   catCats{i} = categories(tbl.(catVars{i}));
               end
               obj.Categories = cell2struct(catCats, catVars);
           end
       end
       function D = designMatrix(obj, varargin)
           % Return either the fixed or random effects design matrix
           type = 'fixed';
           if (nargin > 1)
               if ~(isstring(varargin{1}) || ischar(varargin{1}))
                   error('designMatrix: second argument must be character type');
               end
               type = char(varargin{1});
           end
           if (strcmpi(type, 'fixed'))
               D = obj.FixedEffectsMatrix;
           elseif (strcmpi(type, 'random'))
               D = full(obj.RandomEffectsMatrix);
           else
               error('designMatrix: unrecognized option for second argument: %s', ...
                   type);
           end
       end
       function nms = randomEffectsNames(obj)
           % Extract the random effects names from a
           % fnirs1.mixed_effects_design object
           nms = cell(size(obj.RandomEffectsNameInfo, 1), 1);
           for i = 1:numel(nms)
               nms{i} = sprintf('%s:%s %s', ...
                   obj.RandomEffectsNameInfo.Group{i}, ...
                   obj.RandomEffectsNameInfo.Level{i}, ...
                   obj.RandomEffectsNameInfo.Name{i});
           end
       end
       function obj = reformat_base_condition(obj, varargin)
           % Reformats the fixed effect matrix base condition in a model
           % with categorical predictors (so that the base is the first
           % condition in the first appearing categorical main effect).
           % Intended to be kind of a work-around for Matlab's difficulty 
           % supporting ANOVA-type coded designs
           %
           
           relevelCondition = '';
           if (nargin > 1)
               try
                   relevelCondition = char(varargin{1});
               catch ME
                   error(ME.identifier, '%s', ME.message);
               end
           end
           
           X = obj.FixedEffectsMatrix;
           dummyColumns = logical(size(X, 2));
           for j = 1:size(X, 2)
               u = unique(X(:, j));
               dummyColumns(j) = all((u == 1) | (u == 0));
           end
           if ~isempty(relevelCondition)
               dummyColumns = find(dummyColumns & obj.MainEffect & ...
                   (fnirs1.utils.regexpl(obj.CoefficientNames, ...
                   relevelCondition) | ...
                   fnirs1.utils.regexpl(obj.CoefficientNames, ...
                   '(Intercept)')));
           else
               dummyColumns = find(dummyColumns & obj.MainEffect);
           end
           % From above: note that we're not considering any interaction
           % terms
           if isempty(dummyColumns)
               warning('reformat_base_condition:CatNotFound', ...
                   'Categories not found');
           elseif (numel(dummyColumns) > 1)
               if (obj.HasIntercept)
                   baseConditionIndex = dummyColumns(2);
               else
                   baseConditionIndex = dummyColumns(1);
               end
               baseConditionName = strsplit(...
                   obj.CoefficientNames{baseConditionIndex}, '_');
               baseConditionName = sprintf('%s_%s', ...
                   baseConditionName{1}, ...
                   obj.Categories.(baseConditionName{1}){1});
               X(:, 1) = ones(size(X, 1), 1) - ...
                   sum(X(:, dummyColumns(2:end)), 2);
               obj.CoefficientNames{1} = baseConditionName;
               obj.FixedEffectsMatrix = X;
           end
       end
   end
end
