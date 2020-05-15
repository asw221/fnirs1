
classdef mixed_effects_design
   properties (SetAccess = private) 
       Categories;           %
       CoefficientNames;     % Names of fixed effects coefficients
       FixedEffectsMatrix;   % Fixed effects design matrix
       Formula;              % Formula the object was constructed with
       MainEffect;           % Logical (fixed) main effects index
       RandomEffectsMatrix;  % Random effects design matrix
       VariableInfo;         % Table of variable information
   end
   properties (GetAccess = private)
       RandomEffectsNameInfo;
   end
   methods
       function obj = mixed_effects_design(tbl, formula, varargin)
           % Constructor for fnirs1.mixed_effects_design objects. Arguments
           % tbl should be an object of class 'table', and formula should
           % be a valid model formula character string. Requires the
           % Statistics and Machine Learning Toolbox
           
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
           nms = cell(size(obj.RandomEffectsNameInfo, 1), 1);
           for i = 1:numel(nms)
               nms{i} = sprintf('%s:%s %s', ...
                   obj.RandomEffectsNameInfo.Group{i}, ...
                   obj.RandomEffectsNameInfo.Level{i}, ...
                   obj.RandomEffectsNameInfo.Name{i});
           end
       end
       function obj = reformat_base_condition(obj)
           X = obj.FixedEffectsMatrix;
           dummyColumns = logical(size(X, 2));
           for j = 1:size(X, 2)
               u = unique(X(:, j));
               dummyColumns(j) = all((u == 1) | (u == 0));
           end
           dummyColumns = find(dummyColumns & obj.MainEffect);
           if (numel(dummyColumns) > 1)
               baseConditionIndex = dummyColumns(2);
               baseConditionName = strsplit(...
                   obj.CoefficientNames{baseConditionIndex}, '_');
               baseConditionName = sprintf('%s_%s', ...
                   baseConditionName{1}, ...
                   obj.Categories.(baseConditionName{1}){1});
               X(:, 1) = X(:, 1) - X(:, baseConditionIndex);
               obj.CoefficientNames{1} = baseConditionName;
               obj.FixedEffectsMatrix = X;
           end
       end
   end
end
