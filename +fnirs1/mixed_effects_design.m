
classdef mixed_effects_design
   properties (SetAccess = private) 
       FixedEffectsMatrix;   % Fixed effects design matrix
       RandomEffectsMatrix;  % Random effects design matrix
   end
   methods
       function obj = mixed_effects_design(tbl, formula)
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
           
           sset = statset('LinearMixedModel');
           sset.MaxIter = 0;
           
           % Turn warnings off before "fitting" with fitlme
           warning('off', 'all');
           dummyLme = fitlme(tbl, formula, 'OptimizerOptions', sset);
           warning('on', 'all');
           
           obj.FixedEffectsMatrix = designMatrix(dummyLme, 'Fixed');
           obj.RandomEffectsMatrix = designMatrix(dummyLme, 'Random');
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
   end
end
