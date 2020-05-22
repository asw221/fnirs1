
function etbl = expand_table_conditions(tbl, nConditions, varargin)
% FNIRS1.EXPAND_TABLE_CONDITIONS  Add a categorical column to a table and
% duplicate rows for each category
%
% E = FNIRS1.EXPAND_TABLE_CONDITIONS(tbl, n) appends a categorical column
% named 'Condition' to tbl with n categories. Duplicates the rows of tbl n
% times for each category in 'Condition' so that the result, E, has
% size(tbl, 1) * n rows and size(tbl, 2) + 1 columns.
% 
% E = FNIRS1.EXPAND_TABLE_CONDITIONS(tbl, n, 'ConditionColumnName')
% performs the same function, but allows for the appended categorical
% column to have a different name.
%
% Examples:
%   ID = [1:4]'
%   Var = rand(size(ID))
%   tbl = table(ID, Var)
%   E = fnirs1.expand_table_conditions(tbl, 2)
%
% See also
% table

conditionName = 'Condition';
if ~istable(tbl)
    error('fnirs1.expand_table_conditions: first argument must be a table');
end
if ~(isnumeric(nConditions) && isscalar(nConditions))
    error('fnirs1.expand_table_conditions: nConditions must be a numeric scalar');
else
    nConditions = floor(abs(nConditions));
end
if (nargin > 2)
    if (ischar(varargin{1}) || isstring(varargin{1}))
        conditionName = varargin{1};
    else
        error('fnirs1.expand_table_conditions: condition name argument must be character');
    end    
end
N = size(tbl, 1);
etblIndex = reshape(repmat(1:N, nConditions, 1), nConditions * N, 1);
newCondition = repmat((1:nConditions)', N, 1);
etbl = [tbl(etblIndex, :), ...
    table(categorical(newCondition), 'VariableNames', {conditionName})];
end
