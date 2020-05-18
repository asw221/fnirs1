
% ETBL = FNIRS1.EXPAND_TABLE_CONDITIONS(TBL, N);
% ETBL = FNIRS1.EXPAND_TABLE_CONDITIONS(TBL, N, 'ConditionColumnName');
function etbl = expand_table_conditions(tbl, nConditions, varargin)
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
