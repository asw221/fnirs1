
function B = is_column_categorical(tbl)
if ~istable(tbl)
    error('fnirs1.is_column_categorical: input must be a table');
end
M = size(tbl, 2);
B = false(1, M);
colnms = tbl.Properties.VariableNames;
for j = 1:numel(colnms)
    B(j) = iscategorical(tbl.(colnms{j}));
end
end
