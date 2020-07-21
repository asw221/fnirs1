
function T = center(T)
% Mean center non-categorical columns of a table, T
if ~istable(T)
    error('fnirs1.utils.center - argument should be a table')
end

nms = T.Properties.VariableNames;
for j = 1:numel(nms)
    if ( ~iscategorical(T.(nms{j})) && ~isbinary(T.(nms{j})) )
        T.(nms{j}) = T.(nms{j}) - nanmean(T.(nms{j}));
    end
end
end

function B = isbinary(A)
B = all(A == 0 | A == 1);
end
