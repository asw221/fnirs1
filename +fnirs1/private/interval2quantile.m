
function Q = interval2quantile(w)
% Convert interval widths to probabilities for quantile bounds
if (~isnumeric(w) || ~isvector(w))
    error('interval2quantile only for numeric vector input');
end
if (size(w, 1) > size(w, 2))
    n = numel(w);
    w = reshape(w, [1 n]);
end
w = max(min(w, 1), 0);
Q = [ (1 - w)/2; 0.5 + w/2 ];
end
