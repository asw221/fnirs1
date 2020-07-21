
function ord = sort_order(A, varargin)
% Return the order of A resulting from a sort
%
% Example usage:
%   A = rand(10, 1);
%   all( A(fnirs1.utils.sort_order(A)) == sort(A) )  % returns true
%
[~, ord] = sort(A, varargin{:});
end
