
function H = head(A, varargin)
% Return the beginning of a vector or matrix-like object, A
H = [];
if ( nargin < 2 )
    n = 6;
else
    n = int64(varargin{1});
end
singletons = size(A) == 1;
if ( sum(singletons) == 1 )
    if ( abs(n) > length(A) )
        n = sign(n) * length(A);
        if ( n < 0 && abs(n) == length(A))
            n = 0;
        end
    end
    if ( n > 0 )
        H = A(1:n);
    elseif ( n < 0 )
        H = A(1:(end + n));
    end
else
    if ( abs(n) > size(A, 1) )
        n = sign(n) * size(A, 1);
        if ( n < 0 && abs(n) == size(A, 1) )
            n = 0;
        end
    end
    if ( n > 0 )
        H = A(1:n, :);
    elseif ( n < 0 )
        H = A(1:(end + n), :);
    end
end
end
