
function AC = acorr(x, varargin)
% Compute the linear auto-correlation of a univariate signal using the fast
% fourier transform. Applied column-wise for matrices
if ~isnumeric(x)
    error('acorr only defined for numeric inputs');
end

trim = true;
if (nargin > 1)
    try
        trim = logical(varargin{1}(1));
    catch ME
        warning(ME.identifier, '%s', ME.message);
        error('acov: trim argument should convert to logical');
    end
end

AC = fnirs1.utils.acov(x, trim);
if (ismatrix(x) && ~isvector(x))
    AC = (AC' ./ AC(1, :)')';
else
    AC = AC / AC(1);
end
end
