
function AC = acov(x, varargin)
% Compute the auto-covariance of a univariate signal using the fast fourier
% transform. Applied column-wise for matrices
if ~isnumeric(x)
    error('acov only defined for numeric inputs');
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

if (ismatrix(x) && ~isvector(x))
    N = size(x, 1);
else
    N = numel(x);
end

NN = 2^ceil(log2(N));
cx = complex((x' - mean(x)')');
z = fft(cx, NN);
AC = real(ifft(z .* conj(z)));
if (trim)
    if (ismatrix(x) && ~isvector(x))
        AC = AC(1:N, :);
    else
        AC = AC(1:N);
    end
end
end

%     N <- length(x)
%     NN <- 2^ceiling(log2(N))
%     x0 <- c(x - mean(x, na.rm = TRUE), rep(0, NN - N))
%     z <- fft(x0)
%     ac <- Re(fft(z * Conj(z), inverse = TRUE))/NN
%     if (.trim) 
%         ac[1:N]
%     else ac
