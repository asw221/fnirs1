
function N = ess(x)
% Univarite effective sample sizes using the method of Geyer
if (ismatrix(x) && ~isvector(x))
    P = size(x, 1);
else
    P = numel(x);
    x = reshape(x, P, 1);
end
rho = fnirs1.utils.acorr(x);
method = 'geyer';
if strcmpi(method, 'geyer')
    sumRho = zeros(1, size(x, 2));
    for j = 1:size(x, 2)
        condition = false;
        M = 0;
        while (~condition && (2 * M + 2) <= P)
            condition = 0 >= (rho(2 * M + 1, j) + rho(2 * M + 2, j));
            M = M + 1;
        end
        M = 2 * max(M - 1, 1) + 2;
        if (M > 1)
            sumRho(j) = sum(rho(2:M, j));
        end
    end
    N = P ./ (1 + 2 * sumRho);
elseif strcmpi(method, 'spectral')
    psd = real(fft(rho));
    if ismatrix(x)
        psd0 = abs(psd(1, :));
    else
        psd0 = abs(psd(1));
    end
    N = P * var(x) ./ psd0;
end
N(isnan(N)) = 0;
N = max(N, 1);
end
