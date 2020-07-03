
function P = explicit_path(pth)
% Return full path string
if (nargin == 0)
    pth = pwd;
end
P = char(pth);
nfs = numel(filesep);
if (isempty(P) || strcmp(P, ['.', filesep]))
    P = pwd;
end
if (numel(P) >= (2 + nfs) && strcmp(P(1:(2 + nfs)), ['..', filesep]))
    if (numel(P) > (2 + nfs))
        P = fullfile(fileparts(pwd), P((3 + nfs):end));
    else
        P = fileparts(pwd);
    end
end
if (numel(P) >= (1 + nfs) && strcmp(P(1:(1 + nfs)), ['.', filesep]))
    if (numel(P) > (1 + nfs))
        P = fullfile(pwd, P((2 + nfs):end));
    else
        P = pwd;
    end
end
if ~strcmp(P(1:nfs), filesep)
    P = fullfile(pwd, P);
end
if (strcmp(P((numel(P) - nfs + 1):end), filesep))
    P = P(1:(numel(P) - nfs));
end
end
