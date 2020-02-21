
function coreInfo = core_info()
% FNIRS1.CORE_INFO returns a struct with information about the number of
% physical and logical CPUs available to the MATLAB session
%
persistent CoreInfo;
if (isempty(CoreInfo))
    info = evalc('feature(''numcores'')');
    cl = strsplit(info, ' ');
    cores = nan(4, 1);
    count = 0;
    for i = 1:length(cl)
        n = sscanf(cl{i}, '%i');
        if ~isempty(n)
            count = count + 1;
            if (count <= length(cores))
                cores(count) = n;
            end
        end
    end
    if (count == 0)
        error('Could not detect core info');
    end
    CoreInfo = struct(...
        'physical', cores(1), ...
        'logical', cores(2), ...
        'assigned', cores(3), ...
        'using', cores(4));
end
coreInfo = CoreInfo;
end
