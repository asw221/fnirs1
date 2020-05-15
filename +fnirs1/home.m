
function p = home
% Return the fnirs1 install path
info = what('fnirs1');
if ~isempty(info)
    p = info.path;
else
    p = '';
end
end
