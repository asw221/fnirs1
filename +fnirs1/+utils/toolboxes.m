
function tbx = toolboxes
% Return cellstr of available toolbox names
v = ver;
tbx = cellstr(char(v.Name));
end
