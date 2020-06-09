
function p = home
% Returns the fnirs1 install path
persistent fp_;
if isempty(fp_)
    fp_ = fileparts(mfilename('fullpath'));
end
p = fp_;
end
