
function path = last_path(varargin)
persistent lp_;
if isempty(lp_)
    lp_ = './';
end
if (nargin == 1)
    lp_ = char(varargin{1});
elseif (nargin > 0)
    error('Too many input arguments');
end
path = lp_;
end
