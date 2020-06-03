
classdef contrast
    properties
        Names;
        Vectors;
    end
    methods
        function obj = contrast(vec, varargin)
            obj.Vectors = vec;
            if (nargin > 1)
                obj.Names = varargin{1};
                if (numel(varargin{1}) < size(obj.Vectors, 1))
                    warning('Contrast names shorter than contrast vector');
                end
            end
        end
        function c   = char(obj)
            % Convert contrast to formatted character string
            c = char(string(obj));
        end
        function cs  = cellstr(obj)
            % Convert contrast to formatted cellstr
            cs = cellstr(string(obj));
        end
        function obj = set.Names(obj, names)
            % Set method for Names property. Names must be string-like
            % (convertable to cellstr)
            try
                obj.Names = cellstr(names);
            catch ME
                error('contrast ''names'' must be convertable to cellstr');
            end
            if (size(obj.Names, 2) > size(obj.Names, 1))
                obj.Names = obj.Names';
            end
        end
        function obj = set.Vectors(obj, vec)
            if ~(isnumeric(vec))
                error('contrast vector should be of numeric type');
            end
            obj.Vectors = vec;
            if (isvector(vec) && ...
                    size(obj.Vectors, 2) > size(obj.Vectors, 1))
                obj.Vectors = obj.Vectors';
            end
        end
        function s = string(obj)
            % Convert contrast to formatted string
            nv = size(obj.Vectors, 1);
            nn = numel(obj.Names);
            if (nn < nv)
                nms = cellstr("X" + (1:(nv - nn)));
                nms(1:nn) = obj.Names;
                nms = string(nms);
            else
                nms = string(obj.Names(1:nv));
            end
            s = repmat("", size(obj.Vectors, 2), 1);
            for j = 1:size(obj.Vectors, 2)
                w = obj.Vectors(logical(obj.Vectors(:, j)), j);
                connms = nms(logical(obj.Vectors(:, j)));
                [w, ord] = sort(w, 'descend');
                connms = connms(ord);
                for i = 1:numel(w)
                    if (abs(w(i)) ~= 1)
                        connms(i) = string(abs(w(i))) + "*" + connms(i);
                    end
                    if (i > 1)
                        if ((w(i) < 0)  &&  (w(i - 1) > 0))
                            connms(i) = " > " + connms(i);
                        elseif ((w(i) < 0  &&  w(i - 1) < 0) || ...
                                (w(i) > 0  &&  w(i - 1) > 0))
                            connms(i) = " + " + connms(i);
                        end
                    end
                end
                s(j) = strjoin(connms, "");
            end
        end
    end
end

