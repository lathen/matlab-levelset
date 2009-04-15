function [C,h] = plot(ls, varargin)
% PLOT  Plots the zero level set and various attributes of a levelset2D
%       object.
%   [C,h] = PLOT(LS, ...) plots the zero level set of LS. Optional output
%   argument C contains points along the contour, and h gives a handle to
%   the graphics object. A number of optional input arguments can be given:
%       'phi'           displays the level set function
%       'narrowband'    plots the narrowband
%       'contour'       plots the zero level set (default if no other
%                       parameters are given)
%       'gradient'      plots the gradient of the level set function

%   Author: Gunnar Läthén (gunnar.lathen@itn.liu.se)
%   $Date: 2007/10/17

if isempty(varargin)
    [C,h] = contour(ls.phi,[0 0],'r');
else
    hold on;

    if isenabled('phi',varargin{:})
        max_abs = max(abs(min(ls.phi(:))), abs(max(ls.phi(:))));
        imagesc(ls.phi, [-max_abs, max_abs]);
        
        % Create colormap (green outside, red inside)
        map_size = 256;
        R = 1 :-1/(map_size/2-1) : 0;
        G = 0 : 1/(map_size/2-1) : 1;
        map = zeros(map_size,3);
        map(1:map_size/2,1) = R;
        map(map_size/2+1:map_size,2) = G;
        colormap(map);
    end

    if isenabled('narrowband',varargin{:})
        param = getparameters('narrowband', varargin{:});
        delta = 1;
        x = str2double(param);
        if ~isnan(x)
            delta = x;
        end
        [X,Y] = ind2sub(size(ls), narrowband(ls));
        plot(Y(1:delta:end),X(1:delta:end),'g.');
    end
    
    if isenabled('contour',varargin{:})
        param = getparameters('contour', varargin{:});
        if isempty(param)
            [C,h] = contour(ls.phi,[0 0],'r');
        else
            [C,h] = contour(ls.phi,[0 0],param);
        end
    end

    if isenabled('gradient',varargin{:})
        param = getparameters('gradient', varargin{:});
        delta = 1;
        x = str2double(param);
        if ~isnan(x)
            delta = x;
        end            
        ind = uint32(1:numel(ls.phi));
        [Dx,Dy] = ls.diff_central(ls.phi, ind);
        [rows,cols] = size(ls.phi);
        Dx = reshape(Dx,rows,cols)*delta;
        Dy = reshape(Dy,rows,cols)*delta;
        quiver(1:delta:cols,1:delta:rows,Dx(1:delta:end,1:delta:end),Dy(1:delta:end,1:delta:end), 0, 'b');
    end
    
end

axis ij;
axis image;
[rows, cols] = size(ls.phi);
axis([0 cols 0 rows]);
end


function tf = isenabled(mode, varargin)

for i = 1:length(varargin)
    if strfind(varargin{i}, mode)
        tf = true;
        return;
    end
end
    
tf = false;
end

function param = getparameters(mode, varargin)

for i = 1:length(varargin)
    if strfind(varargin{i}, mode)
        param = strrep(varargin{i}, mode, '');
        param = strtrim(param);
        return;
    end
end
    
param = [];
end
