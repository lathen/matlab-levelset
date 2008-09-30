function C = plot(ls, varargin)

if isempty(varargin)
    [C,H] = contour(ls.phi,[0 0],'r');
else
    hold on;

    if isenabled('phi',varargin{:})
        imagesc(ls.phi, [min(ls.phi(:)), max(ls.phi(:))]);
        colormap(gray);
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
            contour(ls.phi,[0 0],'r');
        else
            contour(ls.phi,[0 0],param);
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