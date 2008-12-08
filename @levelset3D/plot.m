function plot(ls, varargin)

if isempty(varargin)
    %isosurface(ls.phi,0);
    [F,V,N] = triangulate(ls.phi, ls.band);
    patch('Faces',F, 'Vertices',V, 'VertexNormals', N, 'FaceColor','red','EdgeColor','none');
else

    if isenabled('contour',varargin{:})
        %isosurface(ls.phi,0);
        [F,V,N] = triangulate(ls.phi, ls.band);
        patch('Faces',F, 'Vertices',V, 'VertexNormals', N, 'FaceColor','red','EdgeColor','none');
    end

    hold on;

    if isenabled('narrowband',varargin{:})
        param = getparameters('narrowband', varargin{:});
        delta = 1;
        x = str2double(param);
        if ~isnan(x)
            delta = x;
        end
        [X,Y,Z] = ind2sub(size(ls), narrowband(ls));
        plot3(Y(1:delta:end),X(1:delta:end),Z(1:delta:end),'g.');
    end

    if isenabled('gradient',varargin{:})
        param = getparameters('gradient', varargin{:});
        delta = 1;
        x = str2double(param);
        if ~isnan(x)
            delta = x;
        end            
        ind = uint32(1:numel(ls.phi));
        [Dx,Dy,Dz] = ls.diff_central(ls.phi, ind);
        [rows,cols,slices] = size(ls.phi);
        Dx = reshape(Dx,rows,cols,slices)*delta;
        Dy = reshape(Dy,rows,cols,slices)*delta;
        Dz = reshape(Dz,rows,cols,slices)*delta;
        Dx = reducevolume(Dx,delta);
        Dy = reducevolume(Dy,delta);
        Dz = reducevolume(Dz,delta);
        quiver3(delta, Dx, Dy, Dz, 0, 'b');
    end
end

[rows, cols, slices] = size(ls.phi);
axis equal;
axis([0 cols 0 rows 0 slices]);
axis vis3d;
grid on;
camlight;
lighting gouraud;
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
