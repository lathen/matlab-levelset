function ls = levelset3D(varargin)
%LEVELSET3D Create a levelset3D object.
%   TODO: Comment the code
%
%   See also PROPAGATE, PLOT.

%   Author: Gunnar Läthén (gunnar.lathen@itn.liu.se)
%   $Date: 2007/10/17

% Set standard arguments
ls.phi = [];
ls.band = [];
ls.bandwidth = Inf;
ls.integrate = @euler;
ls.diff_central = @diff_central_order2;
ls.diff_upwind = @diff_upwind_order1;
ls.diff2 = @diff2_order2;
ls.reinitialize = @reinitialize_PDE;

if nargin == 0
    ls = class(ls,'levelset3D');
elseif isa(varargin{1},'levelset3D')
    ls = varargin{1};
else
    if isnumeric(varargin{1}) && ndims(varargin{1}) == 3
        ls.phi = varargin{1};
        ls.band = uint32(1:numel(varargin{1}));
    else
        error('First argument is not a numeric 3 dimensional array');
    end

    if nargin >= 2
        if isnumeric(varargin{2}) && isscalar(varargin{2})
            ls.bandwidth = varargin{2};
            ls.band = uint32(find(abs(ls.phi) <= ls.bandwidth))';
        else
            error('Second argument is not a scalar');
        end
    end
    
    if nargin >= 3
        switch lower(varargin{3})
            case 'euler'
                % Do nothing, already set as standard
            otherwise
                error('Third argument is not a valid integrator');
        end
    end

    if nargin >= 4
        switch lower(varargin{4})
            case 'firstorder'
                % Do nothing, already set as standard
            otherwise
                error('Fourth argument is not a valid spatial discretization scheme');
        end
    end
    
    if nargin >= 5
        switch lower(varargin{5})
            case 'pde'
                % Do nothing, already set as standard
            case 'fastmarching'
                ls.reinitialize = @reinitialize_fastmarching_driver;
            otherwise
                error('Fifth argument is not a valid reinitialization routine');
        end
    end

    ls = class(ls,'levelset3D');
end
