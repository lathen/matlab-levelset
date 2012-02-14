function ls = levelset2D(varargin)
%LEVELSET2D  Create a levelset2D object.
%   LS = LEVELSET2D(A, [nband, temporal, spatial, reinit]) creates a 
%   levelset2D object with level set function A. This does not generally
%   require A to be a signed distance function.
%
%   If nband is specified, a narrowband is used for faster computations.
%   All points x in A such that |x| < nband are included in the narrowband.
%   Default value is Inf.
%
%   The temporal parameter can be specified to set an integration method
%   for the temporal discretization. Default is 'Euler'. Possible values
%   are: 'Euler'
%
%   The spatial parameter can be specified to set a finite difference
%   scheme for the spatial discretation. Default is 'FirstOrder'. Possible
%   values are: 'FirstOrder', 'WENO'
%
%   The reinit parameter can be specified to set a method for
%   reinitilization (resetting the level set function to a signed distance
%   function and adjusting the narrowband). Default is 'FastSweeping'.
%   Possible values are: 'FastSweeping', 'FastMarching', 'PDE'
%
%   See also PROPAGATE, PLOT.

%   Author: Gunnar Läthén (gunnar.lathen@itn.liu.se)
%   $Date: 2007/10/17

% Set some standard arguments
ls.phi = [];
ls.band = [];
ls.bandwidth = Inf;
ls.integrate = @euler;
ls.diff_central = @diff_central_order2;
ls.diff_upwind = @diff_upwind_order1;
ls.diff2 = @diff2_order2;
ls.reinitialize = @reinitialize_fastsweeping_driver;
ls.curvature = @curvature_simple;

if nargin == 0
    ls = class(ls,'levelset2D');
elseif isa(varargin{1},'levelset2D')
    ls = varargin{1};
else
    if isnumeric(varargin{1}) && ndims(varargin{1}) == 2
        ls.phi = varargin{1};
        ls.band = uint32(1:numel(varargin{1}));
    else
        error('First argument is not a numeric 2 dimensional array');
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
            case 'weno'
                ls.diff_upwind = @diff_upwind_WENO;
            otherwise
                error('Fourth argument is not a valid spatial discretization scheme');
        end
    end
    
    if nargin >= 5
        switch lower(varargin{5})
            case 'fastsweeping'
                % Do nothing, already set as standard
            case 'fastmarching'
                ls.reinitialize = @reinitialize_fastmarching_driver;
            case 'pde'
                ls.reinitialize = @reinitialize_PDE;
            otherwise
                error('Fifth argument is not a valid reinitialization routine');
        end
    end

    if nargin >= 6
        switch lower(varargin{6})
            case 'simple'
                % Do nothing, already set as standard
            case 'gaussianderivatives'
                ls.curvature = @curvature_gaussian_derivatives;
            case 'divergencenormals'
                ls.curvature = @curvature_divergence_normals_driver;
            otherwise
                error('Sixth argument is not a curvature estimation routine');
        end
    end
    
    ls = class(ls,'levelset2D');
end
