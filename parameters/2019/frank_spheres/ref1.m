%
% Frank sphere exact test on 80x80 mesh.
%

% Neumann boundary means no-flux ('wall'), Dirichlet is a sink/source.
P.boundary = 'Neumann';         % Outer boundary condition type
P.bc = -1;                      % Boundary condition (if P.boundary = 'Dirichlet')

% If the initial interface is too jagged, some smoothing could be warranted.
P.lset_smooth = 'heat';         % Smooth initial level set by solving heat eq.
P.lset_smooth_val = 0 * 1e-3;   % Smoothing time step (0 means no smoothing)

P.dth = 1e-1;                   % Heat time step
P.dta = P.dth;                  % Advection time step

P.dh = 1.0;                     % Heat diffusion rate
P.L = 1.0;                      % Stefan constant
P.eps_C = 0.0;                  % Surface tension coefficient
P.source = 0;                   % Constant heat source in water

P.nIter = 20;                   % Number of main loop iterations

%
% Advanced / unused / unimplemented.
%

P.ptype = 'stefan';         % Process type
P.hexa_nodes = 100;         % Resolution of hexa template, if applicable.
P.Nref = 1;                 % Number of mesh refinement rounds
P.u1 = -1;                  % Initial heat (domain 1) 
P.u2 = 0;                   % Initial heat (domain 2)
P.ifc = 0;                  % Interface condition
P.eps_V = 0.000;            % Molecular kinetic coefficient (NOT IMPLEMENTED)
P.frank_rad = 0.25;         % Initial circle radius for Frank's spheres

%
% Visualization
%

P.draw_step = 2;                  % Plot every nth iteration

P.color_interface = [1.0 0.0 0.0];
P.color_EDJ       = [0.23922 0.74902 0.98824];
P.color_dentin    = [0.5 0.5 0.5];
P.color_enamel    = [0.86275 0.86275 0.86275];

P.width_interface = 1;

ramp = [0:1/63:1]';
P.cmap = ramp * [1 1 1];        % Energy/nutrient visualization

P.caxis = [-1 0];               % Colorbar caxis
P.colorbar = [-1.0 -0.5 0.0];   % Colorbar YTick
P.view = [180, 90];             % View orientation

