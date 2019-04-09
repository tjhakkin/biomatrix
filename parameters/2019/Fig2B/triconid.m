P.template = '../data/2019/Fig2B/tri200fin.mat';

% Neumann boundary means no-flux ('wall'), Dirichlet is a sink/source.
P.boundary = 'Neumann';         % Outer boundary condition type
P.bc = -1;                      % Boundary condition (if P.boundary = 'Dirichlet')

% If the initial interface is too jagged, some smoothing could be warranted.
P.lset_smooth = 'heat';         % Smooth initial level set by solving heat eq.
P.lset_smooth_val = 0 * 1e-4;   % Smoothing time step (0 means no smoothing)

P.dth = 1e-3;                   % Nutrient diffusion time step
P.dta = 1e-3;                   % Interface advection time step

P.dh = 1.0;                     % Nutrient diffusion rate
P.L = 1.5;                      % Stefan constant
P.eps_C = 0.002;                % Surface tension coefficient
P.source = 40;                  % Constant nutrient source in liquid

P.nIter = 9;                    % Number of main loop iterations

%
% Advanced / unused / unimplemented.
%

P.ptype = 'matrix';         % Process type
P.hexa_nodes = 100;         % Resolution of hexa template, if applicable.
P.Nref = 0;                 % Number of mesh refinement rounds
P.u1 = 1;                   % Initial nutrient (domain 1) 
P.u2 = 0;                   % Initial nutrient (domain 2)
P.ifc = 0;                  % Interface condition
P.eps_V = 0.000;            % Molecular kinetic coefficient (NOT IMPLEMENTED)

%
% Visualization
%

P.draw_step = 1;                % Plot every iteration

P.color_interface = [0.0 0.0 0.0];
P.color_EDJ       = [0.23922 0.74902 0.98824];
P.color_dentin    = [0.5 0.5 0.5];
P.color_enamel    = [0.98 0.98 0.98];

P.width_interface = 1;

P.cmap = [0:1/63:1]' * [0.75294 0.75294 0.75294];  % Nutrient visualization

P.caxis = [0 1];                % Colorbar caxis
P.colorbar = [0.0 0.5 1.0];     % Colorbar YTick
P.view = [180, 90];             % View orientation

