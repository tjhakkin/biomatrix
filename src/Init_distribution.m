%
% Initialized nutrient/heat distribution.
%
% Input:
% P             - Parameters struct
% lset          - Level set object
% save_images   - Folder to store the images, or []
%
% Output:
% u             - Initial nutrient/heat distribution
%

function [u] = Init_distribution( P, lset, save_images )

    mesh = lset.mesh_if;

    if (strcmp(P.ptype, 'stefan'))          % Frank-sphere exact test
        u = Frank_init( lset, P );
        k_full = Curvature_interpolate( lset );
        u = Interface_condition(lset, k_full, P, u, 1);
    elseif (strcmp(P.ptype, 'matrix'))      % matrix deposion
        u = P.u1 * ones(size(mesh.p,2), 1);
        if (strcmp(P.boundary, 'Mixed'))
            % Obtain approximately steady state solution for the initial nutrient 
            % distribution.
            P0 = struct(P);
            P0.init_heat = 1;
            P0.dth = 1e-2;       % TODO: 2/dx
            ndelta = 1;
            fprintf('Smoothing initial distribution.\n');
            while (ndelta > 1e-4)
                [u, ndelta] = Solve_heat( lset, P0, u );
                fprintf('   norm(u-un): %f\n', ndelta);
            end
        end

        t_solid = find(mesh.t(4,:) == 2);   % solid phase triangles
        n_solid = unique(mesh.t(1:3, t_solid));
        u(n_solid) = P.u2;
        k_full = Curvature_interpolate( lset );
        u = Interface_condition(lset, k_full, P, u, 1);
        % Smooth the initial temperature just created by the amount equal to one
        % iteration of the heat solver in the main loop:
        [u,~] = Solve_heat( lset, P, u );
    else
        fprintf('Error: Unknown process type "%s".\n', P.ptype);
    end
    
end
