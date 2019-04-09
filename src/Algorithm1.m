%
% Solves diffusion-limited matrix deposion on a triangulated mesh. Alternatively,
% can be used to solve the classical Stefan problem, which is equivalent to the
% matrix deposition with zero background source and process direction reversed
% (diffusion out of the interface rather than into the interface).
%
% Simulates general matrix deposition if P.ptype = 'matrix', and Stefan problem 
% for Frank-sphere test if P.ptype = 'stefan'.
%
% NOTE: Call main_matlab() to start the simulation.
% 
% Terminology: 
% - 'boundary' refers to the outer edges of the triangulated domain, 'interface'
%   refers to the phase transition boundary between liquid and solid.
% - 'nutrient' is the diffusing element in the matrix deposition simulation, 
%   'heat' the diffusing element in Stefan problem.
%

function [] = Algorithm1( P, save_images )

fprintf('Initializing mesh, nutrient distribution...'); tic;
lset = Init_mesh(P);
u = Init_distribution(P, lset, save_images);
fprintf(' done after %.2f seconds.\n', toc);

% Store the initial interface, plot initial conditions:
if0 = lset.if_nodes;
p0 = lset.mesh_if.p;
Plot_step(lset, 0, u, p0, if0, P, [], save_images);

% For Frank-sphere test, initialize an output file for storing radiuses.
if (strcmp(P.ptype, 'stefan'))
    s = sprintf('%s/ref_circle.txt', save_images);
    fref = fopen(s, 'w');
    fprintf(fref, 'Time\tRadius\tExact_radius\n');
    fprintf(fref, '1.000000\t%f\t%f\n', P.frank_rad, P.frank_rad);
end

fprintf('Entering main loop.\n');

for i = 1:P.nIter
    fprintf('Iteration %d / %d.\n', i, P.nIter);

    fprintf('   * Computing interface velocity...'); tic;
    Ve = Stefan_condition( lset, u', P );
    fprintf(' done after %.2f seconds.\n', toc);

    fprintf('   * Extending the velocity...'); tic;
    F = Velocity_extension_Helmholtz( lset, Ve );
    fprintf(' done after %.2f seconds.\n', toc);

    fprintf('   * Updating the interface...'); tic;
    lset = Interface_growth( lset, F, P.dta );
    fprintf(' done after %.2f seconds.\n', toc);

    fprintf('   * Constructing interface-conforming mesh...'); tic;
    lset.create_interface_mesh();
    fprintf(' done after %.2f seconds.\n', toc);

    fprintf('   * Reinitializing level set...'); tic;
    lset.update_level_set_new();
    fprintf(' done after %.2f seconds.\n', toc);

    fprintf('   * Computing curvature...'); tic;
    k_full = Curvature_interpolate( lset );
    fprintf(' done after %.2f seconds.\n', toc);

    fprintf('   * Updating nutrient distribution...'); tic;
    u = Interface_condition( lset, k_full, P, u, i );
    [u,~] = Solve_heat( lset, P, u );
    fprintf(' done after %.2f seconds.\n', toc);

    if ((mod(i,P.draw_step) == 0))
        fprintf('   * Plotting, writing data...'); tic;
        ref_p = [];
        if (strcmp(P.ptype, 'stefan'))
            t_solid = find(lset.mesh_if.t(4,:) == 2);
            [ref_p, L2, area, r] = Fit_circle( lset, t_solid );
            fprintf('\n     Fitted circle radius: %f\n', r);
            time = 1 + i*P.dth;
            r_exact = P.frank_rad * sqrt(time);   % Frank-sphere exact radius
            fprintf(fref, '%f\t%f\t%f\n', time, r, r_exact);
        end
        Plot_step( lset, i, u, p0, if0, P, ref_p, save_images );
        fprintf('     done after %.2f seconds.\n', toc);
    end 
    
    fprintf('\n');
end

end     % End Stefan().
