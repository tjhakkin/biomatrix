%
% Solves one step of heat equation with a source.
%
% NOTE: Interface_condition() must be called before this function to set the 
% node values at the interface.
%
% Input:
% lset      - Level set object
% P         - Parameters struct
% u         - Nutrient/heat distribution at time t
%
% Output:
% un        - Updated nutrient/heat distribution
% ndelta    - norm(u - un)
%

function [un, ndelta] = Solve_heat( lset, P, u )
    
    if (isfield(P, 'init_heat'))
        mesh = lset.mesh;
        F = [];
    else
        mesh = lset.mesh_if;
        F = unique(lset.if_nodes)';
    end

    N = size(mesh.p, 2);
    I = setdiff(1:N, F);                % nodes in which to solve

    asm = Assembler();
    asm.set_mesh(mesh);
    weak = @(u,v,ux,uy,vx,vy,h)(u.*v);
    M = asm.assemble_bilin(weak);
    weak = @(u,v,ux,uy,vx,vy,h)(ux.*vx + uy.*vy);
    K = asm.assemble_bilin(weak);
    
    s = zeros(N,1);
    ext = find(lset.phi_if > 0);
    s(ext) = P.source;                  % background production
    S = asm.assemble_lin(s, @(v,vx,vy,h)(v));

    un = u;

    %
    % Three boundary condition modes supported: 
    % Pure Dirichlet, pure Neumann, mixed Dirichlet/Neumann (default).
    %
    if (strcmp(P.boundary, 'Dirichlet')) 
        b = find(mesh.e2t(2,:) == 0);
        B = mesh.edges(:,b);
        B = unique(B(:))';              % boundary nodes
        un(B) = P.bc;                   % Dirichlet boundary condition
        B = union(B, F);
    elseif (strcmp(P.boundary, 'Neumann'))
        I = setdiff(1:N, F);            % non-interface nodes
        B = F;
    else
        % Mixed Dirichlet, Neumann boundaries. 
        % Currently only supports having Dirichlet at the yMax boundary, Neumann 
        % elsewhere.
        Bmax = find(lset.mesh.p(2,:) >= max(max(lset.mesh.p(2,:)))-eps);
        un(Bmax) = P.bc;
        Bmin = find(lset.mesh.p(2,:) <= min(min(lset.mesh.p(2,:)))+eps);
        un(Bmin) = P.u1;

        B = union([Bmax, Bmin], F);     % Dirichlet boundary + interface nodes
        I = setdiff(1:N, B);
    end

    % Implicit Euler / IMEX
    a = P.dh;
    dt = P.dth;
    un(I) = (M(I,I) + a*dt*K(I,I)) \ (M(I,I)*u(I) + dt*(S(I) - a*K(I,B)*u(B)));

    ndelta = norm(u - un);

end

