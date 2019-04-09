%
% Computes interface velocity V = sgn*[[\grad u \cdot n]] / L on given interface,
% where u is nutrient/heat distribution, n interface normal and L Stefan constant.
% The sign of the velocity is determined by process type: For classical Stefan 
% problem the sign is negative, for matrix deposition position.
%
% Input: 
% lset      - Level set object
% u         - Nutrient/heat distribution
% P         - Parameters struct
%
% Output:
% Ve        - Interface edges velocities
%

function [Ve] = Stefan_condition( lset, u, P )

    % Get interface triangles in both liquid and solid.
    mesh = lset.mesh_if;
    t_solid = find(mesh.t(4,:) == 2);       % solid phase triangles
    t_liquid = find(mesh.t(4,:) == 1);      % liquid phase triangles
    if_edges = intersect( mesh.t2e(:,t_solid), mesh.t2e(:,t_liquid) );
    tris = mesh.e2t(:,if_edges);
    tris = tris(:);
    i_liquid = find( mesh.t(4,tris) == 1 );
    i_solid = setdiff( 1:length(tris), i_liquid );
    tif_liquid = tris(i_liquid);
    tif_solid = tris(i_solid);

    curv = Curvature();
    [normals, ~] = curv.edge_normals( mesh, lset.phi_if, if_edges );
    asm = Assembler();
    asm.set_mesh(mesh);
    grad = asm.grad(u);

    % By default, compute the interface velocity with reverse jump condition
    % (negative sign) for the classical Stefan condition; for matrix secretion 
    % reverse to positive.
    sgn = -1;
    if (strcmp(P.ptype, 'matrix'))
        sgn = 1;
    end
    Ve = sgn * (dot(grad(:,tif_liquid), normals) - dot(grad(:,tif_solid), normals)) / P.L;

end


