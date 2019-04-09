%
% Extends given interface edge velocities to the whole domain by solving the 
% Helmholtz equation.
%
% NOTE: In the document we refer to the extended velocity as 
% 
% F = \nabla u \cdot (\nabla phi / |\nabla phi|)
%
% whereas this function returns \nabla u as the extended velocity. Later in the
% advection matrix assembly (Assembler.m/assemble_advection()) we compute
% \nabla u \cot \nabla phi, which equals F|\nabla phi|.
%
% Input:
% lset      - Level set object 
% V         - Interface edge velocities
% 
% Output:
% F         - Extended velocity field (unscaled).
%

function [F] = Velocity_extension_Helmholtz( lset, V )
    
    %
    % Assemble (V,v)_\Gamma for velocity V.
    %
    
    mesh = lset.mesh_if;
    t_liquid = find(mesh.t(4,:) == 1);      % liquid phase triangles    
    t_solid = find(mesh.t(4,:) == 2);       % solid phase triangles
    if_edges = intersect( mesh.t2e(:,t_solid), mesh.t2e(:,t_liquid) );
    if_nodes = unique( mesh.edges(:,if_edges) );
    Vv = sparse( size(mesh.p,2), 1 );
    
    % Interface edge lengths.
    edges = mesh.edges(:,if_edges);
    l = mesh.p(:,edges(1,:)) - mesh.p(:,edges(2,:));
    l = sqrt( sum(l.^2) );

    % Assembling the linear form over the interface edges.
    for k = 1:size(edges,2)
        nodes = edges(:,k) ;
    
        % Velocity V constant over element. Quadrature point 0.5, hence basis
        % functions 1-x, x both equal to 0.5. Quadrature weight 1.0.
        % So qint = qw * phi(j) * abs(det(A)) * f(k), or
        qint = 0.5 * l(k) * V(k);

        Vv(nodes) = Vv(nodes) + [qint; qint];
    end
 
    %
    % Extend the velocity to the whole domain.
    %
 
    asm = Assembler();
    weak = @(u,v,ux,uy,vx,vy,x,y,h)( (ux.*vx + uy.*vy) + u.*v );
    n1 = unique( mesh.t(1:3, t_solid) );      % all solid nodes + interface
    n2 = unique( mesh.t(1:3, t_liquid) );

    % V extension in solid.  
    meshp = Mesh( mesh.p, mesh.t(:,t_solid) );
    asm.set_mesh( meshp );
    A = asm.assemble_bilin( weak );

    b1 = zeros( size(mesh.p,2), 1 );
    b1(n1) = A(n1,n1) \ Vv(n1);

    asm.set_mesh(mesh);
    grad_solid = asm.grad_nodes( b1 );
    grad_solid(:,n2) = 0;                     % zero in liquid

    % V extension in liquid.
    meshp = Mesh( mesh.p, mesh.t(:,t_liquid) );
    asm.set_mesh( meshp );
    A = asm.assemble_bilin( weak );

    b2 = zeros( size(mesh.p,2), 1 );
    b2(n2) = -A(n2,n2) \ Vv(n2);

    asm.set_mesh(mesh);
    grad_liquid = asm.grad_nodes( b2 );
    grad_liquid(:,n1) = 0;                   % zero in solid

    N = length( lset.phi );
    F = [ grad_solid(1,1:N)' + grad_liquid(1,1:N)'; 
          grad_solid(2,1:N)' + grad_liquid(2,1:N)' ]; 

end
