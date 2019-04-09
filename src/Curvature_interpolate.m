%
% Computes curvate from the level set function at the static mesh nodes, 
% interpolates the curvature at the interface nodes.
%
% Input:
% lset      - Level set object
%
% Output:
% k_full    - Curvature at mesh nodes and interpolated at interface nodes
%

function [k_full] = Curvature_interpolate( lset )
    
    %
    % Compute curvature at all mesh nodes.
    %

    asm = Assembler();
    asm.set_mesh( lset.mesh );
    grad = asm.grad_nodes( lset.phi );

    % Curvature as k = \nabla \cdot (\nabla \phi / |\nabla \phi|)
    gn = grad ./ sqrt( grad(1,:).^2 + grad(2,:).^2 );
    gn_x = asm.grad_nodes( gn(1,:) );
    gn_y = asm.grad_nodes( gn(2,:) );
    k = gn_x(1,:) + gn_y(2,:);

    % 
    % Interpolate curvature at interface nodes.
    %

    [~, n1, n2] = lset.edges_crossing_interface();
    u1 = lset.phi( n1 );
    u2 = lset.phi( n2 );
    t = u1 ./ (u1 - u2);
    k1 = k( n1 );
    k2 = k( n2 );
    k_interp = k2.*t' + k1.*(1-t');

    k_full = [k k_interp];

end
