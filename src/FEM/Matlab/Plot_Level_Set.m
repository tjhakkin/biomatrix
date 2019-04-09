%
% Input:
% mesh      - Mesh.
% phi       - Level set values, used for both z and colors by default.
% if_nodes  - Interface node pairs.
% approx_IF - Approximated IF nodes (optinal)
% values    - Alternative values used for colors (optional)
%

function [] = Plot_Level_Set( mesh, phi, if_nodes, if_nodes_new, values )
    figsize = [1300, 1300];
    delta = [500,50,0,0];
    f = figure('Position', [0, 0, figsize]+delta)
    hold;

    if (length(values) > 0) 
        trisurf(mesh.t(1:3,:)', mesh.p(1,:), mesh.p(2,:), phi, values);
    else
        trisurf(mesh.t(1:3,:)', mesh.p(1,:), mesh.p(2,:), phi);
    end
    
    axis equal;
    view(0,90);
    x = mesh.p(1,:);
    y = mesh.p(2,:);
    H = line( x(if_nodes), y(if_nodes) ); 
    set(H,'color',[1 0 0]);
    set(H,'linewidth',2);

    if (length(if_nodes_new) > 0)
        p = unique(if_nodes_new);
        x = mesh.p(1,:);
        y = mesh.p(2,:);
        plot3( x(p), y(p), ones(length(p)), 'w.', 'MarkerSize', 15)

        H = line( x(if_nodes_new), y(if_nodes_new), ones(size(if_nodes_new)) ); 
        set(H,'color',[1 1 1]);
        set(H,'linewidth',2);
    end
end

