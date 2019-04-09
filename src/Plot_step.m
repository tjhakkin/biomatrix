%
% Plots nutrient/heat distribution and both the current and initial phase 
% interface.
%
% By default only writes .png images. Uncomment the print command at the end
% of the file to write .eps.
%
% Input:
% lset          - Level set object
% i             - Current iteration
% u             - Nutrient/heat distribution
% p0            - Initial interface vertex coordinates
% if0           - Initial interface edges as node pairs
% P             - Parameters struct
% ref_p         - Comparison object node coordinates as 2xN matrix
% save_images   - Folder to store the images, or []
%

function [] = Plot_step( lset, i, u, p0, if0, P, ref_p, save_images )

    mesh = lset.mesh_if;
    if_nodes = lset.if_nodes;

    % fig = figure('Position', [0, 0, 1200, 1200]);
    fig = figure( 'visible', 'off' );
    colormap(P.cmap);

    % Nutrient/heat distribution.
    trisurf( mesh.t(1:3,:)', mesh.p(1,:)', mesh.p(2,:)', -1*ones(1,length(u)), u );
    view(P.view);
    shading interp;

    % Fill interface interior with constant color.
    if (isfield(P, 'color_enamel') && length(P.color_enamel) == 3)
        clusters = Cluster_interface_nodes(if_nodes);
        for j = 1:size(clusters,2)
            c = clusters{j};
            patch( mesh.p(1,c(:))', mesh.p(2,c(:))', P.color_enamel, ... 
                   'EdgeColor', 'none' );
        end
    end

    % Fill initial interior with constant color.
    if (isfield(P, 'color_dentin') && length(P.color_dentin) == 3)
        clusters = Cluster_interface_nodes(if0);
        for j = 1:size(clusters,2)
            c = clusters{j};
            patch( p0(1,c(:))', p0(2,c(:))', P.color_dentin, ...
                   'EdgeColor', 'none' );
        end
    end

    % Plot initial boundary.
    if (isfield(P, 'color_EDJ') && length(P.color_EDJ) == 3)
        x0 = p0(1,:);
        y0 = p0(2,:);
        H = line( x0(if0), y0(if0), ones(size(if0)) ); 
        set(H, 'color', P.color_EDJ);
        set(H, 'linewidth', P.width_interface);      
    end

    % Add reference object, if given.
    if (length(ref_p) ~= 0)
        H = line( ref_p(1,:), ref_p(2,:), ones(size(ref_p)) );
        set(H, 'color', [0 1 0]);
        set(H, 'linewidth', P.width_interface);      
    end

    % Plot current interface.
    if (i > 0 && isfield(P, 'color_interface'))
        x = mesh.p(1,:);
        y = mesh.p(2,:);
        H = line( x(if_nodes), y(if_nodes), ones(size(if_nodes)) );        
        set(H, 'color', P.color_interface);        
        set(H, 'linewidth', P.width_interface);      
    end

    % xlim([-0.8 0.8]);
    % ylim([-0.8 0.8]);
    xlim( [min(mesh.p(1,:)), max(mesh.p(1,:))] );
    ylim( [min(mesh.p(2,:)), max(mesh.p(2,:))] );
    c = colorbar;    
    caxis(P.caxis);
    set(c, 'YTick', P.colorbar);    
    axis equal;
    drawnow;
    s = sprintf('Iteration %d', i);
    title(s);

    if (save_images)
        % s = sprintf('%s/i_%.5d_heat.eps', save_images, i);
        % print(s, '-depsc', '-r600');     % Uncomment to write .eps image.
        s = sprintf('%s/i_%.5d_%s.png', save_images, i, P.runID);
        print(s, '-dpng', '-r300'); 
        % saveas(fig, s);
    end
    
    close(fig);

end
