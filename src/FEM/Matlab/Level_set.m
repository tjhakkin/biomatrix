%
% Implementation of level sets for 2D triangular meshes.
%
classdef Level_set < handle
    % Representation of level set function with associated methods.
    
    properties
        phi       % Level set values at orginal mesh nodes (Nx1)
        phi_if    % Level set at original mesh nodes + interface (zeros by def.).
        
        mesh      % Initial mesh.
        
        mesh_if   % Mesh with interface.
        if_nodes  % Interface node pairs for mesh_if.
    end

    
    
    methods (Access = public)
    
    function [obj] = Level_set( varargin )
        % LEVEL_SET Constructor for level set object.
        %   LEVEL_SET( MESH, IF_EDGES, DOMAIN ) computes the initial level set 
        %   from MESH nodes to the nodes at the edge nodes in IF_EDGES. Nodes
        %   belonging to triangles with domain definition DOMAIN are assigned
        %   negative distance. Note that as the initial zero level set runs along 
        %   the MESH edges, the level set must be advected before calling methods
        %   such as create_interface_mesh() that assume the interface does not
        %   coincide with mesh edges. 
        %
        %   LEVEL_SET( MESH, PHI ) sets the initial mesh as MESH and the level 
        %   set function as PHI. MESH does not need to contain correct domain
        %   definitions, as no level set computation is performed. 
        
        mesh = varargin{1};
        obj.mesh = mesh;

        if (nargin == 2)
            obj.phi = varargin{2};
            obj.phi_if = varargin{2};
            return
        end

        if_edges = varargin{2};
        domain = varargin{3};

        phi = zeros( size(mesh.p,2), 1 );
        phi_if = zeros( size(mesh.p,2), 1 );

        if (size(mesh.t,1) ~= 4)
            disp('ERROR: ''mesh.t'' must be of size 4xN.')
            obj.phi = phi;
            obj.phi_if = phi_if;
            return;
        end
        
        edges = mesh.edges(:, if_edges);     
        if_p = unique(edges(:)');
        if (isempty(if_p))
            obj.phi = phi;
            obj.phi_if = phi_if;
            return;
        end
        
        % Compute the distance between the interface nodes and others.
        ind = setdiff( 1:size(mesh.p,2), if_p );
        
        nx = mesh.p(1,ind);
        ny = mesh.p(2,ind);
        ix = mesh.p(1,if_p);
        iy = mesh.p(2,if_p);
        
        t1 = find( mesh.t(4,:) == domain );
        t1_p = unique( mesh.t(1:3, t1) );
        
        phi(ind) = min( sqrt(bsxfun(@minus,nx,ix').^2 + ... 
                             bsxfun(@minus,ny,iy').^2) );
        phi(t1_p) = -phi(t1_p);
        
        obj.phi = phi;
        obj.phi_if = phi;
        
    end     % END Level_set().
    
    
    
    function [] = set_phi( obj, phi )
        % SET_PHI   Sets the current phi.
        %   SET_PHI( PHI ) set PHI as the current level set function. PHI may
        %   either be with or without the interface nodes.

        if (length(obj.phi) ~= length(phi) && length(obj.phi_if) ~= length(phi))
            fprintf( 'Error: phi does not correspond to the current mesh.\n');
            return;
        end
               
        if (length(phi) == length(obj.phi_if))
            obj.phi_if = phi;
            obj.phi = phi( 1:length(obj.phi) );
        else
            N = length(obj.phi_if) - length(phi);
            obj.phi_if = [phi; zeros(N,1)];
            obj.phi = phi;
        end
    
    end     % END set_phi().
    


    function [e_if, n1, n2] = edges_crossing_interface( obj )
        % EDGES_CROSSING_INTERFACE  Returns edges and the associated nodes 
        % crossing the zero level set.
        %   [e_if, n1, n2] = EDGES_CROSSING_INTERFACE(), where e_if are the edge
        %   IDs, n1 and n2 the edge nodes.

        edges   = obj.mesh.edges;
        phi     = obj.phi;
        e1      = find( [(phi(edges(1,:)) <= 0.0) & (phi(edges(2,:)) >= 0.0)] );
        e2      = find( [(phi(edges(2,:)) <= 0.0) & (phi(edges(1,:)) >= 0.0)] );       
        e_if    = union(e1,e2);
        n1      = edges(1,e_if);
        n2      = edges(2,e_if);
    
    end     % END edges_crossing_interface()


    
    %
    % ** OBSOLETE - use update_level_set_new() instead. 
    %
    function [] = update_level_set( obj, mesh, if_nodes, p_if )
        % UPDATE_LEVEL_SET  Resets the level set to a signed distance function.
        %
        % Computes direct node-to-node distances for nodes away from the
        % interface, and projection distances to the interface in the
        % neighborhood of the interface for improved accuracy.
        %
        % Note: The pre-reset level set function phi is needed to determine the 
        % correct signs after recomputing the distances. It follows this 
        % method applies only when the level set changes are sufficiently small
        % such that there is no incorrect sign changes anywhere within the
        % domain. In other words, distances can be incorrect but signs must be 
        % correct in the pre-reset phi.
        %
        % Input:
        % mesh      - Original mesh.
        % if_nodes  - Interface nodes.
        % p_if      - Interface node positions.
        %

        N = size(mesh.p, 2);
        phi_new = zeros(N,1);
   
        if_n = unique(if_nodes);
        if (isempty(if_n))
            fprintf('Error: List of interface nodes empty.\n');
            return;
        end
        
        % Compute the distance between the interface nodes and others.
        % Splitting the set of nodes to avoid using excessive amounts of memory.
        % TODO: Adaptively split the set of nodes rather than using a fixed number 
        % of sets.

        ix = mesh.p(1,if_n);
        iy = mesh.p(2,if_n);
        if (size(p_if,2) > 0)
            ix = p_if(1,:);
            iy = p_if(2,:);        
        end
        
        nset = 256;               % Number of sets.
        bsize = ceil(N/nset);    % Size of the largest set of nodes.

        for k = 1:nset
            i = 1+(k-1)*bsize;
            j = min( k*bsize, N );
        
            ind = setdiff( i:j, if_n );
            nx = mesh.p(1,ind);
            ny = mesh.p(2,ind);

            % Request single precision to save memory.
            a = single( bsxfun(@minus,nx,ix') );
            b = single( bsxfun(@minus,ny,iy') );
            phi_new(ind) = min( sqrt(a.^2 + b.^2) );
        end
                                
        % Inner nodes have negative distance.
        idx = find( obj.phi < 0 );
        phi_new(idx) = -phi_new(idx);
        
        obj.phi_if = phi_new;
        obj.phi = phi_new( 1:length(obj.phi) );
        
    end     % END update_level_set().



  function [] = update_level_set_new( obj )
        % UPDATE_LEVEL_SET_NEW  Resets the level set to a signed distance function.
        %
        % Note: The pre-reset level set function phi is needed to determine the 
        % correct signs after recomputing the distances. It follows this 
        % method applies only when the level set changes are sufficiently small
        % such that there is no incorrect sign changes anywhere within the
        % domain. In other words, distances can be incorrect but signs must be 
        % correct in the pre-reset phi.
        %

        if_n = unique( obj.if_nodes );
        if (isempty(if_n))
            fprintf('Error: List of interface nodes empty.\n');
            return;
        end

        mesh = obj.mesh;
        mesh_if = obj.mesh_if;
        phi_new = obj.phi_if;

        ix = mesh_if.p(1,if_n);
        iy = mesh_if.p(2,if_n);
        
        % Do not update the positions of interface-crossing edge nodes yet
        % to avoid unnecessary interface moving:
        [~, n1, n2] = obj.edges_crossing_interface();
        ife_nodes = unique( [n1, n2] );     % nodes to exclude from update

        N = size( mesh.p, 2 );      % now 1:N is all non-interface nodes
        nset = 256;                 % number of computation sets
        bsize = ceil(N/nset);       % number of nodes in the largest set

        for k = 1:nset
            i = 1+(k-1)*bsize;
            j = min( k*bsize, N );
        
            ind = setdiff( i:j, ife_nodes );
            nx = mesh.p(1,ind);
            ny = mesh.p(2,ind);

            a = bsxfun( @minus, nx, ix' );
            b = bsxfun( @minus, ny, iy' );
            phi_new(ind) = min( sqrt(a.^2 + b.^2) );
        end

        % Inner nodes have negative distance.
        idx = setdiff( find(obj.phi < 0), ife_nodes );
        phi_new(idx) = -phi_new(idx);
        
        % Reset interface nodes to zero.
        phi_new(if_n) = 0;
        
        obj.phi_if = phi_new;
        obj.phi = phi_new( 1:length(obj.phi) );
        
        %
        % For each mesh node p and interface edge node pair (a,b), compute
        %
        % k = ((a-b).p + (a-b).b) / |a-b|^2
        %
        % where '.' is the dot product. Then the intersection point c with min.
        % distance to n is 
        %
        % c = k(a-b) + b
        %

        % NOTE: The projection distances are computed for nodes within distance
        % 10*mean(hK), where mean(hK) is the mean element excircle diameter.
        [~, hK] = getExCircle(mesh);
        bwidth = 10 * mean(hK);        % 10 is arbitrary, 'wide enough' usually
        upd_n = find( abs(obj.phi) < bwidth );

        p = mesh.p( :, upd_n );
        a = mesh_if.p( :, obj.if_nodes(1,:) );
        b = mesh_if.p( :, obj.if_nodes(2,:) );
        ab = a - b;

        abp = bsxfun( @times, p(1,:)', ab(1,:) ) + ...
              bsxfun( @times, p(2,:)', ab(2,:) );       % dot product (a-b).p
        abb = bsxfun( @times, b(1,:), ab(1,:) ) + ...
              bsxfun( @times, b(2,:), ab(2,:) );        % dot product (a-b).b

        k = (abp - abb) ./ sum( ab.^2 );

        % Replace node-node distances for those nodes to which there exists
        % a projection distance to an interface edge.
        % TODO: vectorize
        phi_new = obj.phi_if;

        for ii = 1:length(upd_n)

            j = intersect( find( k(ii,:) >= 0 ), find( k(ii,:) <= 1 ) );
            if (isempty(j))
                continue;
            end

            % projected nodes on edges
            c = [k(ii,j); k(ii,j)] .* ab(:,j) + b(:,j);
            d = (p(:,ii) - c).^2;
            [d_min, ~] = min( sqrt(sum(d)) );

            n = upd_n(ii);          % global node index
            if (phi_new(n) < 0)
                d_min = -d_min;
            end

            if (abs(phi_new(n)) > abs(d_min))
                phi_new(n) = d_min;
            end
        end

        obj.phi_if = phi_new;
        obj.phi = phi_new( 1:length(obj.phi) );

    end     % END update_level_set_new().

    
    
    function [rv] = create_interface_mesh( obj )
        % CREATE_INTERFACE_MESH  Constructs a level set conforming mesh.
        %   [mesh_if, phi_if, if_nodes] = CREATE_INTERFACE_MESH() expands
        %   the initial triangular mesh (defined at the constructor) to include
        %   the nodes at the zero level set.
        %
        %   Returns 0 if no errors, -1 if errors.
        %
        %   Limitations:
        %   - phi is assumed piecewise linear.
        %   - Cannot handle situations where zero level set follows an existing edge.
        %
        
        rv = 0;
        mesh = obj.mesh;
        phi = obj.phi;
    
        % First checking if the interface runs along the original mesh edges;
        % if so, then nothing to do here.
        edges = mesh.edges;
        e1 = find( [(phi(edges(1,:)) == 0.0) & (phi(edges(2,:)) == 0.0)] );
        nodes = unique( mesh.edges(:,e1) );
        if_nodes_lset = find( obj.phi_if == 0 );
        if (length(nodes) ~= 0 && length(nodes) == length(if_nodes_lset))
            obj.if_nodes = mesh.edges(:,e1);
            obj.mesh_if = mesh.copy();
            return;       
        end
    
        % Edges crossing the zero level set and the associated nodes.
        [e_if, if_n1, if_n2] = obj.edges_crossing_interface();

        % Triangles involved with the interface.
        e_is = zeros(1, size(edges,2));
        e_is(e_if) = 1;
        active_t = find( sum(e_is(mesh.t2e)) > 0 );

        % Find the position of zero level set at the edges.
        u1 = phi( if_n1 );
        n1(1,:) = mesh.p(1, if_n1);
        n1(2,:) = mesh.p(2, if_n1);

        u2 = phi( if_n2 );
        n2(1,:) = mesh.p(1, if_n2);
        n2(2,:) = mesh.p(2, if_n2);

        % Compute zero level set vertices.
        t = (0 - u1)./(u2-u1);
        for i = 1:2
            p_is(i,:) = n2(i,:).*t' + n1(i,:).*(1-t');
        end
        p_new = [mesh.p, p_is];

        t_nans = find( isnan(t(i)) );
        t_zeros = find(t==0);
        if (length(t_nans) || length(t_zeros))
            fprintf('Warning: Zero level set coincides with the orginal mesh! Behaviour undefined.\n');
        end
    
        % Exclude the triangles that will be split from the new triangles list.
        idx = setdiff(1:size(mesh.t,2), active_t);
        t_new = mesh.t(:, idx);

        N = size(mesh.p,2);
        if_nodes = [];    

        % Construct the new triangles.
        for i = 1:size(active_t,2)
        
            % Given ith triangle, get the indices of the edges involved with
            % the interface.
            t_edges = mesh.t2e(:,active_t(i));
            idx = [ find( e_if == t_edges(1) ), ...
                    find( e_if == t_edges(2) ), ... 
                    find( e_if == t_edges(3) ) ];
            
            % Retriangulating the ith triangle by creating an edge between
            % the if nodes, and another edge between one of the if nodes
            % and the triangle node of the non-involved edge.
        
            edge_free = setdiff(t_edges, e_if(idx));    % Non-involved edge.
            node_free = setdiff( mesh.edges(:, edge_free), ... 
                                 mesh.edges(:, e_if(idx(1))) );
        
            % Each existing triangle becomes three triangles. If the level
            % set value of the non-interface node is positive the new
            % triangle is marked as domain 1, else domain 2.

            % Triangle 1.
            node_common = intersect( mesh.edges(:, edge_free), ...
                                     mesh.edges(:, e_if(idx(1))) );
            domain = 2;
            if (phi(node_common) > 0)
                domain = 1;
            end
            t_new = [t_new [idx(1)+N, node_free, node_common, domain]'];
        
            % Triangle 2.
            if_nodes = [if_nodes, [idx(1)+N, idx(2)+N]'];
            domain = 2;
            if (phi(node_free) > 0)
                domain = 1;
            end
            t_new = [t_new [idx(1)+N, node_free, idx(2)+N, domain]'];
        
            % Triangle 3.
            node_common = intersect( mesh.edges(:, e_if(idx(1))), ...
                                     mesh.edges(:, e_if(idx(2))) );
            domain = 2;
            if (phi(node_common) > 0)
                domain = 1;
            end
            t_new = [t_new [idx(1)+N, idx(2)+N, node_common, domain]'];            
           
        end
     
        mesh_if = Mesh( p_new, t_new );

        % The last entries in phi are the new zero level set.
        N = size(mesh_if.p,2) - size(phi,1);

        % Level set zero at the interface by def.
        phi_if = [phi; zeros(N,1)];
    
        % Update the domains for all triangles in the new mesh.
        phis = sum( phi_if(mesh_if.t(1:3,:)) );
        phis = find( phis > 0 );
        mesh_if.t(4,:) = 2;
        mesh_if.t(4,phis) = 1;

        obj.phi = phi;
        obj.phi_if = phi_if;
        obj.mesh_if = mesh_if;
        obj.if_nodes = if_nodes;
            
    end     % END Level_set_2D().
    
    
    
    function [] = plot( obj )
        % PLOT  Plots the current level set function. 
        %
        
        trisurf( obj.mesh.t(1:3,:)', obj.mesh.p(1,:), obj.mesh.p(2,:), obj.phi );

    end     % END plot()
    
    end     % END PUBLIC METHODS.



    methods (Access = private)        
    end     % END PRIVATE METHODS.
    
end

