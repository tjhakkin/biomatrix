%
% Computes 2D interface curvature. Call get_curvature().
%
classdef Curvature < handle
    
    properties ( Access = private )
    end
    
    
    
    methods ( Access = public )

    function [obj] = Curvature()
        % CURVATURE Constructor for interface curvature computation object.
    end


    
    function [normals nodes k] = get_curvature( obj, mesh, phi, if_edges )
        % GET_CURVATURE Computes curvature for 2D interface.
        %   [NORMALS, NODES, K] = GET_CURVATURE( MESH, PHI, IF_EDGES ) returns
        %   curvature K at NODES with vertex NORMALS.
        %
        %   Reference:
        %   Meyer et. al. Discrete Differential-Geometry Operators for
        %   Triangulated 2-Manifolds (2001).
        %
        %   Note: If for whatever reason the curvature cannot be computed at
        %   a node (e.g. degenerate interface), NORMALS, NODES and K will be
        %   set to zero. User should therefore check for zero nodes and normals.
    
        [normals, nodes] = obj.vertex_normals(mesh, phi, if_edges);           
        node_pairs = mesh.edges(:, if_edges);

        N = length(nodes);
        k = zeros(1,N);
        
        for i = 1:N
            node = nodes(i);
            if (node==0 || norm(normals(:,i))==0)
                continue;   % Warning has already been displayed, no need here.
            end
            
            % Find the interface edges that contain the node. Must get exactly two edges.
            [~,J] = find(node_pairs == node);
            if (length(J)<2) 
                s = sprintf('Warning: Node %d has only one edge. Skipping.', node);
                disp(s);
                continue;
            end
            if (length(J)>2)
                s = sprintf('Warning: Node %d has more than two edges. Skipping.', node);
                disp(s);
                continue;
            end
            
            % Get the neighbour nodes, indices of the neighbour normals.
            nbours = setdiff(node_pairs(:,J), node);
            
            % Neighbour vectors & distances.
            r = [ mesh.p(:,node) - mesh.p(:,nbours(1)), ...
                  mesh.p(:,node) - mesh.p(:,nbours(2)) ];
            d = [ norm(r(:,1)), norm(r(:,2)) ];
            if (d(1) == 0 || d(2) == 0)
                s = sprintf('Warning: Division by zero when computing curvature!');
                disp(s);
            end
            
            k1 = 2*(dot(normals(:,i), r(:,1))) / d(1)^2;
            k2 = 2*(dot(normals(:,i), r(:,2))) / d(2)^2;
            k(i) = (k1+k2) / 2;
        end
    end     % END get_curvature().
    

    
    function [normals, nodes] = vertex_normals(obj, mesh, phi, if_edges)
        % VERTEX_NORMALS    Computes interface vertex normals.
        %   [NORMALS, NODES] = VERTEX_NORMALS( MESH, PHI, IF_EDGES ) returns
        %   the interface vertex NORMALS at NODES.
        % 
        %   Vertex normals are computed as weighted averages of edge normals.
        
        [n_edges,~] = obj.edge_normals(mesh, phi, if_edges);
    
        edge_nodes = mesh.edges(:,if_edges);
        u_nodes = unique(edge_nodes);
        
        N = length(u_nodes);
        nodes = zeros(1,N);
        normals = zeros(2,N);
        
        % Loop through all nodes associated with the edges.
        for i = 1:N
            node = u_nodes(i);
            
            % Find the edges that contain the node. Must get exactly two edges.
            [~,J] = find(edge_nodes==node);
            if (length(J)<2) 
                s = sprintf('Warning: Node %d has only one edge. Skipping.', node);
                disp(s);
                continue;
            end
            if (length(J)>2)
                s = sprintf('Warning: Node %d has more than two edges. Skipping.', node);
                disp(s);
                continue;
            end
            
            nodes(i) = node;

            % Get the neighbour nodes, calculate edge lengths.
            nbours = setdiff(edge_nodes(:,J), node);
            if (length(nbours) < 2)
                s = sprintf( '%s(): Node %d has other than two neighbours (%d).', ...
                             mfilename(), node, length(nbours) );
                disp(s);
                edge_nodes(:,J)
                continue;
            end
            
            d1 = norm(mesh.p(:,node) - mesh.p(:,nbours(1)));
            d2 = norm(mesh.p(:,node) - mesh.p(:,nbours(2)));
            if (d1+d2 == 0) 
                s = sprintf('Warning: Edges %d, %d have zero length. Skipping.', ...
                            if_edges(J(1)), if_edges(J(2)));
                disp(s);
                continue;
            end
            
            % Interpolate the vertex normal as an average weighted by the edge lengths.
            t = d1/(d1+d2);
            
            % For surface tension using weights by lengths breaks things.
            n = (1-t)*n_edges(:,J(1)) + t*n_edges(:,J(2));    % Weight by inverse lengths.
            % n = t*n_edges(:,J(1)) + (1-t)*n_edges(:,J(2));      % Weight by lengths.
            n = n/norm(n);
            normals(:,i) = n;
        end
    end     % END vertex_normals().
        

    
    function [n_vectors, n_origins] = edge_normals( obj, mesh, phi, if_edges )
        % EDGE_NORMALS  Computes interface edge normals.
        %   [NORMALS, NODES] = EDGE_NORMALS( MESH, PHI, IF_EDGES ) returns
        %   the interface edges NORMALS with normal origin NODES corresponding
        %   the to edge mid points.
        
        N = length(if_edges);
        n_origins = zeros(2,N);
        n_vectors = zeros(2,N);
        
        % Needed for normal reorienting:
        asm = Assembler();
        asm.set_mesh(mesh);
        normals = asm.grad_nodes(phi);

        % Compute the edge normals.
        for i = 1:N
            edge = if_edges(i);
            p1 = mesh.p(:, mesh.edges(1,edge));
            p2 = mesh.p(:, mesh.edges(2,edge));

            % Normal origins.
            n_origins(:,i) = (p1-p2)/2 + p2; 
            
            % Non-oriented normals.
            vec = [0, -1; 1 0]*(p1-p2);
            n_vectors(:,i) = vec./norm(vec);
            
            % Align the normals so that they always point towards the positive level set.
            % The correct orientation is determined by looking at the value of phi at the
            % nodes belonging to the edge triangles.
            %
            % NOTE: Impossible to resolve the orientation if both of the nodes adjacent to
            % interface edge have the same sign. This can happen if the interface domain is 
            % too thin (mainly so that phi(a)=phi(b)=0).
            tris = mesh.e2t(:,edge);
            nodes = setdiff(mesh.t(1:3, tris), mesh.edges(:,edge));
            np = n_origins(:,i) + n_vectors(:,i);
            val = phi(nodes);
            if (sign(val(1))==sign(val(2)))
                s = sprintf('Warning: Ambiguous interface normal orientation at edge %d!', ...
                            edge);
                disp(s);
            end

            % Reorient the normal to point towards positive phi.
            % TODO: What we do here is to use the level set -derived normals
            % to reorient the 'exact' edge normals. Would be nice to just use
            % the level set normals.
            en = mesh.edges(:, edge);
            en = normals(:, en(1)) + normals(:, en(2));
            en = en ./ norm(en);
            if (dot(n_vectors(:,i), en) < 0)
                n_vectors(:,i) = -n_vectors(:,i);
            end

            % quiver( n_origins(1,:), n_origins(2,:), n_vectors(1,:), n_vectors(2,:) );

        end
    end     % END edge_normals().


    end % END PUBLIC METHODS.
    
end

