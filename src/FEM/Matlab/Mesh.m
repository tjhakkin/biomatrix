classdef Mesh < handle
    % A representation of triangular mesh with associated methods.

    properties ( Access = public )
        p           % Vertex coordinates, 2xN.
        t           % Triangles, 3xN or 4xN with domain definitions.
        edges       % Edges, 2xN.
        t2e         % Edges of triangles, 3xN.
        e2t         % Triangles of edges, 2xN.
    end
    
    
    methods ( Access = public )

    function [obj] = Mesh(p, t, nref)
        % MESH  Constructor for a mesh. 
        %   MESH = Mesh(P,T) initializes and returns a mesh object given vertex 
        %   positions P and triangles T.
        %
        %   MESH = Mesh(P,T,Nref) initializes and returns a mesh object that
        %   has been refined Nref times.
        %
        %   Adapted from inittri.m by Hannukainen & Juntunen (2008). 
        
        % Retain domain indicators if present.
        mat = [];
        if (size(t,1) == 4)
            mat = t(4,:);
        end
        
        t = sort( t(1:3,:) );
        t = [t; mat];

        % Initalize size variables
        Nt = size(t,2);
        e = [1 2; 2 3; 1 3]';

        % Construct edges, triangles to edges (t2e).
        edges = [];
        for i = 1:size(e,2)
            sort_edges = sort( [t(e(1,i),:); t(e(2,i),:)], 1 );
            edges = [ edges sort_edges ];
        end
        
        [edges,~,t2e] = unique( sort(edges,1)', 'rows' );            
        edges = edges';
        t2e = reshape( t2e(:)', Nt, 3 )';

        % Edges to triangles (e2t).
        e = [t2e(1,:) t2e(2,:)  t2e(3,:)];
        tt = repmat([1:Nt],1,3);

        [ef, If] = unique(e,'first');
        [el, Il] = unique(e,'last');

        e2t(1,ef) = tt(If);
        e2t(2,el) = tt(Il);

        find_edges = find( (e2t(1,:) - e2t(2,:)) == 0 );
        e2t(2, find_edges) = 0;
        
        % Add to object properties.
        obj.p = p;
        obj.t = t;
        obj.edges = edges;
        obj.e2t = e2t;
        obj.t2e = t2e;
        
        % Refines the mesh.
        if (nargin == 3)
            for i = 1:nref
                obj.refine();
            end
        end
    end


    
    function [mesh] = copy(obj)
        % COPY  Returns a deeps copy of the current mesh.
        
        mesh = Mesh(obj.p, obj.t, 0);    
    end

    

    function [] = refine(obj)
        % REFINE    Uniformly refines the mesh.
        %
        %   Adapted from refine_tri_IP.m by Antti Hannukainen (2008).
        
        t = obj.t;
        p = obj.p;
        e = obj.edges;
        t2e = obj.t2e;
        
        % Split edges.
        for n = 1:2
            e_nodes(n,:) = sum( [p(n,e(1,:)) ; p(n,e(2,:))] ) / 2;
        end          
        rmesh.p = [p  e_nodes];
        
        % Create new elements. Edges as n1->n2, n2->n3, n1->n3.  
        eb = size(p,2);
        
        if (size(t,1) == 4)
            new_t = [t(1,:) ; t2e(1,:)+eb ; t2e(3,:)+eb ; t(4,:)];
            rmesh.t = [new_t];
            
            new_t = [t(2,:) ; t2e(1,:)+eb ; t2e(2,:)+eb ; t(4,:)];
            rmesh.t = [rmesh.t  new_t];
            
            new_t = [t(3,:) ; t2e(3,:)+eb ; t2e(2,:)+eb ; t(4,:)];
            rmesh.t = [rmesh.t  new_t];
            
            new_t = [t2e(1,:)+eb ; t2e(2,:)+eb ; t2e(3,:)+eb ;  t(4,:)];
            rmesh.t = [rmesh.t  new_t];
        else
            new_t = [t(1,:) ; t2e(1,:)+eb ; t2e(3,:)+eb];
            rmesh.t = [new_t];
            
            new_t = [t(2,:) ; t2e(1,:)+eb ; t2e(2,:)+eb];
            rmesh.t = [rmesh.t  new_t];
            
            new_t = [t(3,:) ; t2e(3,:)+eb ; t2e(2,:)+eb];
            rmesh.t = [rmesh.t  new_t];
            
            new_t = [t2e(1,:)+eb ; t2e(2,:)+eb ; t2e(3,:)+eb];
            rmesh.t = [rmesh.t  new_t];
        end
        
        % Populate the refined mesh, update current mesh properties.
        rmesh = Mesh( rmesh.p, rmesh.t );
        obj.p = rmesh.p;
        obj.t = rmesh.t;
        obj.edges = rmesh.edges;
        obj.t2e = rmesh.t2e;
        obj.e2t = rmesh.e2t;
    end
   
    end     % END PUBLIC METHODS.
    
    
    
    methods ( Access = private )

    end % END PRIVATE METHODS.
        
end
