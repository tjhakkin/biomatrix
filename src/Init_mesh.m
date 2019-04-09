%
% Initializes the computational domain by constructing a triangular mesh and 
% a level set function for the initial interface.
%
% Input:
% P         - Parameters struct
%
% Output:
% lset      - Levet set object
%

function [lset] = Init_mesh(P)

    if (strcmp(P.ptype, 'stefan'))          % Frank-sphere exact test          
        lset = init_circle(P);
    elseif (strcmp(P.ptype, 'matrix'))      % matrix deposion
        p = [];
        t = [];
        phi = [];
        load(P.template);
        if (~isempty(p) && ~isempty(t))   % note: phi optional
            lset = init_template( p, t, phi, P );
        else
            fprintf('Error: Unknown template file format.\n');
            exit
        end
    else
        fprintf('Error: Unknown process type "%s".\n', P.ptype);
    end

end


%
% LOCAL FUNCTIONS
% 


%
% Creates an evenly spaces triangulated square mesh in [-2,2]^2, with a circle 
% of radius P.frank_rad centered at the origin.
%
% The mesh has (80x80)*ref nodes for ref>0. E.g. ref=2 gives 160x160 mesh.
%
function [lset] = init_circle( P )

    % Create nodes
    ref = P.Nref;
    if (ref < 1)
        ref = 1;
    end
    n = 80*ref - 1;
    [X,Y] = meshgrid( 0:1:n, 0:1:n );
    X = (2*X(:)-n) ./ n;
    Y = (2*Y(:)-n) ./ n;
    p = [2*X 2*Y]';

    % Triangulate
    t = delaunay(X,Y);

    % Node distances to the origin, circle interior nodes
    d = sqrt( p(1,:).^2 + p(2,:).^2 );
    circle_n = find( d < P.frank_rad );

    % Gather circle interior triangles
    circle_t = [];
    for i = 1:length(circle_n)
        [k,~] = find( t == circle_n(i) );
        for j = 1:length(k)
            nint = intersect( t(k(j),:), circle_n );
            if (length( nint ) == 3)
                circle_t = [ circle_t, k(j) ];
            end
        end
    end
    circle_t = unique( circle_t );

    % Assign approximate domains. Note that these are not accurate since the 
    % interface does not run along mesh edges! 
    domains = ones( 1, length(t) );
    domains( circle_t ) = 2;
    t = [t'; domains];

    % Create the mesh
    mesh = Mesh( p, t, 0 );

    % For Frank's spheres exact test, the initial interface is a circle 
    % centered at the origin.
    d = sqrt( mesh.p(1,:).^2 + mesh.p(2,:).^2 )';
    lset = Level_set( mesh, d - P.frank_rad );

    % Construct level set conforming mesh.
    lset.create_interface_mesh();

end



%
% Initializes interface when mesh nodes p and triangles t are given, and 
% optionally a precomputed level set function phi.
% 
function [lset] = init_template( p, t, phi, P )

    mesh = Mesh(p, t, P.Nref);

    if (~isempty(phi))
        mesh.t = [mesh.t(1:3,:); ones(1,length(t))];   % temporary domain defs
        mesh.t = double(mesh.t);
        mesh.edges = double(mesh.edges);
        lset = Level_set( mesh, phi );
    else
        % Constructing phi using the triangle domain definitions.
        t_liquid = find( mesh.t(4,:) == 1 );
        t_solid = find( mesh.t(4,:) == 2 );
        if_edges = intersect( mesh.t2e(:,t_liquid), mesh.t2e(:,t_solid) );
        lset = Level_set( mesh, if_edges, 2 );

        % The interface now runs along mesh edges. A smoothing is required to 
        % push the zero of phi off the mesh nodes, either by solving the heat 
        % equation or by adding a constant.
        phi = lset.phi;
        if (strcmp(P.lset_smooth, 'heat'))      % smoothen by solving heat eq.
            asm = Assembler();
            asm.set_mesh(mesh);
            weak = @(u,v,ux,uy,vx,vy,h)(u.*v);
            M = asm.assemble_bilin(weak);
            weak = @(u,v,ux,uy,vx,vy,h)(ux.*vx + uy.*vy);
            K = asm.assemble_bilin(weak);
            phi = (M + P.lset_smooth_val*K) \ (M*phi);
        else                                    % add a constant
            phi = phi + P.lset_smooth_val;          
        end

        lset.set_phi(phi);
    end
    
    % Construct level set conforming mesh.
    lset.create_interface_mesh();

end
