%
% Fits a circle of given area to the zero level set, computes the L2 error to 
% the given interface.
%
% Calls local function Domain_area().
%
% Input:
% lset      - Level set object
% t_in      - List of triangles on which to fit the circle.
%
% Output:
% p         - Circle node coordinates as 2xN matrix.
% L2        - L2 error.
% area      - Total area of triangles t_in.
% r         - Fitted circle radius
%

function [p, L2, area, r] = Fit_circle( lset, t_in )

    mesh        = lset.mesh_if;
    if_nodes    = lset.if_nodes;
    area        = Domain_area( mesh, t_in );

    % Reference circle radius.
    r = sqrt(area/pi);
    
    p1 = mesh.p(:,if_nodes(1,:));
    p2 = mesh.p(:,if_nodes(2,:));
    
    % Edge lengths.
    d = p1-p2;
    l = sqrt( d(1,:).^2 + d(2,:).^2 );
    
    % Edge mid points.
    em = (p1+p2) / 2;

    % Center the reference circle at the level set minimum.
    % [~,i] = min( lset.phi_if );
    % p0 = mesh.p(:,i);
    p0 = mean( [p1,p2], 2 );

    % Trapezoidal rule for integration over the interface.
    d = em - repmat(p0, 1, length(em));
    d = sqrt( d(1,:).^2 + d(2,:).^2 );
    
    % Distances to the circle.
    d = d - repmat(r, 1, length(d));
    
    % L2 error. Values of d are distances from quadrature points to the
    % reference circle, l contains the edge lengths.
    L2 = sqrt( d.^2 * l' );
       
    % Reference circle coordinates.
    ang = [0:0.01:2*pi];
    p = [ r*cos(ang) + p0(1); r*sin(ang) + p0(2) ];

end



%
% Returns the area of given list of triangles.
%
function [area] = Domain_area( mesh, t_in )
    p = mesh.p;
    t = mesh.t;

    % Affine map matrix B (i.e. F(x) = Bx + b) for all elements.
    B = cell(2,2);
    B{1,1} = p(1,t(2,:))' - p(1,t(1,:))';
    B{1,2} = p(1,t(3,:))' - p(1,t(1,:))';
    B{2,1} = p(2,t(2,:))' - p(2,t(1,:))';
    B{2,2} = p(2,t(3,:))' - p(2,t(1,:))';

    % Determinants.
    detB = B{1,1}.*B{2,2} - B{1,2}.*B{2,1};

    % Triangle areas.
    At = abs(detB)/2;
    
    area = sum( At(t_in) );
end
