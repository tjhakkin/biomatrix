% Mika Juntunen / 14.2.2012
%
% This function computes the diameter of the circle around each triangle in
% the mesh.
%
% If the sides of the triangle are a,b, and c, then the radius of this
% circle is a*b*c/sqrt( 16*s*(s-a)*(s-b)*(s-c) ) in which s = (a+b+c)/2 
%
% Input:
% mesh - mesh.p has nodes, mesh.t has triangles
%
% Output:
% HMax - maximum of the diameters
% HK   - diameter of each triangle
%
% Note: If you want only HMax:
%       HMax = getExCircle( mesh ) ;
%       If you want only diameter in each triangle:
%       [~,HK] = getExCircle( mesh ) ;
%

function [HMax, HK] = getExCircle( mesh )

% Get sides
a = sqrt( sum( (mesh.p(:, mesh.t(1,:)) - mesh.p(:, mesh.t(2,:))).^2 ) ) ; 
b = sqrt( sum( (mesh.p(:, mesh.t(2,:)) - mesh.p(:, mesh.t(3,:))).^2 ) ) ; 
c = sqrt( sum( (mesh.p(:, mesh.t(3,:)) - mesh.p(:, mesh.t(1,:))).^2 ) ) ; 

% Compute s
s = (a+b+c)/2 ;

% Compute diameters
HK = 2*a.*b.*c./sqrt( 16*s.*(s-a).*(s-b).*(s-c) ) ;

% Get max
HMax = max( HK ) ;

HK=HK(:);
end
