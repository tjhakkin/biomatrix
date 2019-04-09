%
% Sets values of nutrition/heat u at the interface nodes, needed by both 
% Stefan_condition() and Solve_heat().
%
% Input:
% lset      - Level set object
% k         - Curvature
% P         - Parameters struct
% u         - Nutrient/heat distribution
% i         - Current iteration
%
% Output:
% un        - Updated nutrient/heat distribution.
%

function [un] = Interface_condition( lset, k, P, u, i )

    % Initialize un on static mesh + interface nodes
    N = size(lset.mesh_if.p, 2);
    un = zeros(N, 1);

    % Copy non-interface values to un
    N = size(lset.mesh.p, 2);
    un(1:N) = u(1:N);

    % By default, set the sign of the interfacial tension for the classical
    % Stefan problem; for matrix secretion reverse to positive.
    sgn = -1;
    if (strcmp(P.ptype, 'matrix'))
        sgn = 1;
    end

    % Assign interface node values
    F = unique(lset.if_nodes)';                 % interface nodes
    if (length(P.eps_C) > 1)
        un(F) = P.ifc + sgn*P.eps_C(i) * k(F);  % surface tension array
    elseif (P.eps_C ~= 0)
        un(F) = P.ifc + sgn*P.eps_C * k(F);     % fixed surface tension
    else
        un(F) = P.ifc;                          % no tension (fixed temp.)
    end

end
