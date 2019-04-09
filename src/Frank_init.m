%
% Fills the initial heat distribution at time t=1 for Frank-sphere exact test.
%

function [u] = Frank_init( lset, P )

    f = @(t)exp(-t) ./ t;
    s0 = P.frank_rad;
    Fs0 = integral( f, s0^2/4, inf );
    Tinf = -0.05709187113307;
    T = @(r)( Tinf * (1 - (integral( f, r^2/4, inf ) / Fs0)) );

    p = lset.mesh.p;
    circle_n = find( lset.phi < 0 );

    u = zeros( size(p,2), 1 );
    ext = setdiff( 1:size(p,2), circle_n );
    d = sqrt( p(1,ext).^2 + p(2,ext).^2 );
    for i = 1:length(ext)
        u(ext(i)) = T(d(i));
    end

end
