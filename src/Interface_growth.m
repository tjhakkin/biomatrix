%
% Updates the interface position by solving advection equation on the level set
% function.
%
% Input:
% lset      - Level set object
% F         - Extended velocity in the form \nabla u
% dta       - Advection time step
%
% Output:
% lset      - Updated level set object
%

function [lset] = Interface_growth( lset, F, dta )
    
    asm = Assembler();
    asm.set_mesh( lset.mesh );
    K = asm.assemble_advection( F );
    weak = @(u,v,ux,uy,vx,vy,h)( v.*u );
    M = asm.assemble_bilin( weak );
         
    % Move interface by solving advection eq. for level set b:
    % Mb. + Kb = 0
    % db/dt = -inv(M)*K*b
    % b(t+1) = b(t) - dt*inv(M)*K*b(t+1)
    % b(t+1)(1 + dt*inv(M)*K) = b(t)
    % b(t+1)(M + dt*K) = M*b(t)    
    phi = (M + dta*K) \ (M * lset.phi);
    lset.set_phi( phi );
        
end

