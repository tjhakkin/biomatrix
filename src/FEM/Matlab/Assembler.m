%
% 2D FEM matrix assembler for P1 triangular elements.
%
% ** Basic usage **
%
% - Construct the assembler:
% asm = Assembler().
%
% - Set the current mesh:
% asm.set_mesh( mesh );
%
% - Assemble FEM matrix for e.g. (\nabla u, \nabla v):
% weak = @(u,v,ux,uy,vx,vy,h)(ux.*vx + uy.*vy);
% K = asm.assemble_bilin( weak );
%
% - Compute gradient of V at all elements:
% grad = asm.grad( V );
%

classdef Assembler < handle

    properties ( Access = private )
        qp          % Quadrature points.
        qw          % Quadrature weights.
        phi         % Basis function values at quadrature points.
        gphi        % Basis function gradients at quadrature points.
        gphi_global % Global basis function gradients.
        
        mesh        % Mesh object containing: p, t, e2t, t2e, edges.
        B           % Reference -> global element affine map matrices.
        detB        % Determinants of B.
        invBt       % Inverse transposes of B.
    end


        
    methods ( Access = public )
    
    function [obj] = Assembler() 
        % ASSEMBLER Constructor for FEM matrix assembler.
        %   asm = Assembler() initializes and returns the FEM matrix assembler.
        
        % Quadrature points.
        obj.qp = [1/2, 1/2, 0;
                  0, 1/2, 1/2];          
          
        % Quadrature weights.
        obj.qw = [1/6, 1/6, 1/6];
        
        % Basis function values at quadrature points.
        obj.phi{1} = [1,1,1] - obj.qp(1,:) - obj.qp(2,:);
        obj.phi{2} = obj.qp(1,:);
        obj.phi{3} = obj.qp(2,:);
        
        % Basis function gradients at quadrature points.
        obj.gphi{1} = [-1, -1, -1;
                       -1, -1, -1];
        obj.gphi{2} = [1, 1, 1;
                       0, 0, 0];
        obj.gphi{3} = [0, 0, 0;
                       1, 1, 1];
    
    end     % END Assembler()
    
    
    
    function [] = set_mesh( obj, mesh )
        % SET_MESH  Sets the current mesh used by assembly functions.
        
        p = mesh.p;
        t = mesh.t(1:3,:);

        % Affine map matrix B (i.e. F(x) = Bx + b) for all elements.
        B = cell(2,2);
        B{1,1} = p(1,t(2,:))' - p(1,t(1,:))';
        B{1,2} = p(1,t(3,:))' - p(1,t(1,:))';
        B{2,1} = p(2,t(2,:))' - p(2,t(1,:))';
        B{2,2} = p(2,t(3,:))' - p(2,t(1,:))';

        % Determinants.
        detB = B{1,1}.*B{2,2} - B{1,2}.*B{2,1};
            
        % Inverse transposes.
        invBt = cell(2,2);
        invBt{1,1} = B{2,2} ./ detB;
        invBt{2,1} = -B{1,2} ./ detB;
        invBt{1,2} = -B{2,1} ./ detB;
        invBt{2,2} = B{1,1} ./ detB;
    
        % Global basis function gradients.
        gphi_global = cell(3,2);
        for j = 1:3
           gphi_global{j,1} = invBt{1,1}*obj.gphi{j}(1,:) + ...
                              invBt{1,2}*obj.gphi{j}(2,:);
           gphi_global{j,2} = invBt{2,1}*obj.gphi{j}(1,:) + ...
                              invBt{2,2}*obj.gphi{j}(2,:);
        end
        
        obj.B = B;
        obj.detB = detB;
        obj.invBt = invBt;
        obj.gphi_global = gphi_global;  
        obj.mesh = mesh;    
              
    end     % END set_mesh()
    
    
    
    function [l] = assemble_lin( obj, d, f )
        % ASSEMBLE_LIN assembles linear form vector.
        %   L = ASSEMBLE_LIN(D,F) returns the linear form vector assembled for
        %   for function handle F with parameters (v,vx,vy,h) and node values D.
    
        % Mesh parameter.
        h = repmat( sqrt(abs(obj.detB)), 1, size(obj.qp,2) );
              
        mesh = obj.mesh;
        nv = size(mesh.p,2);
        nt = size(mesh.t,2);
        l = zeros(nv,1);
        
        % d contains contains node values. Interpolating at quadrature points.
        V12 = d( mesh.t([1,2],:) )';
        V23 = d( mesh.t([2,3],:) )';
        V31 = d( mesh.t([3,1],:) )';
        V = [ mean(V12,2) mean(V23,2) mean(V31,2) ];       
    
        for i = 1:3
            % Test function at quadrature points.
            v = repmat(obj.phi{i}, nt, 1);
            
            % B_K^{-T} \nabla \phi_i for the test function.
            vx = obj.gphi_global{i,1};
            vy = obj.gphi_global{i,2};
            
            % Compute quadrature integrals over all elements.
            qint = V .* f(v,vx,vy,h) * obj.qw' .* abs(obj.detB);
            
            l = l + sparse( obj.mesh.t(i,:), ones(1,nt), qint, nv, 1 );
            
        end
    
    end     % END assemble_lin().
    
    
    
    function [A] = assemble_bilin( obj, f )
        % ASSEMBLE_BILIN Assembles bilinear form matrix.
        %   A = ASSEMBLE_BILIN(F) returns the bilinear form matrix assembled
        %   for function handle F with parameters (u,v,ux,uy,vx,vy,h).
        %   
        %   NOTE: Set the mesh by calling set_mesh() before calling this!

        % For code readability it would be nice to e.g. explicitly loop over
        % the elements, but it destroys performance in Matlab. Instead the 
        % computations will be performed over all elements using vectorization.
    
        % Mesh parameter.
        h = repmat( sqrt(abs(obj.detB)), 1, size(obj.qp,2) );
              
        np = size(obj.mesh.p,2);
        nt = size(obj.mesh.t,2);
        A = sparse(np,np);
    
        for j = 1:3
            % B_K^{-T} \nabla \phi_j for function u over all elements at 
            % the quadrature points.
            ux = obj.gphi_global{j,1};  % x components of B_K^{-T} \nabla \phi_j
            uy = obj.gphi_global{j,2};  % y components of B_K^{-T} \nabla \phi_j

            % Function values at quadrature points.
            u = repmat(obj.phi{j}, nt, 1);
            
            for i = 1:3
                % B_K^{-T} \nabla \phi_i for the test function.
                vx = obj.gphi_global{i,1};
                vy = obj.gphi_global{i,2};
                
                % Test function at quadrature points.
                v = repmat(obj.phi{i}, nt, 1);
                
                % Compute quadrature integrals over all elements.
                qint = f(u,v,ux,uy,vx,vy,h) * obj.qw' .* abs(obj.detB);
                
                % Assign the integration values to the nodes associated with
                % the elements.
                A = A + sparse(obj.mesh.t(i,:), obj.mesh.t(j,:), qint, np, np);
            end
        end
    
    end     % END assemble_heat().

 
 
    function [A] = assemble_advection( obj, V )
        % ASSEMBLE_ADVECTION Assembles advection term V(ux.*v, uy.*v), 
        % where V is the velocity.
        %
        % TODO: Combine this function with the general bilinear assembler
        % assemble_bilin().
           
        mesh = obj.mesh;
        nv = size(mesh.p,2);
        nt = size(mesh.t,2);
        A = sparse(nv,nv);
        
        % Velocity either as 2nv x 1 array for per node velocities, or
        % as 2nt x 1 for per element velocities.
        if (size(V,1) == 2*nv)
            % Computing the contributions of the node velocities at quadrature 
            % points:
            V_x = [ V( mesh.t(1,:) ) * obj.phi{1} + ... 
                    V( mesh.t(2,:) ) * obj.phi{2} + ... 
                    V( mesh.t(3,:) ) * obj.phi{3} ];
                    
            V_y = [ V( mesh.t(1,:) + nv ) * obj.phi{1} + ... 
                    V( mesh.t(2,:) + nv ) * obj.phi{2} + ... 
                    V( mesh.t(3,:) + nv ) * obj.phi{3} ];
  
        elseif (size(V,1) == 2*nt)
            % Setting velocities at the three quadrature points separately.
            V = [V V V];
            V_x = V( 1:nt, : );
            V_y = V( nt+1:2*nt, : );
            
        else
            s = sprintf( ['Error: Incorrent number of entries in V (%d:%d);', ...
                          ' must be either (%d:1) or (%d:1)'], ...
                         size(V,1), size(V,2), 2*nv, 2*nt );
            disp(s);
            return;
        end
        
        for j = 1:3
            % B_K^{-T} \nabla \phi_j for function u over all elements at 
            % the quadrature points.
            ux = obj.gphi_global{j,1};  % x components of B_K^{-T} \nabla \phi_j
            uy = obj.gphi_global{j,2};  % y components of B_K^{-T} \nabla \phi_j
            
            for i = 1:3
                % Test function at quadrature points.
                v = repmat(obj.phi{i}, nt, 1);
                
                % Compute quadrature integrals over all elements.
                qint = (V_x.*ux.*v + V_y.*uy.*v) * obj.qw' .* abs(obj.detB);
                
                % Assign the integration values to the nodes associated with
                % the elements.
                A = A + sparse(mesh.t(i,:), mesh.t(j,:), qint, nv, nv);
            end
        end
    end     % END assemble_advection().
 
    
 
    function [g] = grad( obj, u )
        % GRAD Returns the gradient of u over all triangles.
        %
        % Note: Applies to P1 elements only.

        if (issparse(u))
            u = full(u);
        end
        if (size(u,1) > size(u,2))
            u = u';
        end

        % Assuming here the gradients are constants over the elements.
        B1 = [ obj.gphi_global{1,1}(:,1) obj.gphi_global{1,2}(:,1) ]';
        B2 = [ obj.gphi_global{2,1}(:,1) obj.gphi_global{2,2}(:,1) ]';
        B3 = [ obj.gphi_global{3,1}(:,1) obj.gphi_global{3,2}(:,1) ]';

        t = obj.mesh.t(1:3,:);
        a1 = zeros(size(B1));
        a2 = zeros(size(B1));
        a3 = zeros(size(B1));
        a1 = [ u(t(1,:)); u(t(1,:)) ];
        a2 = [ u(t(2,:)); u(t(2,:)) ];
        a3 = [ u(t(3,:)); u(t(3,:)) ];
        
        % Gradients.
        g = a1.*B1 + a2.*B2 + a3.*B3;
        
    end     % END grad_u().
    
    
    
    function [g_nodes] = grad_nodes( obj, u )
        % GRAD_NODES Returns the gradient of u interpolated at the nodes.
        %   Gradient interpolation is done by first finding all the triangles
        %   associated with a node, then setting the node gradient as the 
        %   average of the triangle gradients weighted by triangle areas.
        
        detB = abs(obj.detB);       % Note abs() for triangle areas.
        mesh = obj.mesh;

        g_tri = obj.grad(u);      
        gw = g_tri .* [detB'; detB'];   % Weighted triangle gradients.
        
        % Sort the triangle nodes based on their occurrence ranks.
        t = mesh.t(1:3,:);
        [a,b] = sort(t(:));
        
        % Indices of first and last occurrences for nodes.
        [~,j] = unique(a);
        i_first = j;
        i_last = [j(2:end)-1; length(a)]; 
        
        g_nodes = zeros( size(mesh.p) );
        for i = 1:length(mesh.p)
            j = b( i_first(i) : i_last(i) );
            j = ceil(j./3);     % Get the indices in 3xN array.
            g_nodes(:,i) = sum(gw(:,j),2) ./ sum(detB(j));
        end       
    end     % END grad_nodes().


    end     % END public methods.
    
end

