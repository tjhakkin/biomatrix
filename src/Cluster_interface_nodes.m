%
% Given interface edge nodes pairs, forms sorted clusters of node pairs such 
% that disjoint interface loops belong to different clusters.
%

function [clusters] = Cluster_interface_nodes( if_nodes )

    clusters = {};
    ifn = if_nodes;                     % ifn will be modified
    
    while (size(ifn,2) > 0)
        edge = 1;                       % edge currently being processed
        nodes = sort(ifn(:,edge));        
        node0 = nodes(1);               % first node of the current loop

        % Process non-closed interface segments first, if any; these will be 
        % closed (see below) by connecting the orphan nodes to form a new edge.
        [n, nodes] = hist(ifn(:), unique(ifn(:)));
        orphans = find(n == 1);
        if (~isempty(orphans))
            node0 = nodes(orphans(1));
            [~, edge] = find(ifn == node0);
        end

        c = zeros(size(ifn));           % edge node pairs belonging to a cluster
        excl = zeros(size(ifn,2), 1);   % edges to be removed from ifn
        n = 0;                          % number of elements in excl, c
        node = setdiff(ifn(:,edge), node0);

        while (node ~= node0)           % while interface loop not closed
            n = n+1;
            c(:,n) = sort(ifn(:,edge));     % store current edge nodes
            excl(n) = edge;                 % mark edge for removal
            [~,j] = find(ifn == node);  
            edge = setdiff(j, edge);        % next edge associated with node
            if (isempty(edge))
                % Unclosed loop with no connecting edge -> connect the first and
                % the last nodes to form a closed circle.
                n = n+1;
                c(:,n) = [node0, node];
                break;
            end
            node = setdiff(ifn(:,edge), node);
        end
        c = c(:,1:n);
        c(:,end) = flip(c(:,end));
        clusters{end+1} = c;

        excl = nonzeros(excl);
        ifn(:,excl) = [];               % remove current cluster nodes
    end

end
