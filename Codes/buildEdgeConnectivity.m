function edgeConnectivity = buildEdgeConnectivity(connectivity)
% buildEdgeConnectivity - Generate global edge connectivity for a Q4 mesh
%
% This function constructs a topological map of the mesh by identifying 
% which elements share which edges. This is a strict prerequisite for 
% computing Verfürth a posteriori error estimators, which require the 
% evaluation of jump residuals (e.g., traction or flux differences) 
% strictly across internal element boundaries.
%
% INPUT:
%   connectivity : [nelem x 4] global element connectivity matrix
%
% OUTPUT:
%   edgeConnectivity : [nedges x 4] array mapping global edges to elements.
%       Columns: [elemL, localEdgeL, elemR, localEdgeR]
%       locL, locR in {1,2,3,4} following the standard Q4 convention.
%       For external boundary edges, elemR and localEdgeR are set to 0.
%
% Q4 local edge convention:
%   edge 1: nodes 1-2
%   edge 2: nodes 2-3
%   edge 3: nodes 3-4
%   edge 4: nodes 4-1

nelem = size(connectivity,1);

% -------------------------------------------------------------------
% Hash-Map Initialization
% A dictionary/hash-map is utilized to achieve near O(N) algorithmic 
% complexity when searching for shared edges, avoiding the prohibitive 
% O(N^2) cost of brute-force global matrix searches.
% -------------------------------------------------------------------
edgeMap = containers.Map('KeyType','char','ValueType','any');

% Local edge nodes for Q4 isoparametric elements
localEdgeNodes = [1 2; 2 3; 3 4; 4 1];

for e = 1:nelem
    for loc = 1:4
        % Extract the global node IDs for the current local edge
        nodes = connectivity(e, localEdgeNodes(loc,:));
        
        % CRITICAL: Sort the nodes to ensure direction-independence.
        % Edge A-B must generate the exact same string key as Edge B-A 
        % so they are correctly identified as the same shared interface.
        nodes = sort(nodes); 
        key = sprintf('%d_%d', nodes(1), nodes(2));
        
        if isKey(edgeMap,key)
            % If the key exists, this is an interior edge shared by 2 elements.
            % Append the current element (e) and its local edge ID (loc).
            edgeMap(key) = [edgeMap(key); e loc]; 
        else
            % If the key does not exist, initialize it with the first element found.
            edgeMap(key) = [e loc];
        end
    end
end

% -------------------------------------------------------------------
% Construct the Output Array (edgeConnectivity)
% Iterate through the hash-map to classify edges as interior or boundary.
% -------------------------------------------------------------------
edgeConnectivity = [];
keysEdge = keys(edgeMap);

for k = 1:length(keysEdge)
    val = edgeMap(keysEdge{k});
    
    if size(val,1) == 1
        % Boundary Edge: Only one element contains this node pair.
        % Pad the right-side element and local edge variables with zeros.
        edgeConnectivity = [edgeConnectivity; val(1) val(2) 0 0];
        
    elseif size(val,1) == 2
        % Interior Edge: Two adjacent elements share this node pair.
        % Map the left (L) and right (R) elements for jump evaluations.
        eL = val(1,1); locL = val(1,2);
        eR = val(2,1); locR = val(2,2);
        edgeConnectivity = [edgeConnectivity; eL locL eR locR];
        
    else
        % Topological Sanity Check: A 2D manifold cannot have more than 
        % two elements sharing a single edge. 
        error('Topological error: More than 2 elements share the same edge. Invalid mesh.');
    end
end
end