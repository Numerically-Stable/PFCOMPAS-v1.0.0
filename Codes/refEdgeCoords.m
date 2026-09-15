function [xi, eta] = refEdgeCoords(loc, s)
% ============================================================
% refEdgeCoords - 1D to 2D Reference Domain Mapping
%
% COMPUTATIONAL MECHANICS RELEVANCE:
% When evaluating contour integrals (such as the traction vectors for 
% macroscopic reaction forces or the flux jumps for the Verfürth error 
% estimators), the numerical integration is performed over a 1D edge. 
% 
% The Gauss1D routine provides the integration points 's' strictly in the 
% 1D reference domain [-1, 1]. However, to evaluate the shape functions 
% (N) and their spatial derivatives (B-matrix), we require the 2D parent 
% coordinates (xi, eta) within the bi-unit square [-1, 1] x [-1, 1].
%
% This utility function rigidly maps the 1D coordinate 's' to its exact 
% (xi, eta) location on the perimeter of the Q4 reference element, 
% depending on which local edge ('loc' 1 through 4) is being integrated.
% ============================================================

switch loc
    case 1
        % Bottom Edge (Nodes 1-2): eta is fixed at -1, xi varies with s
        xi = s;  
        eta = -1;
        
    case 2
        % Right Edge (Nodes 2-3): xi is fixed at +1, eta varies with s
        xi = 1;  
        eta = s;
        
    case 3
        % Top Edge (Nodes 3-4): eta is fixed at +1, xi varies with s
        xi = s;  
        eta = 1;
        
    case 4
        % Left Edge (Nodes 4-1): xi is fixed at -1, eta varies with s
        xi = -1; 
        eta = s;
        
    otherwise
        error('Invalid local edge index. Must be 1, 2, 3, or 4.');
end

end