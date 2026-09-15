function vtkCellID = vtkCellTypeFromElemType(elemType)
% ============================================================
% vtkCellTypeFromElemType - VTK Topological Dictionary Mapping
%
% PHASE-FIELD VISUALIZATION RELEVANCE:
% Standard finite element codes utilize internal nomenclature to define 
% element topologies (e.g., 'Q4' for a 4-node quadrilateral, 'H8' for an 
% 8-node hexahedron). However, the VTK XML Unstructured Grid format (.vtu) 
% relies on a strict, standardized C++ enumerated dictionary of integer IDs 
% to instruct ParaView on how to render the geometric connections.
%
% This helper function acts as the translation layer between the internal 
% PF-COMPAS mesh struct and the external VTK rendering engine.
%
% Standard VTK Cell Types:
%   VTK_TRIANGLE   = 5
%   VTK_QUAD       = 9
%   VTK_TETRA      = 10
%   VTK_HEXAHEDRON = 12
% ============================================================

% Convert the input string to lowercase to ensure case-insensitive matching
switch lower(elemType)
    case {'q4','quad4'}
        % 4-Node Bilinear Quadrilateral
        vtkCellID = 9;    
        
    case {'t3','tri3'}
        % 3-Node Linear Triangle
        vtkCellID = 5;    
        
    case {'tet4'}
        % 4-Node Linear Tetrahedron
        vtkCellID = 10;   
        
    case {'h8'}
        % 8-Node Linear Hexahedron
        vtkCellID = 12;   
        
    otherwise
        % If a user tries to export an unsupported higher-order element 
        % (like Q16 or H20) without proper VTK high-order mapping, 
        % the export is safely aborted to prevent ParaView rendering crashes.
        error('Unsupported element type for VTK export: %s', elemType);
end

end