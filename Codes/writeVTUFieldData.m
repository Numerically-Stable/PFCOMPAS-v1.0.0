function writeVTUFieldData(fid,sectionName,dataStruct,nEntities)
% ============================================================
% writeVTUFieldData - VTK XML Field Array Generation
%
% PHASE-FIELD HPC & VISUALIZATION RELEVANCE:
% A phase-field simulation exports a mixture of scalar fields (damage phi, 
% strain history H, error estimators Eta) and vector/tensor fields 
% (displacement u, Cauchy stress sigma). 
%
% Instead of hard-coding separate XML write blocks for every physical 
% variable, this subroutine utilizes a dynamic structure parser. It iterates 
% through the provided dataStruct, automatically querying the field dimensions 
% and writing the appropriate <DataArray> block. 
%
% CRITICAL SAFETY FEATURE:
% VTK parsers are notoriously unforgiving. If a field array has 999 entries 
% but the mesh declares 1000 nodes, ParaView will either silently corrupt 
% the visualization or crash entirely. This routine strictly validates 
% the row count of every field against nEntities before writing a single 
% byte to the disk, safeguarding the integrity of the .vtu file.
% ============================================================

%% ------------------------------------------------------------
% Empty data check
% ------------------------------------------------------------
% If no cell data or point data was requested for this step, exit cleanly.
if nargin < 3 || isempty(dataStruct)
    return;
end

%% ------------------------------------------------------------
% Open XML Section (e.g., <PointData> or <CellData>)
% ------------------------------------------------------------
fprintf(fid, '<%s>\n',sectionName);

%% ------------------------------------------------------------
% Dynamic Field Iteration
% ------------------------------------------------------------
names = fieldnames(dataStruct);
for i = 1:numel(names)
    name = names{i};
    data = dataStruct.(name);
    
    %% --------------------------------------------------------
    % Validate Dimensions and Formatting
    % --------------------------------------------------------
    % Force purely 1D arrays (like phi or Eta) into explicit column vectors
    if isvector(data)
        data = data(:);
    end
    
    % STRICT VALIDATION: The number of rows must perfectly match the 
    % number of declared entities (nNodes for PointData, nElem for CellData).
    if size(data,1) ~= nEntities
        error( ...
            ['Field "%s" in %s has %d rows, ' ...
             'but %d entities are expected.'], ...
            name, sectionName, size(data,1), nEntities);
    end
    
    % The number of columns defines the NumberOfComponents for the VTK array.
    % Scalar = 1, Vector = 2 or 3, Tensor (Voigt) = 3 or 6.
    nComp = size(data,2);
    
    %% --------------------------------------------------------
    % Data Type Verification
    % --------------------------------------------------------
    % Prevent the accidental writing of chars or logicals which break 
    % the Float64 ASCII formatting expectation.
    if ~isnumeric(data)
        error('Field "%s" must be numeric.', name);
    end
    
    %% --------------------------------------------------------
    % Write DataArray Header
    % --------------------------------------------------------
    fprintf(fid, ...
        '<DataArray type="Float64" Name="%s" NumberOfComponents="%d" format="ascii">\n', ...
        name,nComp);
        
    %% --------------------------------------------------------
    % Write Values
    % --------------------------------------------------------
    % High-precision ASCII export (16 significant digits) is utilized to 
    % ensure the infinitesimal variations in the phase-field order parameter 
    % (especially near the phi ~ 0 bounds) are preserved for contouring.
    for k = 1:nEntities
        fprintf(fid,'%.16g ',data(k,:));
        fprintf(fid,'\n');
    end
    
    fprintf(fid, '</DataArray>\n');
end

%% ------------------------------------------------------------
% Close XML Section
% ------------------------------------------------------------
fprintf(fid, '</%s>\n',sectionName);

end