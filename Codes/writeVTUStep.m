function writeVTUStep(filename,mesh,pointData,cellData)
% ============================================================
% writeVTUStep - ParaView/VTK Unstructured Grid Export
%
% PHASE-FIELD HPC & VISUALIZATION RELEVANCE:
% Simulating phase-field fracture naturally generates massive arrays of data 
% across multiple physics fields: displacement vectors (u), damage scalars (phi), 
% stress tensors (sigma), and error estimators (Eta). Rendering highly refined 
% 2D or 3D meshes natively in MATLAB via patch() or trisurf() is notoriously 
% slow and severely bottlenecks the post-processing phase.
%
% To guarantee scalable, high-fidelity visualization, the PF-COMPAS framework 
% offloads all rendering to ParaView by exporting converged states directly 
% into the VTK XML Unstructured Grid (.vtu) format.
%
% MEMORY MANAGEMENT PARADIGM:
% This subroutine is intentionally "timestep-local". It writes strictly the 
% current converged state directly to the hard drive. By doing so, the solver 
% avoids holding entire, massive time-history arrays in RAM, preventing 
% out-of-memory crashes during simulations with thousands of adaptive load steps.
% ============================================================

%% ============================================================
% GEOMETRY EXTRACTION
% ============================================================
nodes = mesh.nodes.coord;
conn  = mesh.elem.conn;
nNodes = mesh.nodes.n;
nElem  = mesh.elem.n;
nen = size(conn,2);

%% ------------------------------------------------------------
% Convert 2D coordinates to 3D VTK coordinates
% ------------------------------------------------------------
% The VTK standard strictly assumes a 3D Cartesian space for spatial points, 
% even for purely 2D plane-strain simulations. We must pad the z-coordinate 
% with zeros to ensure ParaView parses the geometry correctly.
if size(nodes,2) == 2
    nodes3 = [nodes,zeros(nNodes,1)];
elseif size(nodes,2) == 3
    nodes3 = nodes;
else
    error('Unsupported coordinate dimension: %d.', size(nodes,2));
end

%% ============================================================
% VTK CELL TYPE MAPPING
% ============================================================
% Converts FEM topology nomenclature (e.g., 'Q4', 'H8') into standard 
% VTK integer cell types (e.g., 9 for VTK_QUAD, 12 for VTK_HEXAHEDRON).
vtkCellID = vtkCellTypeFromElemType(mesh.meta.elemType);

%% ============================================================
% OPEN FILE WITH SAFE CLEANUP
% ============================================================
fid = fopen(filename,'w');
if fid == -1
    error('Could not open VTU file for writing:\n%s', filename);
end

% The onCleanup object ensures that even if the code crashes or the user 
% aborts mid-write, the file stream is safely closed, preventing file corruption.
cleanupObj = onCleanup(@() fclose(fid));

%% ============================================================
% XML HEADER DEFINITION
% ============================================================
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, '<UnstructuredGrid>\n');
fprintf(fid, '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n', nNodes,nElem);

%% ============================================================
% POINTS (Node Coordinates)
% ============================================================
fprintf(fid,'<Points>\n');
fprintf(fid, '<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n');
for i = 1:nNodes
    fprintf(fid, '%.16g %.16g %.16g\n', nodes3(i,1), nodes3(i,2), nodes3(i,3));
end
fprintf(fid,'</DataArray>\n');
fprintf(fid,'</Points>\n');

%% ============================================================
% CELLS (Element Connectivity and Topology)
% ============================================================
% CRITICAL FIX: Zero-Based Indexing
% MATLAB utilizes 1-based indexing arrays. However, the VTK format is built 
% upon C/C++ architecture, which strictly mandates 0-based indexing for 
% all node references in the connectivity array.
conn0 = conn - 1;

offsets = (1:nElem)*nen;
types = vtkCellID*ones(nElem,1);

fprintf(fid,'<Cells>\n');

%% ------------------------------------------------------------
% Connectivity
% ------------------------------------------------------------
fprintf(fid, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');
for e = 1:nElem
    fprintf(fid,'%d ',conn0(e,:));
    fprintf(fid,'\n');
end
fprintf(fid,'</DataArray>\n');

%% ------------------------------------------------------------
% Offsets
% ------------------------------------------------------------
fprintf(fid, '<DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(fid,'%d ',offsets);
fprintf(fid,'\n</DataArray>\n');

%% ------------------------------------------------------------
% Cell types
% ------------------------------------------------------------
fprintf(fid, '<DataArray type="UInt8" Name="types" format="ascii">\n');
fprintf(fid,'%d ',types);
fprintf(fid,'\n</DataArray>\n');
fprintf(fid,'</Cells>\n');

%% ============================================================
% POINT DATA (Nodal Fields: u, phi)
% ============================================================
% PointData handles continuously interpolated fields at the nodes.
% The displacement vector and phase-field order parameter naturally reside here.
writeVTUFieldData(fid, 'PointData', pointData, nNodes);

%% ============================================================
% CELL DATA (Elemental Fields: sigma, Eta)
% ============================================================
% CellData handles piecewise constant/averaged fields evaluated over the element.
% Stresses evaluated at Gauss points and the Verfürth discretization error 
% estimators are stored here to visualize domain discontinuities.
writeVTUFieldData(fid, 'CellData', cellData, nElem);

%% ============================================================
% CLOSE XML
% ============================================================
fprintf(fid, '</Piece>\n');
fprintf(fid, '</UnstructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');

end