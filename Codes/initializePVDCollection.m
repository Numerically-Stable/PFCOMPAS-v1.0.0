function pvdFile = initializePVDCollection(filename)
% ============================================================
% initializePVDCollection - ParaView Data (PVD) Header Initialization
%
% PHASE-FIELD VISUALIZATION RELEVANCE:
% Phase-field fracture is inherently a time-dependent (or pseudo-time 
% dependent) dissipative process. The simulation outputs a vast time-series 
% of individual .vtu files. Manually loading hundreds of .vtu files into 
% ParaView is tedious and completely loses the associated load parameter 
% (e.g., lambda or displacement delta) for each frame.
%
% A .pvd (ParaView Data) file acts as a "Master Index" or XML wrapper. 
% It binds a specific physical pseudo-time value to a specific .vtu file. 
% When a user opens the single .pvd file in ParaView, the software seamlessly 
% loads the entire timeline, allowing the user to animate the crack propagation 
% mapped correctly to the macroscopic load steps.
%
% HPC I/O RELEVANCE:
% Note that this function opens the file, writes the header, and immediately 
% CLOSES the file identifier (fid). Keeping I/O streams open during a 
% multi-day HPC simulation is highly dangerous; if the cluster crashes or 
% the power fails, an open file stream will corrupt the entire collection. 
% By writing incrementally and closing immediately, the framework is robust 
% against unexpected interruptions.
% ============================================================

fid = fopen(filename,'w');
if fid == -1
    error('Could not create PVD file:\n%s', filename);
end

% ------------------------------------------------------------
% Write XML Header for VTK Collection
% ------------------------------------------------------------
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, '  <Collection>\n');

% ------------------------------------------------------------
% Safely Close File Stream
% ------------------------------------------------------------
% The <Collection> and <VTKFile> tags are intentionally left open. 
% Subsequent converged load steps will append <DataSet> entries to this file. 
% At the end of the simulation, finalizePVDCollection.m will write the 
% closing XML tags.
fclose(fid);

pvdFile = filename;

end