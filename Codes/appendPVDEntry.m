function appendPVDEntry(pvdFile,vtuFile,timeValue)
% ============================================================
% appendPVDEntry - Link Spatial VTU Data to Pseudo-Time
%
% PHASE-FIELD VISUALIZATION RELEVANCE:
% This function dynamically injects the latest converged spatial state 
% (the .vtu file) into the master .pvd collection, binding the discrete 
% finite element mesh to a continuous physical parameter (timeValue). 
% In the context of the PF-COMPAS solver, timeValue typically represents 
% the macroscopic load parameter (lambda) or the prescribed displacement.
%
% DATA PORTABILITY (CRITICAL FEATURE):
% By extracting only the relative filename rather than the full absolute 
% system path (e.g., 'C:\Users\...\results\file.vtu'), the resulting 
% output folder becomes completely portable. Researchers can generate 
% the data on a remote Linux HPC cluster and transfer the folder to a 
% local Windows/Mac workstation. Because the .pvd file relies strictly 
% on relative paths, ParaView will load the animation flawlessly.
% ============================================================

% ------------------------------------------------------------
% Open file in Append Mode ('a')
% ------------------------------------------------------------
% Safely open the existing .pvd file without overwriting the header 
% or previous timestep entries.
fid = fopen(pvdFile,'a');
if fid == -1
    error('Could not open PVD file:\n%s', pvdFile);
end

%% ------------------------------------------------------------
% Path Stripping for Dataset Portability
% ------------------------------------------------------------
% Extract the filename and extension, discarding the parent directories.
[~,name,ext] = fileparts(vtuFile);
relativeName = [name ext];

%% ------------------------------------------------------------
% Append XML DataSet Entry
% ------------------------------------------------------------
% The 'timestep' attribute binds the floating-point load value to the frame.
% ParaView interpolates between these timestep values if requested.
fprintf(fid, ...
    '    <DataSet timestep="%.16g" group="" part="0" file="%s"/>\n', ...
    timeValue, ...
    relativeName);

% Immediately close the file stream to protect against unexpected solver crashes.
fclose(fid);

end