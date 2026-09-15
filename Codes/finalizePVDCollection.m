function finalizePVDCollection(pvdFile)
% ============================================================
% finalizePVDCollection - Close XML Wrapper for ParaView
%
% PHASE-FIELD VISUALIZATION RELEVANCE:
% After the entire quasi-static or dynamic load-stepping timeline has 
% concluded—or if the adaptive solver hits a terminal minimum load step—
% the master .pvd collection must be formally sealed.
%
% VTK files are strictly parsed as XML documents. If the closing tags 
% are missing, many XML parsers (including certain versions of ParaView) 
% will declare the file corrupt and refuse to load the animation, even if 
% all the underlying .vtu files are perfectly intact.
%
% This function safely appends the final closing tags, completing the 
% visualization pipeline and rendering the dataset ready for post-processing.
% ============================================================

fid = fopen(pvdFile,'a');
if fid == -1
    error('Could not open PVD file:\n%s', pvdFile);
end

% ------------------------------------------------------------
% Append Closing XML Tags
% ------------------------------------------------------------
fprintf(fid, '  </Collection>\n');
fprintf(fid, '</VTKFile>\n');

% Safely close the file stream one final time.
fclose(fid);

end