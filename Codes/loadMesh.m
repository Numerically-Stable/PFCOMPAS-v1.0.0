function mesh = loadMesh(filename)
    S = load(filename);
    if ~isfield(S,'mesh')
        error('File does not contain mesh struct');
    end
    mesh = S.mesh;
end