function sigmaAvg = postprocessElementStress( ...
    mesh,u,phi,material,model,quad,parallel)

nelem = mesh.elem.n;

sigmaAvg = zeros(nelem,3);

if parallel.useParpool

    parfor e = 1:nelem

        sigmaAvg(e,:) = ...
            elementAverageStress( ...
            mesh,e,u,phi, ...
            material,model,quad);

    end

else

    for e = 1:nelem

        sigmaAvg(e,:) = ...
            elementAverageStress( ...
            mesh,e,u,phi, ...
            material,model,quad);

    end

end

end