function eta2 = phaseElementResidual2D( ...
    mesh,e,phi,H,material,model,W,Q,nGP,Gc,l0)
% ==============================================================
% One-element phase-field residual contribution.
%
% Returns ONE scalar.
% ==============================================================

conn  = mesh.elem.conn(e,:);
coords = mesh.nodes.coord(conn,:);

phi_e = phi(conn);
H_e   = H(e,:);

residualIntegral = 0;
volume = 0;

for q = 1:nGP

    [N,dNdxi] = ...
        lagrange_basis( ...
        mesh.meta.elemType,Q(q,:));

    J = coords' * dNdxi;

    detJ = det(J);

    if detJ <= 0
        error( ...
            'Non-positive Jacobian: element %d, GP %d.', ...
            e,q);
    end

    phi_gp = N'*phi_e;

    % ----------------------------------------------------------
    % Degradation derivative
    % ----------------------------------------------------------

    [~,dg,~] = ...
        degradationFunction( ...
        phi_gp, ...
        model, ...
        material);

    % ----------------------------------------------------------
    % AT2 residual
    % ----------------------------------------------------------

    Rphi = ...
        dg*H_e(q) ...
        + (Gc/l0)*phi_gp;

    weight = W(q)*detJ;

    residualIntegral = ...
        residualIntegral ...
        + Rphi^2*weight;

    volume = volume + weight;

end

% --------------------------------------------------------------
% Characteristic element size
% --------------------------------------------------------------

hK = sqrt(volume);

eta2 = hK^2*residualIntegral;

end