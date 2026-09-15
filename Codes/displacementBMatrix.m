function B = displacementBMatrix(dNdx,dim,nen)
% ============================================================
% displacementBMatrix - Construct FEM Strain-Displacement Operator
%
% COMPUTATIONAL MECHANICS RELEVANCE:
% The infinitesimal strain tensor is mathematically a symmetric 3x3 matrix.
% However, operating on 4th-order constitutive tensors and 2nd-order 
% strain/stress tensors during numerical integration is computationally 
% inefficient. Finite element codes utilize Voigt notation to map these 
% symmetric tensors into 1D vectors and 2D matrices.
%
% This function constructs the global B-matrix which relates the nodal 
% displacement vector (u_e) to the Voigt strain vector (strain_v):
%       strain_v = B * u_e
%
% IMPORTANT - ENGINEERING SHEAR STRAIN:
% This implementation strictly utilizes the engineering shear strain 
% convention (gamma_xy = 2 * epsilon_xy = du/dy + dv/dx). Any constitutive 
% matrix (C) multiplied by this B-matrix must be formulated to expect 
% engineering shear strains to ensure energy equivalence.
% ============================================================

if dim == 2
    % --------------------------------------------------------
    % 2D B-Matrix Assembly (Plane Strain / Plane Stress)
    % Voigt Ordering: [epsilon_xx; epsilon_yy; gamma_xy]
    % DOFs: [u1, v1, u2, v2, ..., un, vn]
    % --------------------------------------------------------
    B = zeros(3,2*nen);
    
    % Row 1: epsilon_xx = du/dx
    B(1,1:2:end) = dNdx(:,1)';
    
    % Row 2: epsilon_yy = dv/dy
    B(2,2:2:end) = dNdx(:,2)';
    
    % Row 3: gamma_xy = du/dy + dv/dx
    B(3,1:2:end) = dNdx(:,2)'; % du/dy
    B(3,2:2:end) = dNdx(:,1)'; % dv/dx

elseif dim == 3
    % --------------------------------------------------------
    % 3D B-Matrix Assembly
    % Voigt Ordering: [eps_xx; eps_yy; eps_zz; gamma_xy; gamma_yz; gamma_xz]
    % DOFs: [u1, v1, w1, u2, v2, w2, ..., un, vn, wn]
    % --------------------------------------------------------
    B = zeros(6,3*nen);
    
    % Normal Strains
    B(1,1:3:end) = dNdx(:,1)'; % epsilon_xx = du/dx
    B(2,2:3:end) = dNdx(:,2)'; % epsilon_yy = dv/dy
    B(3,3:3:end) = dNdx(:,3)'; % epsilon_zz = dw/dz
    
    % Engineering Shear Strains
    % gamma_xy = du/dy + dv/dx
    B(4,1:3:end) = dNdx(:,2)'; 
    B(4,2:3:end) = dNdx(:,1)'; 
    
    % gamma_yz = dv/dz + dw/dy
    B(5,2:3:end) = dNdx(:,3)'; 
    B(5,3:3:end) = dNdx(:,2)'; 
    
    % gamma_xz = du/dz + dw/dx
    B(6,1:3:end) = dNdx(:,3)'; 
    B(6,3:3:end) = dNdx(:,1)';
else
    error('Unsupported dimension. Only 2D and 3D are supported.');
end

end