function U3 = reshapeDisplacementToNodal(u,dim,nnodes)

% ============================================================
% RESHAPE GLOBAL DISPLACEMENT VECTOR TO NODAL FORMAT
%
% GLOBAL DOF ORDERING
%
% 2D:
%
%   u = [ux1 uy1 ux2 uy2 ... uxn uyn]'
%
% 3D:
%
%   u = [ux1 uy1 uz1 ux2 uy2 uz2 ... uxn uyn uzn]'
%
%
% OUTPUT
%
%   U3(:,1) = ux
%   U3(:,2) = uy
%   U3(:,3) = uz
%
% Size:
%
%   U3 = nnodes x 3
%
% ============================================================


%% ------------------------------------------------------------
% Input checks
% ------------------------------------------------------------

if ~ismember(dim,[2 3])

    error('dim must be either 2 or 3.');

end


u = u(:);


%% ------------------------------------------------------------
% Check displacement vector size
% ------------------------------------------------------------

if numel(u) ~= dim*nnodes

    error( ...
        ['Invalid displacement vector size.\n' ...
         'Expected %d DOFs, received %d.'], ...
        dim*nnodes, ...
        numel(u));

end


%% ------------------------------------------------------------
% Initialize
% ------------------------------------------------------------

U3 = zeros(nnodes,3);


%% ============================================================
% 2D
% ============================================================

if dim == 2

    U3(:,1) = u(1:2:end);   % ux
    U3(:,2) = u(2:2:end);   % uy


%% ============================================================
% 3D
% ============================================================

elseif dim == 3

    U3(:,1) = u(1:3:end);   % ux
    U3(:,2) = u(2:3:end);   % uy
    U3(:,3) = u(3:3:end);   % uz

end

end