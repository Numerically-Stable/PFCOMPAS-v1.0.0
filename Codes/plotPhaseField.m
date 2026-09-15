function plotPhaseField(mesh, phi)
% ------------------------------------------------------------
% Plot nodal phase field on the associated finite-element mesh
%
% Inputs:
%   mesh : mesh structure
%   phi  : nodal phase-field vector
%
% ------------------------------------------------------------

    figure('Color','w');

    patch( ...
        'Faces',    mesh.elem.conn, ...
        'Vertices', mesh.nodes.coord, ...
        'FaceVertexCData', phi(:), ...
        'FaceColor', 'interp', ...
        'EdgeColor', [0.5 0.5 0.5]);

    axis equal tight;
    view(2);

    xlabel('x');
    ylabel('y');

    title('\phi - Initial Phase Field');

    colormap(jet);
    colorbar;

    clim([0 1]);

    set(gca,'FontSize',12);
    box on;

end