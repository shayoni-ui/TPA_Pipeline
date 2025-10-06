function f = contourfUnstructured(x, y, z)
    [x1,y1] = meshgrid(linspace(min(x), max(x), 20), ...
        linspace(min(y), max(y), 20));
    z1 = griddata(x,y,z,x1,y1, "cubic");
    f = figure('DefaultAxesFontSize', 16);
    contourf(x1, y1, z1,25,'LineStyle','None');
    colorbar(); colormap('jet');
    %hold on;
    %scatter(x,y,10, 'k', 'filled');
end