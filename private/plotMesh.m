function [] = plotMesh(engine_contour,Z,R,t4)

    axes('Parent',t4)
    V = ones(length(Z(:,1)),length(Z(1,:)));
    p = pcolor(Z,R,V); p.FaceColor = '#ADD8E6';
    axis equal; 
    hold on; plot(engine_contour(1,:),engine_contour(2,:),'-r')
    xlabel('Z-Axis $[m]$');ylabel('R-Axis $[m]$');
    title('\bf{Mesh Boundaries -- Physical Domain}')
    axis equal
    xlim([0 Z(end,end)]); ylim([0 1.5*max(R(1,:))]);
    grid on; grid minor;
    legend('Mesh','Interior Contour')

end