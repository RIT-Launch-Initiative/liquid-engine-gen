function [] = plotMesh(interiorEngineContour,Z,R,Zeta,Eta,t3)
    interiorEngineContour = convlength(interiorEngineContour,'m','in');
    Z = convlength(Z,'m','in');
    R = convlength(R,'m','in');

    t = tiledlayout(t3,'flow','TileSpacing','tight');

    nexttile(t,[1,2])
    V = ones(length(Z(:,1)),length(Z(1,:)));
    p = pcolor(Z,R,V); p.FaceColor = '#ADD8E6';
    axis equal; 
    hold on; plot(interiorEngineContour(1,:),interiorEngineContour(2,:),'-r')
    xlabel('Z-Axis $[in]$');ylabel('R-Axis $[in]$');
    title('Boundary Mesh -- Physical Domain')
    axis equal
    xlim([0 Z(end,end)]); ylim([0 1.5*max(R(1,:))]);
    grid on; grid minor;
    legend('Physical Mesh','Interior Contour')

    nexttile(t,[1,2])
    computationalInteriorEngineContour = [Zeta(1,:);Eta(1,:)];
    
    V = ones(length(Zeta(:,1)),length(Zeta(1,:)));
    p = pcolor(Zeta,Eta,V); p.FaceColor = '#ADD8E6';
    axis equal;
    hold on; plot(computationalInteriorEngineContour(1,:),computationalInteriorEngineContour(2,:),'-r')
    xlabel('Z-Axis $[-]$');ylabel('R-Axis $[-]$');
    title('Boundary Mesh -- Computational Domain')
    axis equal
    
    grid on; grid minor;
    legend('Computational Mesh','Computational Interior Contour')

end