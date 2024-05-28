function [] = plotGeometry(interiorEngineContour,exteriorEngineContour,t1)

    interiorEngineContour = convlength(interiorEngineContour,'m','in'); % [in]
    exteriorEngineContour = convlength(exteriorEngineContour,'m','in'); % [in]

    rightEdge = [interiorEngineContour(1,end) exteriorEngineContour(1,end);...
                 interiorEngineContour(2,end) exteriorEngineContour(2,end)];
    leftEdge =  [interiorEngineContour(1,1) exteriorEngineContour(1,1);...
                 interiorEngineContour(2,1) exteriorEngineContour(2,1)];

    upperWall = [interiorEngineContour rightEdge fliplr(exteriorEngineContour) leftEdge];
    lowerWall = [upperWall(1,:);-upperWall(2,:)];

    axes('Parent',t1);
    upperWallPatch = patch('XData',upperWall(1,:),'YData',upperWall(2,:),'EdgeColor','k');
    hold on;
    lowerWallPatch = patch('XData',lowerWall(1,:),'YData',lowerWall(2,:),'EdgeColor','k');
    hatchfill(upperWallPatch,'single',60,4);
    hatchfill(lowerWallPatch,'single',-60,4);
    
    plot([0 interiorEngineContour(1,end)],[0 0],'-.k','LineWidth',1); hold off;
    axis equal; grid on; grid minor;
    xlim([interiorEngineContour(1,1),interiorEngineContour(1,end)]); 
    ylim([-1.5*max(interiorEngineContour(2,:)) 1.5*max(interiorEngineContour(2,:))]);

    title('\bf{Thrust Chamber Geometry}');
    xlabel('Z-Axis $[in]$'); ylabel('R-Axis $[in]$');
end

