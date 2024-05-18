function [] = plotGeometry(engine_contour,t1)
    engine_contour = convlength(engine_contour,'m','in'); % [in]

    axes('Parent',t1);
    hold on
    plot(engine_contour(1,:),engine_contour(2,:),'k');
    plot(engine_contour(1,:),-engine_contour(2,:),'k');
    plot([0 engine_contour(1,end)],[0 0],'-.k');
    axis equal; grid on;
    xlim([engine_contour(1,1),engine_contour(1,end)]); 
    ylim([-1.5*max(engine_contour(2,:)) 1.5*max(engine_contour(2,:))]);

    xlabel('Length $[in]$'); ylabel('Radius $[in]$');
    title('\bf{Thrust Chamber Geometry}')

end

