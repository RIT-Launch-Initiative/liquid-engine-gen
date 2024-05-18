function [] = plotRelationships(engine_contour,M,P,T,rho,t2)
    engine_contour_in = convlength(engine_contour,'m','in');

    %relationshipsPlot = figure('Name','Fluid-Area Relationships');
    axes('Parent',t2);
    
    subplot(2,2,1)
    plot(engine_contour_in(1,:),M,'Color','#023E8A')
    xlabel('$L_{e}$ $[in]$');ylabel('$M$ $[-]$');
    set(gca,'YColor','#023E8A')
    title('Mach Number vs. Length');
    xlim([engine_contour_in(1,1) engine_contour_in(1,end)])
    ylim([0 max(M)*1.2])
    grid on

    subplot(2,2,2)
    plot(engine_contour_in(1,:),P,'Color','#08A045')   
    xlabel('$L_{e}$ $[in]$');ylabel('$P$ $[Pa]$');
    set(gca,'YColor','#08A045')
    title('Pressure vs. Length')
    xlim([engine_contour_in(1,1) engine_contour_in(1,end)]);
    ylim([0 max(P)*1.2])
    grid on

    subplot(2,2,3)
    plot(engine_contour_in(1,:),T,'Color','#B53737');
    xlabel('$L_{e}$ $[in]$');ylabel('$T$ $[K]$')
    set(gca,'YColor','#B53737')
    title('Temperature vs. Length')
    xlim([engine_contour_in(1,1) engine_contour_in(1,end)])
    ylim([0 max(T)*1.2])
    grid on

    subplot(2,2,4);
    plot(engine_contour_in(1,:),rho,'Color','#A45729');
    xlabel('$L_{e}$ $[in]$');ylabel('$\rho$ $[\frac{kg}{m^3}]$')
    set(gca,'YColor','#A45729')
    title('Density vs. Length')
    xlim([engine_contour_in(1,1) engine_contour_in(1,end)])
    ylim([0 max(rho)*1.2])
    grid on

    sgtitle('\bf{Fluid-Area Relationships}')

end

