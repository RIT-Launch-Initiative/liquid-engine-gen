function [] = plotFilmCoefficient(engine_contour,h_g,T,t3)
    engine_contour_in = convlength(engine_contour,'m','in');

    %filmCoefficientPlot = figure('Name','Film Coefficient');
    axes('Parent',t3);
    hold on
    yyaxis("left")
    plot(engine_contour_in(1,:),h_g,'Color',"#A2142F");
    ylabel('$h_{g}$ $[\frac{W}{m^{2}K}]$')
    grid on
    set(gca,'YColor',"#A2142F")

    yyaxis("right")
    plot(engine_contour_in(1,:),T,'Color',"#D95319")
    ylabel('$T_g$ $[^\circ{}K]$')
    grid on
    set(gca,'YColor',"#D95319")
    
    title('$h_{g}$ and $T_g$ vs. Length')
    xlabel('$L_{e}$ $[in]$');
    xlim([engine_contour_in(1,1),engine_contour_in(1,end)])    
end

