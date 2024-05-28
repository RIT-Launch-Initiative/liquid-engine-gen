function [] = plotProperties(engine_contour,M,P,T,rho,h_g,t2)
    % Unit conversions
    engineContour_in = convlength(engine_contour,'m','in');
    P = P./10^6; % [MPa]
    h_g = h_g./1000; % [kW/m^2-K]

    t = tiledlayout(t2,'flow','TileSpacing','compact');
    
    nexttile(t)
    plot(engineContour_in(1,:),M,'Color','#023E8A','LineWidth',1)
    xlabel('Z-Axis $[in]$');ylabel('$M$ $[-]$');
    set(gca,'YColor','#023E8A')
    title('Mach Number vs. Length');
    xlim([engineContour_in(1,1) engineContour_in(1,end)])
    ylim([0 max(M)*1.2])
    grid on; grid minor;

    nexttile(t)
    plot(engineContour_in(1,:),P,'Color','#08A045','LineWidth',1)   
    xlabel('Z-Axis $[in]$');ylabel('$P$ $[MPa]$');
    set(gca,'YColor','#08A045')
    title('Pressure vs. Length')
    xlim([engineContour_in(1,1) engineContour_in(1,end)]);
    ylim([0 max(P)*1.2])
    grid on; grid minor;

    nexttile(t)
    plot(engineContour_in(1,:),rho,'Color','#A45729','LineWidth',1);
    xlabel('Z-Axis $[in]$');ylabel('$\rho$ $[\frac{kg}{m^3}]$')
    set(gca,'YColor','#A45729')
    title('Density vs. Length')
    xlim([engineContour_in(1,1) engineContour_in(1,end)])
    ylim([0 max(rho)*1.2])
    grid on; grid minor;

    nexttile(t,[1,3])
    yyaxis('left')
    plot(engineContour_in(1,:),T,'Color','#B53737','LineWidth',1);
    ylabel('$T$ $[K]$'); ylim([0 max(T)*1.2])
    set(gca,'YColor','#B53737')
    yyaxis('right')
     plot(engineContour_in(1,:),h_g,'Color',"#0C374D",'LineWidth',1);
    ylabel('$h_{film}$ $[\frac{kW}{m^{2}K}]$')
    grid on
    set(gca,'YColor',"#0C374D")
    title('Temperature \& Convection Coefficient vs. Length'); xlabel('Z-Axis $[in]$');
    xlim([engineContour_in(1,1) engineContour_in(1,end)])
    grid on; grid minor;

end

