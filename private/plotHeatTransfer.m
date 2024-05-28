function [] = plotHeatTransfer(transientTemperatureData,Z,R,burnTime,t4)
    Z = convlength(Z,'m','in');
    R = convlength(R,'m','in');

    ax1 = axes('Parent',t4);

    s1 = surf(Z,R,transientTemperatureData(:,:,end));
    view(0,90); axis equal;
    xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
    s1.EdgeColor = 'none'; s1.FaceColor = 'interp';

    title('\bf{2-D Transient Heat Diffusion F.E.A.}');
    xlabel('Z-Axis $[in]$');ylabel('R-axis $[in]$');
    subtitle(sprintf('Burn Time: %.2f $[s]$',burnTime),'FontSize',14)
    a1 = colorbar; colormap(ax1, "hot");
    ylabel(a1,'Temperature $[^\circ{} K]$','FontSize',16,'Rotation',270,'Interpreter','latex');
    %set(ax1,'ColorScale','log')
    grid on; grid minor;
   

end

