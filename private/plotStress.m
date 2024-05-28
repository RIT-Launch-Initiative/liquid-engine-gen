function [] = plotStress(Z,R,FOS,vonMisesStress,burnTime,t5)
    Z = convlength(Z,'m','in');
    R = convlength(R,'m','in');

    t = tiledlayout(t5,'flow','TileSpacing','tight');
    
    t1 = nexttile(t,[1,2]);
    s1 = surf(Z,R,vonMisesStress./10^6);
    view(0,90); axis equal;
    xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
    s1.EdgeColor = 'none';s1.FaceColor = 'interp';
    
    title('\bf{von Mises Stress}');
    xlabel('Z-Axis $[in]$');ylabel('R-axis $[in]$');
    subtitle(sprintf('Burn Time: %.2f $[s]$',burnTime),'FontSize',14)
    a1 = colorbar; colormap(t1,"turbo");
    ylabel(a1,'von Mises Stress $[MPa]$','FontSize',16,'Rotation',270,'Interpreter','latex');
    grid on; grid minor;
    
    t2 = nexttile(t,[1,2]);
    s2 = surf(Z,R,FOS);
    view(0,90); axis equal;
    xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
    s2.EdgeColor = 'none';
    s2.FaceColor = 'interp';
    
    title('\bf{F.O.S. to Yield Strength}');
    xlabel('Z-Axis $[in]$');ylabel('R-axis $[in]$');
    subtitle(sprintf('Burn Time: %.2f $[s]$',burnTime),'FontSize',14)
    a2 = colorbar; cmap = colormap(t2,"turbo"); colormap(t2,flipud(cmap)); clim([0, 15]);
    ylabel(a2,'F.O.S. $[-]$','FontSize',16,'Rotation',270,'Interpreter','latex');
    grid on; grid minor;

end
