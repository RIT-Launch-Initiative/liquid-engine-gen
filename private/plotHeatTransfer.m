function [animationObject] = plotHeatTransfer(transientTemperatureData,Z,R,burnTime,dt,t4)
    axes('Parent',t4);

    totalFrames = length(transientTemperatureData(1,1,:));

    targetFPS = 60;
    totalRequiredFrames = targetFPS*burnTime;
    iterSkip = round(totalFrames/totalRequiredFrames);

    hold on
    title('\bf{2D Transient Heat Transfer FEA}');
    xlabel('Z-Axis $[m]$');ylabel('R-axis $[m]$');
    subtitle(sprintf('Burn Time: %.2f $[s]$',burnTime),'FontSize',14)
    a = colorbar;
    colormap(hot)
    ylabel(a,'Temperature $[^\circ{} K]$','FontSize',16,'Rotation',270,'Interpreter','latex');
    set(gca,'ColorScale','log')
    grid on; grid minor;

    s = surf(Z,R,transientTemperatureData(:,:,end));axis equal; 
    xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';

    endColormap = colormap;

    j = 1;
    for i=1:iterSkip:length(transientTemperatureData)
        s = surf(Z,R,transientTemperatureData(:,:,i));axis equal; 
        xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
        colormap(endColormap)
        s.EdgeColor = 'none';
        s.FaceColor = 'interp';

        clc;
        fprintf('frame %.0f of %.0f\n',i/iterSkip,totalFrames/iterSkip)
        animationObject(j) = getframe();
        j = j+1;
    end
    
end

