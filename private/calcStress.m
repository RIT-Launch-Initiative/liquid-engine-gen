function [FOS, vonMisesStress] = calcStress(transientTemperature,Z,R,axialPressure,~,atmosphericPressure,CTE,poissonsRatio)
%% INITIALIZE
axialPressure = axialPressure';

pressureRadialStress = NaN(size(Z));
pressureHoopStress = NaN(size(Z));

thermalStress = NaN(size(Z));
yieldStrength = NaN(size(Z));
%% 304 SS MATERIAL PROPERTIES
% Yield Strength vs. Temperature
% We want to use the yield strength because we want to maintain the shape
% of the thrust chamber throughout the burn so it could potentially be reused,
% rather than the point where it will rupture at the ultimate tensile
% strength.
yieldStrengthTable = [290 182 150 134 125 113 95 68].*10^6; % [Pa]
yieldStrengthTempTable = [27 149 260 371 482 593 704 816] + 273.15; % [K]
yieldFit = fit(yieldStrengthTempTable',yieldStrengthTable','pchipinterp');

% Check if data exceeds beyond literature material properties table. It is
% not recommended to use the stress values for material properties beyond
% the recorded temperature, for accuracy sake.
if(max(max(transientTemperature(:,:,end),[],2)) > yieldStrengthTempTable(end))
    warning('Maximum wall temperature exceeds maximum material properties temperature recorded in literature, exiting.')
end

% Young's Modulus (tensile)
youngsModulusTable = [193 192 187 183 179 177 170 166 160 155 150 145 141 134 125].*10^6; % [Pa]
youngsModulusTempTable = [27 93 149 204 260 316 371 427 482 538 593 649 704 760 816]+273.15; % [K]
youngsModulusFit = fit(youngsModulusTempTable',youngsModulusTable','pchipinterp');

%% OUTER ENGINE CONTOUR
engineContour = [Z(1,:);R(1,:)];
outerEngineContour = [Z(end,:);R(end,:)];

%% AXIAL STRESS
% Axial Stress is constant, and is non-zero only for closed cylindrical
% pressure vessels. Manually set to zero since this is an "open
% pressure-vessel"
pressureAxialStress = (axialPressure.*engineContour(2,:).^2 - atmosphericPressure.*outerEngineContour(2,:).^2)...
    ./(outerEngineContour(2,:).^2-engineContour(2,:).^2); % [Pa]
pressureAxialStress = 0;
%% HOOP/TANGENTIAL & RADIAL STRESS
for i=1:length(R(:,1))
    iterRadius = R(i,:);
    pressureHoopStress(i,:) = pressureAxialStress + engineContour(2,:).^2.*outerEngineContour(2,:).^2.*(axialPressure-atmosphericPressure)...
        ./(iterRadius.^2.*(outerEngineContour(2,:).^2 - engineContour(2,:).^2));
    pressureRadialStress(i,:) = pressureAxialStress - engineContour(2,:).^2.*outerEngineContour(2,:).^2.*(axialPressure-atmosphericPressure)...
        ./(iterRadius.^2.*(outerEngineContour(2,:).^2 - engineContour(2,:).^2));
end

%% THERMAL STRESS
% Change in temperature at the final timestep relative to the initial
% temperature.
deltaT = transientTemperature(:,:,end)-transientTemperature(:,:,1);

% Thermal Strain is equal in r, z, and theta directions
thermalStrain = CTE.*deltaT;

% Hooke's law
% Thermal stress is equal in r, z, and theta directions
for i=1:length(R(:,1))
    for j=1:length(Z(1,:))
        thermalStress(i,j) = youngsModulusFit(transientTemperature(i,j,end))./((1+poissonsRatio)*(1-2*poissonsRatio)).*...
            ((1-poissonsRatio).*thermalStrain(i,j) + poissonsRatio.*(2.*thermalStrain(i,j)));
    end
end
%% TOTAL STRESS
totalAxialStress = pressureAxialStress + thermalStress;
totalRadialStress = pressureRadialStress + thermalStress;
totalHoopStress = pressureHoopStress + thermalStress;

vonMisesStress = sqrt(0.5.*((totalRadialStress-totalHoopStress).^2+(totalHoopStress-totalAxialStress).^2+(totalAxialStress-totalRadialStress).^2));
%% FACTOR OF SAFETY
for i=1:length(R(:,1))
    for j=1:length(Z(1,:))
        yieldStrength(i,j) = yieldFit(transientTemperature(i,j,end));
    end
end

FOS = yieldStrength./vonMisesStress;


% clf;
% s1 = surf(Z,R,thermalStress./10^6);
% view(0,90); axis equal;
% xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
% s1.EdgeColor = 'none';
% s1.FaceColor = 'interp';
% 
% title('\bf{Thermal Stress}');
% xlabel('Z-Axis $[m]$');ylabel('R-axis $[m]$');
% a1 = colorbar; colormap("turbo");
% ylabel(a1,'Axial Stress $[MPa]$','FontSize',16,'Rotation',270,'Interpreter','latex');
% grid on; grid minor;
% 
% % clf;
% % s1 = surf(Z,R,pressureAxialStress./10^6);
% % view(0,90); axis equal;
% % xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
% % s1.EdgeColor = 'none';
% % s1.FaceColor = 'interp';
% % 
% % title('\bf{Axial Stress}');
% % xlabel('Z-Axis $[m]$');ylabel('R-axis $[m]$');
% % a1 = colorbar; colormap("turbo");
% % ylabel(a1,'Axial Stress $[MPa]$','FontSize',16,'Rotation',270,'Interpreter','latex');
% % grid on; grid minor;
% 
% clf;
% s1 = surf(Z,R,pressureRadialStress./10^6);
% view(0,90); axis equal;
% xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
% s1.EdgeColor = 'none';
% s1.FaceColor = 'interp';
% 
% title('\bf{Radial Stress}');
% xlabel('Z-Axis $[m]$');ylabel('R-axis $[m]$');
% a1 = colorbar; colormap("turbo");
% ylabel(a1,'Radial Stress $[MPa]$','FontSize',16,'Rotation',270,'Interpreter','latex');
% grid on; grid minor;
% 
% clf;
% s1 = surf(Z,R,pressureHoopStress./10^6);
% view(0,90); axis equal;
% xlim([0 Z(1,end)]); ylim([0 1.5*max(R(end,:))]);
% s1.EdgeColor = 'none';
% s1.FaceColor = 'interp';
% 
% title('\bf{Hoop Stress}');
% xlabel('Z-Axis $[m]$');ylabel('R-axis $[m]$');
% a1 = colorbar; colormap("turbo");
% ylabel(a1,'Hoop Stress $[MPa]$','FontSize',16,'Rotation',270,'Interpreter','latex');
% grid on; grid minor;

end