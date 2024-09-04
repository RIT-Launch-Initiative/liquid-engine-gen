function [h_g] = calcBartzCorrelation(wallTemp, combustionStartAxialLocation, engine_contour,throatDiameter,throatArea,prandtlNumber,dynamicViscosity,specificHeat,specificHeatRatio,chamberPressure,characteristicVelocity,localMachNumber,flameTemp)
%% TITLE: Calculating Convective Film Coefficient In Terms of Engine Length
% Author: Evan Olsen
% Date: 08/18/2024
% Description: 

%% ASSUMPTIONS
% All determined at chamber conditions aside from Mach number, and Local
% Area.
%% ANALYSIS
% Determining full-length heat transfer curve, assuming combustion @ x=0
localArea = pi*engine_contour(2,:).^2;
throatCurvature = .382*throatDiameter/2;
specificHeat = specificHeat*1000; % [J/kg-K]

m = 0.6; % [-] Exponent of Viscosity Dependance on Temperature

term1 = 0.026/(throatDiameter^(0.2));
term2 = (dynamicViscosity^(0.2)*specificHeat)/(prandtlNumber^(0.6));
term3 = (chamberPressure/characteristicVelocity)^(0.8);
term4 = (throatArea./localArea).^(0.9);
term5 = (throatDiameter/throatCurvature)^(0.1);

term6 = 0.5*(wallTemp/flameTemp) * (1+((specificHeatRatio-1)/2).*(localMachNumber.^2))+ 0.5;
term7 = 1 + ((specificHeatRatio-1)/2).* localMachNumber.^2;

sigma = 1./(term6.^(0.8-0.2*m).*term7.^(0.2*m));
h_g = term1.*term2.*term3.*term4.*term5.*sigma;

% Correcting curve based on combustion start location
[~,ix] = min(abs(engine_contour(1,:)-combustionStartAxialLocation));
h_g(1:ix) = 28; % [W/m^2-K]

%% ANSYS FORMAT OUTPUT
% Inverses the x-values due to ANSYS coordinates being inversed. Change
% this value as represents your coordinate system in ANSYS.
% Creates vertical columns for easy import into ANSYS Transient Thermal
% Analysis.
output_matrix = [engine_contour(1,:);h_g]'; % [m];[W/(m^2*K)]

outputFolder = 'output';
filePath = fullfile(outputFolder,'filmCoefficient.txt');
writematrix(output_matrix,filePath)
% Outputs file path where the 'contour.csv' file was generated
fprintf('''filmCoefficient.txt'' generated at %s\n',filePath);

end

