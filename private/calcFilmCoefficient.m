function [h_film_x] = calcFilmCoefficient(engine_contour,Rc,mdot,visc,cp_tran,prantyl)
%% TITLE: Calculating Convective Film Coefficient In Terms of Engine Length
% Author: Evan Olsen
% Date: 03/26/2024
% Description: Uses pre-calculated interior contour and CEA transport
% properties to calculate heat transfer coefficient, according to the
% equation described in source: https://jffhmt.avestia.com/2015/PDF/006.pdf

%% INPUT DEFINITION

% Lateral surface area of the interior geometry revolved about the z-axis.

f_Eng = @(x) interp1(engine_contour(1,:), engine_contour(2,:), x, 'pchip','extrap');
a = min(engine_contour(1,:));
b = max(engine_contour(1,:));

A_lat = 2*pi*integral(@(x) f_Eng(x), a, b);

% Calculates the local cross-sectional area of the engine based on the
% local radius.

A_x = pi.*engine_contour(2,:).^2; % [m^2]

%%  Eq. Before Throat, Fr. After Throat
% Found to be the most accurate method to produce heat transfer film
% coefficient: Uses equillibrium transport values prior to throat, and
% frozen transport values following the throat.
%
% Linearly interpolates between each of the values to give the best
% interpretation of the values between each point.
temp = engine_contour(2,:)-engine_contour(2,1);
for i=1:length(engine_contour(2,:))
    if(temp(i) < 0)
        indce = i;
        break
    end
end
[~,indt] = min(engine_contour(2,:));
indne = length(engine_contour(1,:));
INJ2CE = indce;
CE2TH = indt-indce;
TH2NE = indne-indt;

% C_p [J/(kg*K)] --- Specific Heat Constant Pressure
Cp_x = [linspace(cp_tran.eql(1),cp_tran.eql(2),INJ2CE),...
    linspace(cp_tran.eql(2),cp_tran.eql(3),CE2TH),...
    linspace(cp_tran.froz(3),cp_tran.froz(4),TH2NE)].*1000;

% Pr [-] --- Prandtl Number
Pr_x = [linspace(prantyl.eql(1),prantyl.eql(2),INJ2CE),...
    linspace(prantyl.eql(2),prantyl.eql(3),CE2TH),...
    linspace(prantyl.froz(3),prantyl.froz(4),TH2NE)];

% Dynamic Viscosity as provided by NASA CEA as a special output.
% Ensure you convert cP to Pa*s before running.
visc_x = [linspace(visc(1),visc(2),INJ2CE),...
    linspace(visc(2),visc(3),CE2TH),...
    linspace(visc(3),visc(4),TH2NE)]; % [Pa*s]
%% CALCULATION
% Solving equation as is described in https://jffhmt.avestia.com/2015/PDF/006.pdf

% Correction factor based on combustion chamber radius.
Zc= pi*Rc^2/A_lat; % [-]
h_film_x = Zc.*(mdot./(2.*A_x)).*Cp_x.*visc_x.^(0.3).*Pr_x.^(2/3); % [W/(m^2*K)]
%% ANSYS FORMAT OUTPUT
% Inverses the x-values due to ANSYS coordinates being inversed. Change
% this value as represents your coordinate system in ANSYS.
% Creates vertical columns for easy import into ANSYS Transient Thermal
% Analysis.
output_matrix = [engine_contour(1,:);h_film_x]'; % [m];[W/(m^2*K)]

outputFolder = 'output';
filePath = fullfile(outputFolder,'filmCoefficient.txt');
writematrix(output_matrix,filePath)
% Outputs file path where the 'contour.csv' file was generated
fprintf('''filmCoefficient.txt'' generated at %s\n',filePath);

end

