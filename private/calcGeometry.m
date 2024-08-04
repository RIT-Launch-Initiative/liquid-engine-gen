%% Function Declaration
function [interiorEngineContour,exteriorEngineContour] = calcGeometry(specificHeatRatio,exitMach,chamberRadius,exitRadius,throatRadius,chamberLength,nozzleType,contourPoints,interiorContractionAngle,exteriorContractionAngle,chamberThickness,nozzleThickness,nozzleFilet)
%% NOZZLE INPUT SELECTION
if nozzleType == 1 % Minimum Length Conical Nozzle
    throatDSFilet = 1.5*throatRadius;
    % 1.5*Rt is used exclusively for conical nozzles in common application.
    % Unlike other nozzle types which have much smaller entrant radii.

    entrantAngle = 15;

elseif nozzleType == 2 % Rao Approximation - 80% Bell Nozzle
    %throatDSFilet = 0.382*throatRadius;
    throatDSFilet = 1.5*throatRadius;
    % Typical value used in industry for any application other than conical
    % nozzles.
    
    entrantAngle = 22;
    % nu_e/2 is typically the value used for minimum entrant angle. The
    % 2/3rds rule is applied here as a "factor of safety" of nozzle
    % expansion to guarantee slow-enough expansion to produce minimal
    % additional inefficiencies.
elseif nozzleType == 3
    fprintf('Perfect Nozzle selected, ensure correct geometry input before use...\n')

    throatDSFilet = 1*throatRadius;
    [~,PMAngle,~] = flowprandtlmeyer(specificHeatRatio,exitMach);
    entrantAngle = PMAngle/2;
    % nu_e/2 is allowed to be used in this case, since Method of
    % Characteristics is a method which employs analyzing the expansion
    % characteristics of the supersonic flow to produce no shocks.

end
%% INTERIOR GEOMETRY

% Generates the downstream throat radius
% Curve from -90 [deg] -> throat exit angle
ds_throat_angles=linspace(-90,-90+entrantAngle,contourPoints);
center_x_ds = 0;
center_y_ds = throatRadius + throatDSFilet;
x_ds = throatDSFilet*cosd(ds_throat_angles)+center_x_ds;
x_ds = x_ds(1:end-1);
y_ds = throatDSFilet*sind(ds_throat_angles)+center_y_ds;
y_ds = y_ds(1:end-1);
% Generates the upstream throat radius
% Curve from throat entrant angle -> -90 [deg]
throatUSFilet = 1.5*throatRadius;
us_throat_angles=linspace(270-interiorContractionAngle,270,contourPoints);
center_x_us = x_ds(1) - throatUSFilet*cosd(270);
center_y_us = y_ds(1) - throatUSFilet*sind(270);
x_us = throatUSFilet*cosd(us_throat_angles) + center_x_us;
x_us = x_us(1:end-1);
y_us = throatUSFilet*sind(us_throat_angles) + center_y_us;
y_us = y_us(1:end-1);

% Generates the combustion chamber filet
y1 = chamberRadius-y_us(1);
x1 = y1/tand(interiorContractionAngle/2);
y2 = x1*tand(90-interiorContractionAngle);

flt_ch = y1+y2;

ch_start_y = chamberRadius - flt_ch*(1 - sind(interiorContractionAngle));
ch_start_x = x_us(1) - (ch_start_y - y_us(1));
ch_angles = linspace(90,90-interiorContractionAngle,contourPoints);
ch_center_y = chamberRadius - flt_ch;
ch_center_x = ch_start_x - flt_ch*cosd(-interiorContractionAngle);
ch_x = flt_ch*cosd(ch_angles) + ch_center_x;
ch_x = ch_x(1:end-1);

ch_y = flt_ch*sind(ch_angles) + ch_center_y;
ch_y = ch_y(1:end-1);

% Generates the length of the combustion chamber.
ch_len_x = linspace(-chamberLength,ch_x(1),contourPoints);
ch_len_x = ch_len_x(1:end-1);
ch_len_y = linspace(chamberRadius,ch_y(1),contourPoints);
ch_len_y = ch_len_y(1:end-1);

%% NOZZLE GENERATION
% ALL MEASUREMENTS IN METERS
if nozzleType == 1 % 15-Degree Half Angle (Conical) Nozzle
    exitAngle=15;
    noz_dx=(exitRadius-throatRadius)/tand(exitAngle);
    
    noz_x=linspace(x_ds(end),x_ds(end)+noz_dx,contourPoints); 
    noz_x=noz_x(2:end);
    noz_y=linspace(y_ds(end),exitRadius,contourPoints);
    noz_y=noz_y(2:end);

    interiorContourX = [ch_len_x,ch_x,x_us,x_ds,noz_x];
    interiorContourY = [ch_len_y,ch_y,y_us,y_ds,noz_y];
elseif nozzleType == 2 % Rao Nozzle
% http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf


    bellPercent = 0.8;
    throatArea=pi*throatRadius^2;
    exitArea=pi*exitRadius^2;
    expAR=exitArea/throatArea;
    exitAngle=13.5;

    % Quadratic Bezier Curve
    N_x = x_ds(end);
    N_y = y_ds(end);
    E_x = bellPercent*((sqrt(expAR) - 1)*throatRadius/tand(15)); 
    E_y = exitRadius;
    m1 = tand(entrantAngle); 
    m2 = tand(exitAngle); 
    C1 = N_y - m1*N_x; 
    C2 = E_y - m2*E_x; 
    Q_x = (C2 - C1)/(m1 - m2); 
    Q_y = (m1*C2 - m2*C1)/(m1 - m2); 
    t=linspace(0,1,contourPoints);
    % Nozzle Curve Generation
    noz_x = (1 - t).^2*N_x + 2*(1 - t).*t*Q_x + t.^2*E_x;
    noz_x = noz_x(2:end);
    noz_y = (1 - t).^2*N_y + 2*(1 - t).*t*Q_y + t.^2*E_y;
    noz_y = noz_y(2:end);

    interiorContourX = [ch_len_x,ch_x,x_us,x_ds,noz_x];
    interiorContourY = [ch_len_y,ch_y,y_us,y_ds,noz_y];

elseif nozzleType == 3

%[noz_x,noz_y] = methodOfCharacteristics(Rt,gam,theta_n);

% Currently requires dependance on NASA's Method of Characteristics script
% to produce nozzle geometry.
load('input/wall_contour_1Rth.mat');
noz_x = wall_contour_1Rth(:,1)'.*throatRadius; noz_y = wall_contour_1Rth(:,2)'.*throatRadius;

interiorContourX = [ch_len_x,ch_x,x_us,noz_x];
interiorContourY = [ch_len_y,ch_y,y_us,noz_y];

end
%% OUTER CONTOUR GENERATION
% Outer Chamber Straight
outerChamberX = linspace(-chamberLength,ch_x(1),contourPoints);
outerChamberY = linspace(chamberRadius+chamberThickness,chamberRadius+chamberThickness,contourPoints);

% outer chamber filet


flt_ch_out = flt_ch+chamberThickness;
ch_angles = linspace(90,90-exteriorContractionAngle,contourPoints);
ch_center_y = chamberRadius+chamberThickness - flt_ch_out;
ch_center_x = outerChamberX(end);
ch_flt_x = flt_ch_out*cosd(ch_angles) + ch_center_x; ch_flt_x = ch_flt_x(2:end);
ch_flt_y = flt_ch_out*sind(ch_angles) + ch_center_y; ch_flt_y = ch_flt_y(2:end);

% nozzle offset
[noz_x_out,noz_y_out] = offsetCurve(noz_x,noz_y,nozzleThickness);
noz_fit = fit(noz_x_out',noz_y_out','pchipinterp');

noz_x_out = linspace(noz_x_out(1),noz_x_out(end),contourPoints);
noz_y_out = noz_fit(noz_x_out)';

% outer chamber contraction angle

str_start_y = ch_flt_y(end); str_start_x = ch_flt_x(end);
str_end_y = throatRadius;
str_end_x = str_start_x + (str_start_y - str_end_y)/tand(exteriorContractionAngle);

str_x = linspace(str_start_x, str_end_x, contourPoints);
str_y = linspace(str_start_y, str_end_y, contourPoints);

% nozzle exterior filet

for i=1:length(noz_x_out)-1

    theta_noz_start = atand((noz_y_out(i+1)-noz_y_out(i))/(noz_x_out(i+1)-noz_x_out(i)));
    
    noz_flt_angles = linspace(-90+theta_noz_start,-90-exteriorContractionAngle,contourPoints);
    
    noz_flt_start_x = noz_x_out(i);
    noz_flt_start_y = noz_y_out(i);
    
    noz_flt_x = nozzleFilet*cosd(noz_flt_angles);
    noz_flt_y = nozzleFilet*sind(noz_flt_angles);
    
    diff_x = abs(noz_flt_x(1) - noz_x_out(i));
    diff_y = abs(noz_flt_y(1) - noz_y_out(i));
    
    noz_flt_x = noz_flt_x + diff_x;
    noz_flt_y = noz_flt_y + diff_y;

    closest = dsearchn([str_x;str_y]',[noz_flt_x(end);noz_flt_y(end)]');
    dist(i) = sqrt((noz_flt_x(end)-str_x(closest))^2+(noz_flt_y(end)-str_y(closest))^2);

end

% locating closest approach
[~,mind] = min(dist);

    theta_noz_start = atand((noz_y_out(mind+1)-noz_y_out(mind))/(noz_x_out(mind+1)-noz_x_out(mind)));
    
    noz_flt_angles = linspace(-90+theta_noz_start,-90-exteriorContractionAngle,contourPoints);
    
    noz_flt_x = nozzleFilet*cosd(noz_flt_angles);
    noz_flt_y = nozzleFilet*sind(noz_flt_angles);
    
    diff_x = abs(noz_flt_x(1) - noz_x_out(mind));
    diff_y = abs(noz_flt_y(1) - noz_y_out(mind));
    
    noz_flt_x = noz_flt_x + diff_x; noz_flt_x = noz_flt_x(2:end);
    noz_flt_y = noz_flt_y + diff_y; noz_flt_y = noz_flt_y(2:end);

% removing unused part of noz geometry
noz_x_out = noz_x_out(mind:end); noz_y_out = noz_y_out(mind:end);

% fixing the chamber straight
str_end_x = noz_flt_x(end);
str_end_y = noz_flt_y(end);

str_x = linspace(str_start_x, str_end_x, contourPoints); str_x = str_x(2:end-1);
str_y = linspace(str_start_y, str_end_y, contourPoints); str_y = str_y(2:end-1);

% extend wall flat to nozzle end
theta_noz_end = atand((noz_y_out(end)-noz_y_out(end-1))/(noz_x_out(end)-noz_x_out(end-1)));

x_end = interiorContourX(end);
y_end = noz_y_out(end) + (x_end-noz_x_out(end))*tand(theta_noz_end);

noz_end_x = linspace(noz_x_out(end),x_end,contourPoints); noz_end_x = noz_end_x(2:end);
noz_end_y = linspace(noz_y_out(end),y_end,contourPoints); noz_end_y = noz_end_y(2:end);

exteriorContourX = [outerChamberX, ch_flt_x, str_x, noz_flt_x, noz_x_out, noz_end_x];
exteriorContourY = [outerChamberY, ch_flt_y, str_y, noz_flt_y, noz_y_out, noz_end_y];

%% SMOOTHING SPLINE
% Relocating chamber end to be at x=0
interiorContourX = interiorContourX - interiorContourX(1);
exteriorContourX = exteriorContourX - exteriorContourX(1);

% Linear curve fit to evenly separate the x-domain datapoints.
[interiorFit,~,~,~] = fit(interiorContourX',interiorContourY','pchipinterp');
x_new_in = linspace(0,interiorContourX(end),contourPoints);
y_new_in = interiorFit(x_new_in)';

% Linear curve fit to evenly separate the x-domain datapoints.
[exteriorFit,~,~,~] = fit(exteriorContourX',exteriorContourY','pchipinterp');
x_new_out = linspace(0,exteriorContourX(end),contourPoints);
y_new_out = exteriorFit(x_new_in)';
%% COMBINED MATRICES INTO FINAL CURVE
interiorEngineContour = [x_new_in;y_new_in];
exteriorEngineContour = [x_new_out;y_new_out];

%%  Export Full Contour Data to '.csv' in [cm]:
% Interior
    z_data=zeros(1,size(interiorEngineContour,2));
    output_matrix=[interiorEngineContour(1,:);interiorEngineContour(2,:);z_data]';
    
    outputFolder = 'output';
    fileName = 'interiorEngineContour.txt';
    filePath = fullfile(outputFolder,fileName);
    writematrix(output_matrix,filePath)
    
    % Outputs file path where the 'contour.csv' file was generated
    fprintf('''%s'' generated at %s\n',fileName,filePath);
% Exterior
    z_data=zeros(1,size(exteriorEngineContour,2));
    output_matrix=[exteriorEngineContour(1,:);exteriorEngineContour(2,:);z_data]';
    
    outputFolder = 'output';
    fileName = 'exteriorEngineContour.txt';
    filePath = fullfile(outputFolder,fileName);
    writematrix(output_matrix,filePath)
    
    % Outputs file path where the 'contour.csv' file was generated
    fprintf('''%s'' generated at %s\n',fileName,filePath);
end

