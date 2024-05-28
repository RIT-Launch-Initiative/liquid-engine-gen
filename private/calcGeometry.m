%% Function Declaration
function [engineContour] = calcGeometry(Rc,Re,Rt,Lc,nozzle_type,points,Me,gam)
%% CHAMBER & THROAT CURVES GENERATION
% Initial Values
if nozzle_type == 1 % Minimum Length Conical Nozzle
    flt_ex = 1.5*Rt;
    % 1.5*Rt is used exclusively for conical nozzles in common application.
    % Unlike other nozzle types which have much smaller entrant radii.

    [~,nu_e,~] = flowprandtlmeyer(gam,Me);
    theta_n = nu_e/3; 
    % nu_e/2 is typically the value used for minimum entrant angle. The
    % 2/3rds rule is applied here as a "factor of safety" of nozzle
    % expansion to guarantee slow-enough expansion to produce minimal
    % additional inefficiencies.
elseif nozzle_type == 2 % Rao Approximation - 80% Bell Nozzle
    flt_ex = 1*Rt;
    % Typical value used in industry for any application other than conical
    % nozzles.
    
    [~,nu_e,~] = flowprandtlmeyer(gam,Me);
    theta_n = nu_e/3;
    theta_n = 23;
    % nu_e/2 is typically the value used for minimum entrant angle. The
    % 2/3rds rule is applied here as a "factor of safety" of nozzle
    % expansion to guarantee slow-enough expansion to produce minimal
    % additional inefficiencies.
elseif nozzle_type == 3
    flt_ex = 1.5*Rt;
    [~,nu_e,~] = flowprandtlmeyer(gam,Me);
    theta_n = nu_e/2;
    % nu_e/2 is allowed to be used in this case, since Method of
    % Characteristics is a method which employs analyzing the expansion
    % characteristics of the supersonic flow to produce no shocks.
end
% Less critical input parameters, not in the main input script due to being
% less impactful on performance, and typical values are inputted here.
theta_ch=45; % deg
flt_ch = 1*Rt;
flt_ent = 1.5*Rt;

% Generates the downstream throat radius
% Curve from -90 [deg] -> throat exit angle
ds_throat_angles=linspace(-90,-90+theta_n,points);
center_x_ds = 0;
center_y_ds = Rt + flt_ex;
x_ds = flt_ex*cosd(ds_throat_angles)+center_x_ds;
x_ds = x_ds(1:end-1);
y_ds = flt_ex*sind(ds_throat_angles)+center_y_ds;
y_ds = y_ds(1:end-1);
% Generates the upstream throat radius
% Curve from throat entrant angle -> -90 [deg]
us_throat_angles=linspace(270-theta_ch,270,points);
center_x_us = x_ds(1) - flt_ent*cosd(270);
center_y_us = y_ds(1) - flt_ent*sind(270);
x_us = flt_ent*cosd(us_throat_angles) + center_x_us;
x_us = x_us(1:end-1);
y_us = flt_ent*sind(us_throat_angles) + center_y_us;
y_us = y_us(1:end-1);
% Generates the combustion chamber entrant radius
ch_start_y = Rc - flt_ch*(1 - sind(theta_ch));
ch_start_x = x_us(1) - (ch_start_y - y_us(1));
ch_angles = linspace(90,90-theta_ch);
ch_center_y = Rc - flt_ch;
ch_center_x = ch_start_x - flt_ch*cosd(-theta_ch);
ch_x = flt_ch*cosd(ch_angles) + ch_center_x;
ch_y = flt_ch*sind(ch_angles) + ch_center_y;
% Generates the constant-angle combustion chamber contraction curve.
ch_contr_x = linspace(ch_start_x,x_us(1),points);
ch_contr_x = ch_contr_x(2:end-1);
ch_contr_y = linspace(ch_start_y,y_us(1),points);
ch_contr_y = ch_contr_y(2:end-1);
% Generates the length of the combustion chamber.
ch_len_x = linspace(-Lc,ch_x(1),points);
ch_len_x = ch_len_x(1:end-1);
ch_len_y = linspace(Rc,ch_y(1),points);
ch_len_y = ch_len_y(1:end-1);
%% NOZZLE GENERATION
% ALL MEASUREMENTS IN METERS
if nozzle_type == 1 % 15-Degree Half Angle (Conical) Nozzle
    theta_e=15;
    noz_dx=(Re-Rt)/tand(theta_e);
    
    noz_x=linspace(x_ds(end),x_ds(end)+noz_dx,points); 
    noz_x=noz_x(2:end);
    noz_y=linspace(y_ds(end),Re,points);
    noz_y=noz_y(2:end);

    x_data = [ch_len_x,ch_x,ch_contr_x,x_us,x_ds,noz_x];
    y_data = [ch_len_y,ch_y,ch_contr_y,y_us,y_ds,noz_y];
elseif nozzle_type == 2 % Bell Nozzle
% http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf


    percBell = 0.8;
    A_star=pi*Rt^2;
    A_ex=pi*Re^2;
    AeAt=A_ex/A_star;
    theta_e=13;

    % Quadratic Bezier Curve
    N_x = x_ds(end);
    N_y = y_ds(end);
    E_x = percBell*((sqrt(AeAt) - 1)*Rt/tand(15)); 
    E_y = Re;
    m1 = tand(theta_n); 
    m2 = tand(theta_e); 
    C1 = N_y - m1*N_x; 
    C2 = E_y - m2*E_x; 
    Q_x = (C2 - C1)/(m1 - m2); 
    Q_y = (m1*C2 - m2*C1)/(m1 - m2); 
    t=linspace(0,1,points);
    % Nozzle Curve Generation
    noz_x = (1 - t).^2*N_x + 2*(1 - t).*t*Q_x + t.^2*E_x;
    noz_x = noz_x(2:end);
    noz_y = (1 - t).^2*N_y + 2*(1 - t).*t*Q_y + t.^2*E_y;
    noz_y = noz_y(2:end);

    x_data = [ch_len_x,ch_x,ch_contr_x,x_us,x_ds,noz_x];
    y_data = [ch_len_y,ch_y,ch_contr_y,y_us,y_ds,noz_y];

elseif nozzle_type == 3

%[noz_x,noz_y] = methodOfCharacteristics(Rt,gam,theta_n);

% Currently requires dependance on NASA's Method of Characteristics script
% to produce nozzle geometry.
load('input/wall_contour_1Rth.mat');
noz_x = wall_contour_1Rth(:,1)'.*Rt; noz_y = wall_contour_1Rth(:,2)'.*Rt;

x_data = [ch_len_x,ch_x,ch_contr_x,x_us,noz_x];
y_data = [ch_len_y,ch_y,ch_contr_y,y_us,noz_y];

end

x_data = x_data - x_data(1);

% Linear curve fit to evenly separate the x-domain datapoints.
[engine_fit,~,~,~] = fit(x_data',y_data','pchipinterp');
x_new = linspace(0,x_data(end),points);
y_new = engine_fit(x_new); y_new = y_new';
%% COMBINED MATRICES INTO FINAL CURVE
engineContour = [x_new;y_new];
%%  Export Full Contour Data to '.csv' in [cm]:
    z_data=zeros(1,size(engineContour,2));
    output_matrix=[engineContour(1,:);engineContour(2,:);z_data]';
    
    outputFolder = 'output';
    fileName = 'interiorEngineContour.txt';
    filePath = fullfile(outputFolder,fileName);
    writematrix(output_matrix,filePath)
    
    % Outputs file path where the 'contour.csv' file was generated
    fprintf('''%s'' generated at %s\n',fileName,filePath);
end

