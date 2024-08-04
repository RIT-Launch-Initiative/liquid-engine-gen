%% TITLE: Heat Transfer FEA Solver
% Description: Utilizes FEA concepts to do a 2D transient heat transfer
% analysis of a liquid rocket engine geometry.
% Developer: Evan Olsen

% Special thanks to Chad Wake for helping me understand the theory and
% provide his code to reference off of.

% Assumptions -> 
% 1) All nodes are relatively evenly spaced
% 2) Adiabatic at the beginning of the chamber and end of the nozzle
% 3) Assuming convection coefficient to air of ~28 W/m^2-K
% 4) Constant material properties

% Future possibilities: Determine more accurate wall-to-air convection
% coefficient using Nusselt correlations, account for non-linear material
% properties with regards to temperature, calculate P and Q laplacians to
% account for the difference in mesh spacing.

function [T,Z,R,burnTime,dt,flagFailedAnalysis] = transientHeatTransfer(Z,R,n_z,n_r,wallThickness,h_g,T_g,burnTime,atmosphericTemperature,thermalDiffusivity,thermalConductivity)
%% 2D TRANSIENT HEAT EQUATION
progressbar('transientHeatTransfer()...')
flagFailedAnalysis = 0;
% TIME_END -> FUNCTION INPUT
% THERMAL_DIFFUSIVITY -> FUNCTION INPUT
% THERMAL_CONDUCTIVITY -> FUNCTION INPUT

% Differentials & Timesteps
dr = mean(mean(R(2:end,:)-R(1:end-1,:)));
dz = mean(mean(Z(:,2:end)-Z(:,1:end-1)));

dt = (1/(2*thermalDiffusivity*(1/dr^2+1/dz^2)))/5;
totalTimesteps = round(burnTime/dt);

itercount = 0;
% Initial & Boundary Conditions
% Inner boundary condition
h_o = h_g; % [W/m^2-K]
if(length(T_g(1,:)) == 1)
    T_g = T_g';
end
T_inf_o = T_g; % [K]

% Outer boundary condition
h_L = linspace(28.372,28.372,n_z); % [W/m^2-K]
% T_inf_L -> FUNCTION INPUT
% Initial temperature condition
T = atmosphericTemperature.*ones(n_r,n_z,totalTimesteps+1);

% Initializing matrices
T_f_o = zeros(1,n_z);
T_f_L = zeros(1,n_z);
Tzz = zeros(n_r,n_z);
Tzn = zeros(n_r,n_z);
Tnn = zeros(n_r,n_z);
%Tz = zeros(n_r,n_z);
%Tn = zeros(n_r,n_z);
laplacianT = zeros(n_r,n_z);

% Initialize eta and zeta x and y-axis transformations
% Only the center ranges when performing centered difference approximation
% to PDEs.
% Matrix Indices
eta_center = 2:n_r-1;
zeta_center = 2:n_z-1;
eta_all = 1:n_r;
zeta_all = 1:n_z;
% Zeta & Eta Transformantion Matrices
% Zeta = repmat(zeta_all,n_r,1);
% Eta = repmat(eta_all,n_z,1)';
% Finite difference approximation for partial derivatives, in the inner
% region of the XY matrices. dx/deta, dy/deta, dx/dzeta, dy/dzeta; deta and
% dzeta will always be 1 if defined in integer values, therefore they can
% be neglected.

% f'(x) = {-3f(x) + 4f(x+dx) - f(x+2*dx)}/(2*dx): (forward difference approximation: dx^2)
% f'(x) = {f(x+dx) - f(x-dx)}/(2*dx): (centered difference approximation: dx^2)
% f'(x) = {3f(x) - 4f(x-dx) + f(x-2*dx)}/(2*dx): (backward difference approximation: dx^2)

% alpha = (dR/deta)^2 + (dz/deta)^2                             (11.45c)
% beta = (dz/dzeta)(dR/deta) + (dR/dzeta)(dz/deta)              (11.45d)
% gamma = (dz/dzeta)^2 + (dR/dzeta)^2                           (11.45e)
% jacobian = (dz/dzeta)(dR/deta) - (dz/deta)(dR/dzeta)          (11.45f)
% laplacianZeta = d^2zeta/dz^2 + d^2zeta/dR^2 + (1/R)*dzeta/dR  (11.45g)
% laplacianEta = d^2eta/dz^2 + d^2eta/dR^2 + (1/R)*deta/dR      (11.45h)

[alpha] = alphaCalc(Z,R,zeta_center,zeta_all,eta_center,eta_all,n_r,n_z);
[beta] = betaCalc(Z,R,zeta_center,zeta_all,eta_center,eta_all,n_r,n_z);
[gamma] = gammaCalc(Z,R,zeta_center,zeta_all,eta_center,eta_all,n_r,n_z);
[J] = jacobianCalc(Z,R,zeta_center,zeta_all,eta_center,eta_all,n_r,n_z);
for t = 1:totalTimesteps
    % Ficticious nodes i = -1, i = M+1
    % Converting convection boundary conditions into equivalent temperature
    % nodes outside of the computational domain
        
        % Ficticious node at the interior, reacts to combustion
        % temperature.
        T_f_o(zeta_center) = T(2,zeta_center,t) - 2./(gamma(1,zeta_center)+eps).*(beta(1,zeta_center).*(T(1,zeta_center+1,t)-T(1,zeta_center-1,t))./2 + (h_o(zeta_center).*(T(1,zeta_center,t)-T_inf_o(zeta_center))).*J(1,zeta_center).*gamma(1,zeta_center).^(0.5)./thermalConductivity);
        
        % Ficticious node at the exterior, reacts to atmospheric
        % temperature.
        T_f_L(zeta_center) = T(end-1,zeta_center,t) - 2./(gamma(end,zeta_center)+eps).*(beta(end,zeta_center).*(T(end,zeta_center+1,t)-T(end,zeta_center-1,t))./2 + (h_L(zeta_center).*(T(end,zeta_center,t)-atmosphericTemperature)).*J(end,zeta_center).*gamma(end,zeta_center).^(0.5)./thermalConductivity);
        
    % (d^2/dzeta^2)T
        % f''(x) = {f(x+dx) - 2f(x) + f(x-dx)}/dx^2: Second-order centered diff.
        Tzz(eta_all,zeta_center) = T(eta_all,zeta_center+1,t)-...
                                   2*T(eta_all,zeta_center,t)+...
                                   T(eta_all,zeta_center-1,t);
        % Left and right side boundaries: adiabatic in the z-direction
        Tzz(eta_all,1) = 0;
        Tzz(eta_all,end) = 0;
    
    % (d/dzeta)(d/deta)T
        % f'(x) = {f(x+dx) - f(x-dx)}/(2*dx): First-order centered diff.
        Tzn(eta_center,zeta_center) = (1/2).*((T(eta_center+1,zeta_center,t)-T(eta_center-1,zeta_center,t)).*... % dT/dzeta
                                      (1/2).*(T(eta_center,zeta_center+1,t)-T(eta_center,zeta_center-1,t)));     % dT/deta
        % Left and right side boundaries: adiabatic in the z-direction
        Tzn(eta_all,1) = 0;
        Tzn(eta_all,end) = 0;
        % Top and bottom boundaries, defined by the ficticious nodes which
        % transform our boundary condition into an equivalent temperature
        % node.
        Tzn(1,zeta_center) = (1/2).*(T(2,zeta_center,t)-T_f_o(zeta_center)).*...
                             (1/2).*(T(1,zeta_center+1,t)-T(1,zeta_center-1,t));
        Tzn(end,zeta_center) = (1/2).*(T_f_L(zeta_center)-T(end-1,zeta_center,t)).*...
                               (1/2).*(T(end,zeta_center+1,t)-T(end,zeta_center-1,t));
    
    % (d^2/deta^2)T
        % f''(x) = {f(x+dx) - 2f(x) + f(x-dx)}/dx^2: Second-order centered
        % diff.
        Tnn(eta_center,zeta_all) = T(eta_center+1,zeta_all,t)-...
                                   2*T(eta_center,zeta_all,t)+...
                                   T(eta_center-1,zeta_all,t);
        % Top and bottom boundaries, defined by the ficticious nodes which
        % transform our boundary condition into an equivalent temperature
        % due to conduction.
        Tnn(1,zeta_all) = T(2,zeta_all,t) - 2*T(1,zeta_all,t) + T_f_o;
        Tnn(end,zeta_all) = T_f_L - 2*T(end,zeta_all,t) + T(end-1,zeta_all,t);
    
    % Temperature Laplacian
    laplacianT(:,zeta_center) = (1./J(:,zeta_center).^2).*(alpha(:,zeta_center).*...
    Tzz(:,zeta_center)-2.*beta(:,zeta_center).*Tzn(:,zeta_center)+...
    gamma(:,zeta_center).*Tnn(:,zeta_center));
    
    T(:,:,t+1) = T(:,:,t)+thermalDiffusivity.*dt.*laplacianT; %Find the next T
    
    if(sum(any(T(:,:,t+1)<0)) > 0 || sum(sum(isnan(T(:,:,t+1))) > 0))
        fprintf('\n')
        warning('Temperature Matrix Contains Negative/NaN values. Please double check mesh generation results.');
        
        progressbar(1);
        flagFailedAnalysis = 1;
        Z = flipud(Z);
        R = flipud(R);
        T = flipud(T);
        T_end = T(:,:,t+1);
        break
    end
    itercount = itercount + 1;
    progressbar(itercount/totalTimesteps);
end
if(flagFailedAnalysis == 0)
    Z = flipud(Z);
    R = flipud(R);
    T = flipud(T);
    T_end = T(:,:,t+1);
end

end

%% FUNCTIONS
function [alpha] = alphaCalc(Z,R,zeta_center,zeta_all,eta_center,~,n_r,n_z)
% ALPHA CALCULATION
alpha = zeros(n_r,n_z);
% alpha -- top row
alpha(1,zeta_all) = ((-3.*Z(1,zeta_all)+4.*Z(2,zeta_all)-Z(3,zeta_all))./2).^2+... % (dx/deta)^2
                  ((-3.*R(1,zeta_all)+4.*R(2,zeta_all)-R(3,zeta_all))./2).^2;      % (dy/deta)^2
% alpha -- middle rows
alpha(eta_center,zeta_center) = ((Z(eta_center+1,zeta_center)-Z(eta_center-1,zeta_center))./2).^2+... % (dx/deta)^2
                  ((R(eta_center+1,zeta_center)-R(eta_center-1,zeta_center))./2).^2;                  % (dy/deta)^2
% alpha -- bottom row
alpha(end,zeta_all) = ((3.*Z(end,zeta_all)-4.*Z(end-1,zeta_all)+Z(end-2,zeta_all))./2).^2+... % (dx/deta)^2
                    ((3.*R(end,zeta_all)-4.*R(end-1,zeta_all)+R(end-2,zeta_all))./2).^2;      % (dy/deta)^2
end
function [beta] = betaCalc(Z,R,zeta_center,~,eta_center,~,n_r,n_z)
% BETA CALCULATION
beta = zeros(n_r,n_z);
% beta -- top row
beta(1,zeta_center) = (R(1,zeta_center+1)-R(1,zeta_center-1)).*...                         % dy/dzeta - centered diff.
                      (-3.*R(1,zeta_center)+4.*R(2,zeta_center)-R(3,zeta_center))./4 + ... % dy/deta - forward diff.
                      (Z(1,zeta_center+1)-Z(1,zeta_center-1)).*...                         % dx/dzeta - centered diff.
                      (-3.*Z(1,zeta_center)+4.*Z(2,zeta_center)-Z(3,zeta_center))./4;      % dx/deta - forward diff.

% beta -- middle rows
beta(eta_center,zeta_center) = ((Z(eta_center+1,zeta_center)-Z(eta_center-1,zeta_center))./2).*... % dx/deta
                 ((Z(eta_center,zeta_center+1)-Z(eta_center,zeta_center-1))./2)+...                % dx/dzeta
                 ((R(eta_center+1,zeta_center)-R(eta_center-1,zeta_center))./2).*...               % dx/deta
                 ((R(eta_center,zeta_center+1)-R(eta_center,zeta_center-1))./2);                   % dx/dzeta
% beta -- bottom row
beta(end,zeta_center) = ((R(end,zeta_center+1)-R(end,zeta_center-1))./2).*...                          % dy/dzeta - centered diff.
                        (3.*R(end,zeta_center)-4.*R(end-1,zeta_center)+R(end-2,zeta_center))./2 + ...  % dy/deta - backward diff.
                        ((Z(end,zeta_center+1)-Z(end,zeta_center-1))./2).*...                          % dx/dzeta - centered diff.
                        (3.*Z(end,zeta_center)-4.*Z(end-1,zeta_center)+Z(end-2,zeta_center))./2;       % dx/deta - backward diff.
end
function [gamma] = gammaCalc(Z,R,zeta_center,~,eta_center,~,n_r,n_z)
% GAMMA CALCULATION
gamma = zeros(n_r,n_z);
% gamma -- top row
gamma(1,zeta_center) = (((Z(1,zeta_center+1)-Z(1,zeta_center-1))./2).^2+... % (dx/dzeta)^2
                ((R(1,zeta_center+1)-R(1,zeta_center-1))./2).^2);           % (dx/dzeta)^2

% gamma -- middle rows
gamma(eta_center,zeta_center) = ((Z(eta_center,zeta_center+1)-Z(eta_center,zeta_center-1))./2).^2+... % (dx/dzeta)^2
                  ((R(eta_center,zeta_center+1)-R(eta_center,zeta_center-1))./2).^2;                  % (dy/dzeta)^2
% gamma -- bottom row
gamma(end,zeta_center) = (((Z(end,zeta_center+1)-Z(end,zeta_center-1))./2).^2+... % (dx/dzeta)^2
                  ((R(end,zeta_center+1)-R(end,zeta_center-1))./2).^2);           % (dy/dzeta)^2
end
function [J] = jacobianCalc(Z,R,zeta_center,~,eta_center,~,n_r,n_z)
% JACOBIAN CALCULATION
J = zeros(n_r,n_z);
% Jacobian -- top row
J(1,zeta_center) = (Z(1,zeta_center+1)-Z(1,zeta_center-1))./2.*...             % dx/dzeta - centered diff.
            (-3.*R(1,zeta_center)+4.*R(2,zeta_center)-R(3,zeta_center))./2-... % dy/deta - forward diff.
            (R(1,zeta_center+1)-R(1,zeta_center-1))./2.*...                    % dy/dzeta - centered diff.
            (-3.*Z(1,zeta_center)+4.*Z(2,zeta_center)-Z(3,zeta_center))./2;    % dx/deta - forward diff.

% Jacobian -- middle rows
J(eta_center,zeta_center) = ((Z(eta_center,zeta_center+1)-Z(eta_center,zeta_center-1))./2).*... % dx/dzeta
              ((R(eta_center+1,zeta_center)-R(eta_center-1,zeta_center))./2)-...                % dx/deta
              ((R(eta_center,zeta_center+1)-R(eta_center,zeta_center-1))./2).*...               % dy/dzeta
              ((Z(eta_center+1,zeta_center)-Z(eta_center-1,zeta_center))./2);                   % dy/deta

% Jacobian -- bottom row
J(end,zeta_center) = (Z(end,zeta_center+1)-Z(end,zeta_center-1))./2.*...                  % dx/dzeta - centered diff.
              (3.*R(end,zeta_center)-4.*R(end-1,zeta_center)+R(end-2,zeta_center))./2-... % dy/deta - backward diff.
              (R(end,zeta_center+1)-R(end,zeta_center-1))./2.*...                         % dy/dzeta - centered diff.
              (3.*Z(end,zeta_center)-4.*Z(end-1,zeta_center)+Z(end-2,zeta_center))./2;    % dx/deta - backward diff.

% Jacobian -- left row
J(eta_center,1) = (-3.*Z(eta_center,1)+4.*Z(eta_center,2)-Z(eta_center,3))./2.*... % dx/dzeta - forward diff.
           (R(eta_center+1,1)-R(eta_center-1,1))./2-...                            % dy/deta - centered diff.
           (-3.*R(eta_center,1)+4.*R(eta_center,2)-R(eta_center,3))./2.*...        % dy/dzeta - forward diff.
           (Z(eta_center+1,1)-Z(eta_center-1,1))./2;                               % dx/deta - centered diff.

% Jacobian -- right row
J(eta_center,end) = (3.*Z(eta_center,end)-4.*Z(eta_center,end-1)+Z(eta_center,end-2))./2.*... % dx/dzeta - backward diff.
             (R(eta_center+1,end)-R(eta_center-1,end))./2-...                                 % dy/deta - centered diff.
             (3.*R(eta_center,end)-4.*R(eta_center,end-1)+R(eta_center,end-2))./2.*...        % dy/dzeta - backward diff.
             (Z(eta_center+1,end)-Z(eta_center-1,end))./2;                                    % dx/deta - centered diff.

% Jacobian -- top right corner
J(1,end) = (3.*Z(1,end)-4.*Z(1,end-1)+Z(1,end-2))./2.*... % dx/dzeta - backward diff.
           (-3.*R(1,end)+4.*R(2,end)-R(3,end))./2-...     % dy/deta - forward diff.
           (3.*R(1,end)-4.*R(1,end-1)+R(1,end-2))./2.*... % dy/dzeta - backward diff.
           (-3.*Z(1,end)+4.*Z(2,end)-Z(3,end))./2;        % dx/deta - forward diff.

% Jacobian -- top left corner
J(1,1) = (-3.*Z(1,1)+4.*Z(1,2)-Z(1,3))./2.*... % dx/dzeta
         (-3.*R(1,1)+4.*R(2,1)-R(3,1))./2-...  % dy/deta
         (-3.*R(1,1)+4.*R(1,2)-R(1,3))./2.*... % dy/dzeta
         (-3.*Z(1,1)+4.*Z(2,1)-Z(3,1))./2;     % dx/deta

% Jacobian -- bottom left corner
J(end,1) = (-3.*Z(end,1)+4.*Z(end,2)-Z(end,3))./2.*...  % dx/dzeta - forward diff.
           (3.*R(end,1)-4.*R(end-1,1)+R(end-2,1))./2-...% dy/deta - backward diff.
           (-3.*R(end,1)+4.*R(end,2)-R(end,3))./2.*...  % dy/dzeta - forward diff.
           (3.*Z(end,1)-4.*Z(end-1,1)+Z(end-2,1))./2;   % dx/deta - backward diff.

% Jacobian -- bottom right corner
J(end,end) = (3.*Z(end,end)-4.*Z(end,end-1)+Z(end,end-2))./2.*... % dx/dzeta
             (3.*R(end,end)-4.*R(end-1,end)+R(end-2,end))./2-...  % dy/deta
             (3.*R(end,end)-4.*R(end,end-1)+R(end,end-2))./2.*... % dy/dzeta
             (3.*Z(end,end)-4.*Z(end-1,end)+Z(end-2,end))./2;     % dx/deta
end