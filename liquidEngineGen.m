%% TITLE: Liquid Bi-propellant Rocket Engine Generator
% Developer: Evan Olsen
% Date: 05-17-2024
% Description: This script is to be used to develop a small liquid
% bi-propellant rocket engine. It has been tested to be operational for the
% design of a 150 [lbf] liquid rocket engine using nitrous oxide and
% isopropanol as the propellants.

% startup() runs a function which sets default figure fonts, sizes, text
% interpreters, etc. for consistent output on all machines.
startup
% deleteOutput() deletes any files lingering in the output directory to be
% replaced with new outputs after running the program again.
deleteOutput(pwd)
%% ASSUMPTIONS
% Isentropic flow, constant specific heat, frozen reactions following the
% throat, equillibrium reactions prior to the throat.
%% INPUTS -- CONSTANTS
INPUT.accelDueToGrav = 9.80665;                                                         % [m/s^2] Acceleration Due to Gravity
INPUT.universalGasConst = 8314.46261815324;                                             % [J/mol-K] Universal Gas Constant
INPUT.atmosphericPressure = 14.7;                                                         % [psia] Atmospheric Pressure
INPUT.atmosphericTemperature = 298.15; % [K] Atmospheric Temperature
%% INPUTS -- GENERAL
INPUT.nominalThrustLbf = 151;                                                        % [lbf] Nominal Thrust

INPUT.characteristicLength = 0;                                                           % [m] Characteristic Length
% Set to '0' to use default aspect ratio of 2.0
% Select from www.braeunig.us/space/propulsion.htm#engine, pick the largest
% L* that is reasonable to account for injector/combustion inefficiencies.
% Use L* for engines of similar size class & propellant combonations.

INPUT.chamberPressure = 430;                                                            % [psia] Chamber Pressure
% Higher chamber pressure -> higher c* -> higher fuel efficiency
% For best performance, increase chamber pressure. Choose the highest
% chamber pressure your system can handle while being reasonably less than
% your feed/tank pressure.

% The following values will be used to iterate CEA. Choose the maximum
% reasonable contraction area ratio, and maximum and minimum allowable O/F
% ratios to be used in the iteration.
INPUT.maxContractionAR = 20;                                                       % [-] Contraction Area Ratio
INPUT.maxExpansionAR = 15;

INPUT.maxOF=20;                                                           % [-] O/F Ratio
INPUT.minOF = 0.1;

% FUEL PARAMETERS
% Used by NASA CEA to determine transport properties.
INPUT.fuelType='C3H8O';
INPUT.fuelEnthalpy=-272.8;                                                       % [kJ/mol] Standard Enthalpy of Formation
                                                                           % Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C67630&Mask=1#Thermo-Gas
INPUT.fuelPercWeight=100;                                                         % [%] Weight Percentage
INPUT.fuelTemperature=298.15;                                                       % [K] Temperature
INPUT.fuelDensity=0.785;                                                      % [g/cc] Density

% OXIDIZER PARAMETERS
% Used by NASA CEA to determine transport properties
INPUT.oxidizerType='N2O';
INPUT.oxidizerEnthalpy=82.05;                                                        % [kJ/mol] Standard Enthalpy of Formation
                                                                           % Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=1
INPUT.oxidizerPercWeight=100;                                                         % [%] Percentage Weight
INPUT.oxidizerTemperature=298.15;                                                       % [K] Temperature
INPUT.oxidizerDensity=1.220;                                                      % [g/cc] Density

INPUT.nozzleTypeInt = 3; 
% 1 = Minimum Length Conical Nozzle
% 2 = Minimum Length Bell Nozzle
% 3 = Method of Characteristics (WIP)

% The number of iterations controls the "resolution" of the matrices used
% by NASA CEA to determine thermophysical properties of the reaction, as
% well as performance values. A default of 100 is used for estimation of
% these properties. Set to a larger value (computation time will scale
% linearly) to determine more accurate values, once "ballpark" answers are
% validated.
INPUT.CEAIterations = 100;
% The number of points that will be used on the final output contour of the
% thrust chamber, as well as the resolution of the subsidiary output
% matrices.
INPUT.contourPoints = 1000;
%% INPUTS -- OVERRIDES
% Set these values equal to zero in order for them not to be applied.
INPUT.overrideOF = 2.092; % [-]
INPUT.overrideChamberRadius = 1; % [in]

%% INPUTS -- EFFICIENCIES
INPUT.characteristicVelocityEfficiency = 0.85; % Lowest value of typical c*-eff range, 92% - 99.5%
% Chamber adiabatic flame temperature correction: ~Cstar_eff^2

INPUT.nozzleEfficiency = 0; % (WIP)
% Set to '1' to consider nozzle inefficiencies
% Set to '0' to not consider nozzle inefficiencies
% Nozzle efficiency correction is dependent on the half-angle between the
% throat and the nozzle exit, and is calculated internally from nozzle
% geometry generation.

%% INPUTS -- BOOLEANS
INPUT.doTempCorrection = 1;
% Set to '1' to enable temperature correction through the c*-efficiency,
% Set to '0' to ignore c*-eff and use ideal combustion temperature values.
% This will effect the heat transfer analysis, if applicable.

INPUT.doFiniteElementAnalysis = 1;
% Set to '1' to enable finite element analysis
% Set to '0' to disable finite element analysis

INPUT.doExportGraphics = 0;
% Set to '1' to enable vector graphics export
% Set to '0' to disable vector graphics export
%% INPUTS -- FINITE ELEMENT ANALYSIS
% Mesh Sizing
INPUT.zDivisions = 1000;
INPUT.rDivisions = 20;

% Wall Thickness
INPUT.chamberThickness_in = 0.15; % [in]
INPUT.throatThickness_in = 0.125; % [in]
INPUT.nozzleThickness_in = 0.0625; % [in] 1/16th

% Transient Analysis Duration
INPUT.burnTime = 5; % [s]

% Material properties
INPUT.thermalDiffusivity = 4.2e-6; % [m^2/s]
INPUT.thermalConductivity = 16.2; % [W/m-K]
INPUT.CTE = 17.3 * 10^-6; % [1/K]
INPUT.poissonsRatio = 0.29; % [-]

%% NASA CEA -- O/F RATIO & CONTRACTION AREA RATIO OPTIMIZATION (MAX C*)
[OUTPUT.CEA,OUTPUT.contractionAR,OUTPUT.expansionAR,OUTPUT.OFRatio] = iterCEA(INPUT.overrideOF,INPUT.CEAIterations, ...
    INPUT.maxContractionAR,INPUT.minOF,INPUT.maxOF,INPUT.chamberPressure,INPUT.fuelType,INPUT.fuelPercWeight, ...
    INPUT.fuelEnthalpy,INPUT.fuelTemperature,INPUT.fuelDensity,INPUT.oxidizerType,INPUT.oxidizerPercWeight,INPUT.oxidizerEnthalpy, ...
    INPUT.oxidizerTemperature,INPUT.oxidizerDensity,INPUT.atmosphericPressure,INPUT.maxExpansionAR);

OUTPUT.dynamicViscosity        = OUTPUT.CEA.output.eql.viscosity;                      % [cP] Viscosity
OUTPUT.idealCharacteristicVelocity = OUTPUT.CEA.output.eql.cstar(3);                       % [m/s] Ideal Combustion Efficiency
OUTPUT.characteristicVelocity       = OUTPUT.idealCharacteristicVelocity*INPUT.characteristicVelocityEfficiency;                   % [m/s] Corrected Combustion Efficiency
OUTPUT.molarWeight          = OUTPUT.CEA.output.eql.mw(3);                          % [g/mol] Molecular Weight
OUTPUT.specificHeatRatio       = OUTPUT.CEA.output.eql.gamma(3);                       % [-] Specific Heat Ratio
OUTPUT.idealAdiabaticFlameTemperature     = OUTPUT.CEA.output.eql.temperature;                    % [K] Ideal Adiabatic Flame Temperature
OUTPUT.adiabaticFlameTemperature           = OUTPUT.idealAdiabaticFlameTemperature*(INPUT.characteristicVelocityEfficiency)^2;               % [K] Corrected Adiabatic Flame Temperature
OUTPUT.density        = OUTPUT.CEA.output.eql.density;                        % [kg/m^3] Density
OUTPUT.idealMachNumber          = OUTPUT.CEA.output.eql.mach;                           % [-] Mach Number
OUTPUT.sonicVelocity      = OUTPUT.CEA.output.eql.sonvel;                         % [m/s] Local Speed of Sound
OUTPUT.specificHeat          = OUTPUT.CEA.output.eql.cp_tran;                        % [J/kg-K] Specific Heat (Constant Pressure)
OUTPUT.prandtlNumber          = OUTPUT.CEA.output.eql.prandtl;                        % [-] Prandtl Number
OUTPUT.localPressureBar           = OUTPUT.CEA.output.eql.pressure;                       % [bar] Local Pressure

if(INPUT.doTempCorrection == 1)
    OUTPUT.analysisTemperature = OUTPUT.adiabaticFlameTemperature;
else
    OUTPUT.analysisTemperature = OUTPUT.idealAdiabaticFlameTemperature;
end
%% PERFORMANCE/CRITICAL GEOMETRY VALUES
% Calculating exit mach number using isentropic flow relations.
% Thermophysical properties are assumed frozen at the throat region, since
% the nozzle geometry depends on the isentropic flow relation for the
% equations to be valid. CEA uses "moving" equilibrium, recalculating the
% thermophysical properties at each region to be used to determine the
% next.
OUTPUT.exitMachNumber = flowisentropic(OUTPUT.specificHeatRatio,OUTPUT.expansionAR,'sup');                % [-] Exit Mach Number (Frozen)
% Convert input parameters into the SI unit system.
% Converting pressure outputted in [bar] to [Pa]
OUTPUT.exitPressure = OUTPUT.localPressureBar(4)*100000;                                            % [Pa] Exit Pressure
OUTPUT.chamberPressure = OUTPUT.localPressureBar(1)*100000;                                            % [Pa] Chamber Pressure
% Converting thrust inputted in [lbf] to [N]
OUTPUT.nominalThrust = convforce(INPUT.nominalThrustLbf,'lbf','N');                         % [N] Nominal Thrust

OUTPUT.idealExitVelocity = OUTPUT.idealMachNumber(4)*OUTPUT.sonicVelocity(4);                           % [m/s] Ideal Exit Velocity
OUTPUT.exitVelocity = OUTPUT.idealExitVelocity*INPUT.characteristicVelocityEfficiency;                               % [m/s] Exit Velocity
% Assuming the flow is perfectly expanded
OUTPUT.idealMassFlowrate = OUTPUT.nominalThrust/OUTPUT.idealExitVelocity;                         % [kg/s] Ideal Total Mass Flowrate
OUTPUT.massFlowrate = OUTPUT.nominalThrust/OUTPUT.exitVelocity;                                     % [kg/s] Total Mass Flowrate
OUTPUT.specificGasConstant = INPUT.universalGasConst/OUTPUT.molarWeight;                                         % [J/kg*K] Specific Gas Constant
OUTPUT.throatArea = calcAstar(OUTPUT.massFlowrate,OUTPUT.chamberPressure,OUTPUT.analysisTemperature(2),OUTPUT.specificGasConstant, ...
                      OUTPUT.specificHeatRatio);                                       % [m^2] Throat Area
OUTPUT.exitArea = OUTPUT.expansionAR*OUTPUT.throatArea;                                         % [m^2] Exit Area
OUTPUT.chamberArea = OUTPUT.contractionAR*OUTPUT.throatArea;                                         % [m^2] Chamber Area

OUTPUT.idealSpecificImpulse = OUTPUT.idealExitVelocity/INPUT.accelDueToGrav;                                % [s] Ideal Specific Impulse
OUTPUT.specificImpulse = OUTPUT.exitVelocity/INPUT.accelDueToGrav;                                            % [s] Specific Impulse
OUTPUT.fuelMassFlowrate = OUTPUT.massFlowrate/(1+OUTPUT.OFRatio);                                 % [kg/s] Fuel Mass Flowrate
OUTPUT.oxidizerMassFlowrate = OUTPUT.massFlowrate - OUTPUT.fuelMassFlowrate;                              % [kg/s] Oxidizer Mass Flowrate

% Check if overridden chamber diameter was inputted, then verify if it is
% greater than the nominal chamber diameter. Performance increases are
% minimal for larger area ratios, so performance changes are unneccessary.
[OUTPUT.chamberArea] = checkOverride(INPUT.overrideChamberRadius,OUTPUT.chamberArea);
OUTPUT.contractionAR = OUTPUT.chamberArea/OUTPUT.throatArea;

OUTPUT.chamberVolume = INPUT.characteristicLength*OUTPUT.throatArea;                                         % [m^3] Chamber Volume

OUTPUT.chamberRadius = sqrt(OUTPUT.chamberArea/pi);                                            % [m] Chamber Radius
OUTPUT.chamberDiameter = 2*OUTPUT.chamberRadius;                                                   % [m] Chamber Diamater
OUTPUT.chamberDiameterInch = convlength(OUTPUT.chamberDiameter,'m','in');                             % [in] Chamber Diameter

if(INPUT.characteristicLength == 0)
    OUTPUT.chamberLength = OUTPUT.chamberDiameter*2;                                               % [m] Chamber Length
    OUTPUT.chamberLengthInch = convlength(OUTPUT.chamberLength,'m','in');                         % [in] Chamber Length
else
    OUTPUT.chamberLength = OUTPUT.chamberVolume/OUTPUT.chamberArea;                                       % [m] Chamber Length
    OUTPUT.chamberLengthInch = convlength(OUTPUT.chamberLength,'m','in');                         % [in] Chamber Length
end

OUTPUT.throatRadius = sqrt(OUTPUT.throatArea/pi);                                            % [m] Throat Radius
OUTPUT.throatDiameter = 2*OUTPUT.throatRadius;                                                   % [m] Throat Diameter
OUTPUT.throatDiameterInch = convlength(OUTPUT.throatDiameter,'m','in');                             % [in] Throat Diameter

OUTPUT.exitRadius = sqrt(OUTPUT.exitArea/pi);                                            % [m] Exit Radius
OUTPUT.exitDiameter = 2*OUTPUT.exitRadius;                                                   % [m] Exit Diameter
OUTPUT.exitDiameterInch = convlength(OUTPUT.exitDiameter,'m','in');                             % [in] Exit Diamater
    
INPUT.chamberThickness = convlength(INPUT.chamberThickness_in,'in','m'); % [m]
INPUT.throatThickness = convlength(INPUT.throatThickness_in,'in','m'); % [m]
INPUT.nozzleThickness = convlength(INPUT.nozzleThickness_in,'in','m'); % [m]
%% THRUST CHAMBER GEOMETRY
OUTPUT.interiorEngineContour = calcGeometry(OUTPUT.chamberRadius,OUTPUT.exitRadius,OUTPUT.throatRadius, ...
    OUTPUT.chamberLength,INPUT.nozzleTypeInt,INPUT.contourPoints,OUTPUT.exitMachNumber,OUTPUT.specificHeatRatio);      % [m] Engine Geometry

%% AREA-MACH RELATIONSHIP
[OUTPUT.localAR,OUTPUT.localMachNumber,OUTPUT.localTemperature,OUTPUT.localPressure,OUTPUT.localDensity] = areaMach( ...
 OUTPUT.interiorEngineContour,OUTPUT.contractionAR,OUTPUT.throatArea,OUTPUT.specificHeatRatio, ...
 OUTPUT.analysisTemperature(2),OUTPUT.chamberPressure,OUTPUT.density(2));

%% FILM COEFFICIENT
[OUTPUT.FEA.h_g] = calcFilmCoefficient(OUTPUT.interiorEngineContour,OUTPUT.chamberRadius, ...
OUTPUT.massFlowrate,OUTPUT.dynamicViscosity,OUTPUT.specificHeat,OUTPUT.prandtlNumber);                              % [W/m-K] Film Coefficient

%% 2D TRANSIENT HEAT TRANSFER FEA
[OUTPUT.FEA.Z,OUTPUT.FEA.R,OUTPUT.FEA.Zeta,OUTPUT.FEA.Eta,OUTPUT.exteriorEngineContour,OUTPUT.wallThickness] = generateWallMesh(OUTPUT.interiorEngineContour,INPUT.chamberThickness,INPUT.throatThickness,INPUT.nozzleThickness,INPUT.zDivisions,INPUT.rDivisions);

if(INPUT.doFiniteElementAnalysis == 1)
    OUTPUT.FEA.interiorFlameTemperature = OUTPUT.localTemperature;
    
    % Run transient heat transfer solver
    [OUTPUT.FEA.transientTemperatureData,OUTPUT.FEA.Z,OUTPUT.FEA.R,OUTPUT.FEA.burnTime,OUTPUT.FEA.dt,OUTPUT.FEA.flagFailedAnalysis] = transientHeatTransfer(OUTPUT.FEA.Z,OUTPUT.FEA.R,INPUT.zDivisions,INPUT.rDivisions,OUTPUT.wallThickness,...
    OUTPUT.FEA.h_g,OUTPUT.FEA.interiorFlameTemperature,INPUT.burnTime,...
    INPUT.atmosphericTemperature,INPUT.thermalDiffusivity,INPUT.thermalConductivity);
    
%% 2D STRUCTURAL AND THERMAL STRESS ANALYSIS
    [OUTPUT.FEA.FOS,OUTPUT.FEA.vonMisesStress] = calcStress(OUTPUT.FEA.transientTemperatureData,OUTPUT.FEA.Z,OUTPUT.FEA.R,OUTPUT.localPressure,OUTPUT.interiorEngineContour,convpres(INPUT.atmosphericPressure,'psi','Pa'),INPUT.CTE,INPUT.poissonsRatio);
end
%% UI OUTPUT
OUTPUT.tgFigure = figure('Name','liquidEngineGen Output','Visible','off','NumberTitle', 'off','Position',[100 100 1000 800]);
OUTPUT.tg = uitabgroup(OUTPUT.tgFigure);
OUTPUT.t1 = uitab(OUTPUT.tg,'Title','Geometry');
OUTPUT.t2 = uitab(OUTPUT.tg,'Title','Fluid-Area Properties');
OUTPUT.t3 = uitab(OUTPUT.tg,'Title','Meshing');

plotGeometry(OUTPUT.interiorEngineContour,OUTPUT.exteriorEngineContour,OUTPUT.t1);
plotProperties(OUTPUT.interiorEngineContour,OUTPUT.localMachNumber,OUTPUT.localPressure,OUTPUT.localTemperature,OUTPUT.localDensity,OUTPUT.FEA.h_g,OUTPUT.t2);
plotMesh(OUTPUT.interiorEngineContour,OUTPUT.FEA.Z,OUTPUT.FEA.R,OUTPUT.FEA.Zeta,OUTPUT.FEA.Eta,OUTPUT.t3)

if(INPUT.doExportGraphics == 1)
    exportgraphics(OUTPUT.t1,'output\thrustChamberGeometry.pdf','ContentType','vector')
    exportgraphics(OUTPUT.t2,'output\fluidAreaProperties.pdf','ContentType','vector')
    exportgraphics(OUTPUT.t3,'output\meshOutput.pdf','ContentType','vector');
end

if(INPUT.doFiniteElementAnalysis == 1)
    if(OUTPUT.FEA.flagFailedAnalysis == 0)
        OUTPUT.t4 = uitab(OUTPUT.tg,'Title','Heat Diffusion Analysis');
        OUTPUT.t5 = uitab(OUTPUT.tg,'Title','Stress and F.O.S. Analysis');
    
        OUTPUT.tg.SelectedTab = OUTPUT.t4;
        plotHeatTransfer(OUTPUT.FEA.transientTemperatureData,OUTPUT.FEA.Z,OUTPUT.FEA.R,OUTPUT.FEA.burnTime,OUTPUT.t4);
        plotStress(OUTPUT.FEA.Z,OUTPUT.FEA.R,OUTPUT.FEA.FOS,OUTPUT.FEA.vonMisesStress,OUTPUT.FEA.burnTime,OUTPUT.t5)

        if(INPUT.doExportGraphics == 1)
            exportgraphics(OUTPUT.t4,'output\heatTransfer.pdf','ContentType','vector')
            exportgraphics(OUTPUT.t5,'output\vonMisesStress.pdf','ContentType','vector')
        end
    end
    OUTPUT.tg.SelectedTab = OUTPUT.t1;
end
OUTPUT.tgFigure.Visible = 'on';
%% COMMAND WINDOW OUTPUT
cprintf('*black','PERFORMANCE OVERVIEW\n')
cprintf('black','Ideal c*:      %.2f [m/s]\nCorrected c*:  %.2f [m/s]\n',OUTPUT.idealCharacteristicVelocity,OUTPUT.characteristicVelocity);
cprintf('black','Efficiency:    %.2f   [%%]\n\n',INPUT.characteristicVelocityEfficiency*100);

cprintf('black','Ideal Ve:      %.2f [m/s]\nCorrected Ve:  %.2f [m/s]\n',OUTPUT.idealExitVelocity,OUTPUT.exitVelocity);
cprintf('black','Efficiency:    %.2f   [%%]\n\n',INPUT.characteristicVelocityEfficiency*100);

cprintf('black','Ideal Tc:      %.2f [K]\nCorrected Tc:  %.2f [K]\n',OUTPUT.idealAdiabaticFlameTemperature(2),OUTPUT.adiabaticFlameTemperature(2))
cprintf('black','%%-Diff.:       %.2f   [%%]\n\n',100-INPUT.characteristicVelocityEfficiency^2*100)

cprintf('black','Ideal Isp:     %.2f  [s]\nCorrected Isp: %.2f  [s]\n',OUTPUT.idealSpecificImpulse,OUTPUT.specificImpulse);
cprintf('black','%%-Diff.:       %.2f   [%%]\n\n',abs(OUTPUT.idealSpecificImpulse - OUTPUT.specificImpulse)/OUTPUT.specificImpulse*100);