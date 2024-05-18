function [CEARUN,AcAt,AeAt,OF] = iterCEA(OF_override,iterations,AcAt_max,OF_min,OF_max,P_ch_psia,FUEL_type,FUEL_wt,FUEL_h,FUEL_T,FUEL_rho,OXID_type,OXID_wt,OXID_h,OXID_T,OXID_rho,Patm,AeAt_max)
AcAt_mat = linspace(1.2,AcAt_max,iterations);
% As AcAt increases toward +inf, c* approaches a theoretical asymptote.
% To design a reasonably sized combustion chamber, the max c* value for the
% maximum AcAt value in the matrix is determined and rounded to the nearest
% 0.01, which significantly reduces the area ratio for the near-exact same
% performance.
AeAt_mat = linspace(1,AeAt_max,iterations);
Patm = Patm*0.0689476; % [bar]

% Initialize the progress bar for runtime estimation
progressbar('iterCEA()...')
itercount = 0;

if(OF_override == 0)
    total_iterations = iterations*3;

    OF_mat = linspace(OF_min,OF_max,iterations);

    % Initialize the output matrix for smaller runtime
    Cstar = NaN(1,iterations);
    for j=1:iterations
        CEARUN=CEA('problem','rocket','equilibrium','fac',...
            'acat',AcAt_mat(1),...
            'o/f',OF_mat(j),...
            'p,psia',P_ch_psia,...
            'reactants','fuel',FUEL_type,'C',3,'H',8,'O',1,...
            'wt%',FUEL_wt,...
            'h,kJ/mol',FUEL_h,...% Standard Enthalpy of Formation from https://webbook.nist.gov/cgi/cbook.cgi?ID=C67630&Mask=1#Thermo-Gas
            't(k)',FUEL_T,...
            'rho,g/cc',FUEL_rho,...
            'oxid',OXID_type,'N',2,'O',1,...
            'wt%',OXID_wt,...
            'h,kJ/mol',OXID_h,... % Standard Enthalpy of Formation from https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=1
            't(k)',OXID_T,...
            'rho,g/cc',OXID_rho,...
            'outp','tran','end');
        Cstar(j) = CEARUN.output.eql.cstar(2);
        itercount = itercount + 1;
        progressbar(itercount/total_iterations);
    end
    % Determine the i,j location of the maximum value within the Cstar output
    % matrix.
    [~,ind] = max(Cstar);
    OF = OF_mat(ind);
else
    total_iterations = iterations*2;
    OF = OF_override;
end

Cstar = NaN(1,iterations);
for i=1:iterations
    CEARUN=CEA('problem','rocket','equilibrium','fac',...
        'acat',AcAt_mat(i),...
        'o/f',OF,...
        'p,psia',P_ch_psia,...
        'reactants','fuel',FUEL_type,'C',3,'H',8,'O',1,...
        'wt%',FUEL_wt,...
        'h,kJ/mol',FUEL_h,...% Standard Enthalpy of Formation from https://webbook.nist.gov/cgi/cbook.cgi?ID=C67630&Mask=1#Thermo-Gas
        't(k)',FUEL_T,...
        'rho,g/cc',FUEL_rho,...
        'oxid',OXID_type,'N',2,'O',1,...
        'wt%',OXID_wt,...
        'h,kJ/mol',OXID_h,... % Standard Enthalpy of Formation from https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=1
        't(k)',OXID_T,...
        'rho,g/cc',OXID_rho,...
        'outp','tran','end');
    Cstar(i) = CEARUN.output.eql.cstar(2);
    itercount = itercount + 1;
    progressbar(itercount/(total_iterations));
end

% Target rounded AcAt for significantly smaller area ratio with near-exact
% performance.
target = floor(Cstar(end)*10)/10;
for i=1:length(Cstar)
    if Cstar(i)-target < 0
        AcAt = AcAt_mat(i);
        break;
    else
        AcAt = AcAt_mat(1);
    end
end

for i=1:iterations
    CEARUN=CEA('problem','rocket','equilibrium','fac',...
        'acat',AcAt,...
        'super',AeAt_mat(i),...
        'o/f',OF,...
        'p,psia',P_ch_psia,...
        'reactants','fuel',FUEL_type,'C',3,'H',8,'O',1,...
        'wt%',FUEL_wt,...
        'h,kJ/mol',FUEL_h,...% Standard Enthalpy of Formation from https://webbook.nist.gov/cgi/cbook.cgi?ID=C67630&Mask=1#Thermo-Gas
        't(k)',FUEL_T,...
        'rho,g/cc',FUEL_rho,...
        'oxid',OXID_type,'N',2,'O',1,...
        'wt%',OXID_wt,...
        'h,kJ/mol',OXID_h,... % Standard Enthalpy of Formation from https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=1
        't(k)',OXID_T,...
        'rho,g/cc',OXID_rho,...
        'outp','tran','end');
    Pe(i) = CEARUN.output.eql.pressure(end);
    itercount = itercount + 1;
    progressbar(itercount/(total_iterations));
    if(Pe(i) - Patm <= 0)
        itercount = total_iterations;
        progressbar(itercount/(total_iterations));
        break
    end
end
AeAt = AeAt_mat(i);

% Run finallized CEA with determined values.
CEARUN=CEA('problem','rocket','equilibrium','fac',...
    'acat',AcAt,...
    'super',AeAt,...
    'o/f',OF,...
    'p,psia',P_ch_psia,...
    'reactants','fuel',FUEL_type,'C',3,'H',8,'O',1,...
    'wt%',FUEL_wt,...
    'h,kJ/mol',FUEL_h,...% Standard Enthalpy of Formation from https://webbook.nist.gov/cgi/cbook.cgi?ID=C67630&Mask=1#Thermo-Gas
    't(k)',FUEL_T,...
    'rho,g/cc',FUEL_rho,...
    'oxid',OXID_type,'N',2,'O',1,...
    'wt%',OXID_wt,...
    'h,kJ/mol',OXID_h,... % Standard Enthalpy of Formation from https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=1
    't(k)',OXID_T,...
    'rho,g/cc',OXID_rho,...
    'outp','tran','end');
end
