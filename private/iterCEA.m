function [CEARUN,AcAt,AeAt,OF] = iterCEA(OF_override,iterations,AcAt_max,OF_min,OF_max,P_ch_psia,FUEL_type,FUEL_wt,FUEL_h,FUEL_T,FUEL_rho,OXID_type,OXID_wt,OXID_h,OXID_T,OXID_rho,Patm,AeAt_max)
AcAt_mat = linspace(1,AcAt_max,iterations);
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
    Tch(j) = CEARUN.output.eql.temperature(2);
    itercount = itercount + 1;
    progressbar(itercount/total_iterations);
end
% Determine the i,j location of the maximum value within the Cstar output
% matrix.
[Cstar_stoich,ind] = max(Cstar);
OF_stoich = OF_mat(ind);
if(OF_override ~= 0)
    OF = OF_override;
else
    OF = OF_stoich;
end

%% For LaTeX Documentation Only
% Cstar_fit = fit(OF_mat',Cstar','pchipinterp');
% Tch_fit = fit(OF_mat',Tch','pchipinterp');

% t = tiledlayout(figure(2),'flow','TileSpacing','compact');
% t1 = nexttile(t);
% p1 = plot(OF_mat,Cstar,'k'); hold on; p2 = plot(OF,Cstar_fit(OF),'*b','LineWidth',1); p3 = plot(OF_stoich,Cstar_stoich,'or','LineWidth',1);
% 
% 
% ylim([Cstar(1) Cstar_stoich*1.25;])
% yticks(0:200:2000)
% x = t1.XAxis; xticklabels([]);
% ylabel('$\mathrm{c^*}$'); ysecondarylabel('$\mathrm{m/s}$')
% legend([p3, p2],'Stoichiometric O/F','Chosen O/F')
% grid on; grid minor;
% 
% t2 = nexttile(t);
% plot(OF_mat,Tch,'k'); hold on; plot(OF,Tch_fit(OF),'*b','LineWidth',1); plot(OF_stoich,Tch_fit(OF_stoich),'or','LineWidth',1)
% 
% l1 = yline(Tch_fit(OF),'-.b');
% legend(l1,'Max. Allowable $\mathrm{T_0}$');
% 
% ylim([Tch(1) Tch_fit(OF_stoich)*1.25;])
% yticks(0:500:4000)
% xlabel('O/F'); ylabel('$\mathrm{T_{0}}$'); ysecondarylabel('$\mathrm{^\circ{}K}$');
% grid on; grid minor;
% 
% linkaxes([t1,t2],'x');
% 
% exportgraphics(t,'output/OFOptimization.pdf','ContentType','vector')

%% Iterate for AcAt
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

% 95% percentile of performance, relative to minimum.
target = Cstar(1) + 0.95*(Cstar(end)-Cstar(1));
AcAt_fit = fit(Cstar',AcAt_mat','pchipinterp');
AcAt = AcAt_fit(target);

%% For LaTeX Documentation Only
% Cstar_fit = fit(AcAt_mat',Cstar','pchipinterp');
% 
% t = figure(2);
% plot(AcAt_mat,Cstar,'k'); hold on; pt1 = plot(AcAt,Cstar_fit(AcAt),'or','LineWidth',1);
% plot([AcAt AcAt],[0 target],'--r');
% 
% pt2 = plot(13.1707,Cstar_fit(13.1707),'*b','LineWidth',1); % Pre-determined value
% plot([13.1707 13.1707], [0 Cstar_fit(13.1707)],'--b');
% 
% legend([pt1,pt2],'95\% Minimum','1 in. Radius Constraint')
% 
% lim_top = Cstar(end) + 0.25*(Cstar(end) - Cstar(1));
% ylim([Cstar(1) lim_top])
% xlabel('$\mathrm{A_c/A_t}$'); ylabel('$\mathrm{c^*}$'); ysecondarylabel('m/s');
% grid on; grid minor;
% 
% exportgraphics(t,'output/ACATOptimization.pdf','ContentType','vector')

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
end
AeAt_fit = fit(Pe',AeAt_mat','pchipinterp');
AeAt = AeAt_fit(Patm);

%% LaTeX Documentation Only
% Pe = convpres(Pe.*100000,'Pa','psi');
% Pe_fit = fit(AeAt_mat',Pe','pchipinterp');
% AeAt_fit = fit(Pe',AeAt_mat','pchipinterp');
% flowSep1 = AeAt_fit(14.7*0.28);
% flowSep2 = AeAt_fit(14.7*0.38);
% 
% leftBound = [flowSep1 flowSep1; 0 14.7*0.28];
% rightBound = [flowSep2 flowSep2; 14.7*0.38 0];
% lowerBound = [flowSep2 flowSep1; 0 0];
% upperBound = [AeAt_fit(linspace(14.7*0.28,14.7*0.38,100))'; linspace(14.7*0.28,14.7*0.38,100)];
% boundary = [leftBound, upperBound, rightBound, lowerBound];
% 
% t = figure(2);
% plot(AeAt_mat,Pe,'k'); hold on; pt1 = plot(AeAt,Pe_fit(AeAt),'*b','LineWidth',1);
% plot([AeAt AeAt],[0 Pe_fit(AeAt)],'--b');
% 
% boundaryPatch = patch('XData',boundary(1,:),'YData',boundary(2,:),'EdgeColor','r');
% boundaryHatch = hatchfill(boundaryPatch,'single',45,3);
% 
% ylim([0 120])
% xlabel('$\mathrm{A_e/A_t}$'); ylabel('$\mathrm{p_e}$'); ysecondarylabel('psi');
% grid on; grid minor;
% 
% LegendStr = {'Perfect expansion at $\mathrm{p_{atm}}$','Flow separation region: 28\% -- 38\% $\mathrm{p_{atm}}$'};
% [~,object_h,~,~] = legendflex(boundaryPatch, LegendStr(2));
% hatchfill(object_h,'single',45,3);
% 
% 
% exportgraphics(t,'output/AEATOptimization.pdf','ContentType','vector')
% 

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
