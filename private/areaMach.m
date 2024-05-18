function [AR_x,M_x,T,P,rho] = areaMach(engine_contour,AcAt,At,gam,Tc,Pc,rho_c)
    Mi = flowisentropic(gam,AcAt,'sub');

    Ae = pi*engine_contour(2,end)^2;
    Me = flowisentropic(gam,Ae/At,'sup');

    M_subsonic = linspace(Mi,1,1000);
    M_supersonic = linspace(1,Me,1000);
    M_supersonic = M_supersonic(2:end);
    
    [~,~,~,~,AR_subsonic] = flowisentropic(gam,M_subsonic);
    [~,~,~,~,AR_supersonic] = flowisentropic(gam,M_supersonic);
    
    curveFit_subsonic = fit(AR_subsonic',M_subsonic','linear');
    curveFit_supersonic = fit(AR_supersonic',M_supersonic','linear');
    
    [~,mind] = min(engine_contour(2,:));
    AR_subsonic_x = (pi.*engine_contour(2,1:mind).^2)./At;
    AR_supersonic_x = (pi.*engine_contour(2,mind+1:end).^2)./At;
    M_subsonic_x = feval(curveFit_subsonic, AR_subsonic_x);
    M_supersonic_x = feval(curveFit_supersonic,AR_supersonic_x);
    AR_x = [AR_subsonic_x,AR_supersonic_x];
    M_x = [M_subsonic_x;M_supersonic_x];

    [~,T,P,rho,~] = flowisentropic(gam,M_x);
    T = T*Tc; P = P*Pc; rho = rho*rho_c;

    output_matrix = [engine_contour(1,:)',T];
    outputFolder = 'output';
    filePath = fullfile(outputFolder,'filmTemp.txt');
    writematrix(output_matrix,filePath)
    % Outputs file path where the 'contour.csv' file was generated
    fprintf('''filmTemp.txt'' generated at %s\n',filePath);

    output_matrix = [engine_contour(1,:)',P];
    outputFolder = 'output';
    filePath = fullfile(outputFolder,'localPressure.txt');
    writematrix(output_matrix,filePath)
    % Outputs file path where the 'contour.csv' file was generated
    fprintf('''localPressure.txt'' generated at %s\n',filePath);
    
end