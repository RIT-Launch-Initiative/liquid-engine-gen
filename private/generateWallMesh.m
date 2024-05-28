function [Z,R,Zeta,Eta,exteriorEngineContour,wallThickness] = generateWallMesh(interiorEngineContour,chamberThickness,throatThickness,nozzleThickness,n_z,n_r)
%% DETERMINE WALL THICKNESS MATRIX
[~,throatIndex] = min(interiorEngineContour(2,:));

temp = interiorEngineContour(2,:)-interiorEngineContour(2,1);
for i=1:length(temp)
    if temp(i) == 0
        endChamberIndex = i;
    end
end

midNozzleIndex = throatIndex + round((2/3)*(length(interiorEngineContour(1,:)) - throatIndex)); 

initialWallThickness = [interiorEngineContour(1,1) interiorEngineContour(1,endChamberIndex) interiorEngineContour(1,throatIndex) interiorEngineContour(1,midNozzleIndex) interiorEngineContour(1,end);...
                        chamberThickness chamberThickness throatThickness nozzleThickness nozzleThickness];

wallThicknessFit = fit(initialWallThickness(1,:)',initialWallThickness(2,:)','pchipinterp');
wallThickness = wallThicknessFit(interiorEngineContour(1,:))';

%% INITIAL MESH
% Import and initialize mesh points
Z = zeros(n_r,n_z); R = Z;
Z(end,:) = interiorEngineContour(1,:); R(end,:) = interiorEngineContour(2,:);

% Determine radial and axial differentials
dz = (Z(end,end)-Z(end,1))/(n_z-1);

% Set up interior curve fit to maintain the correct shape of geometry
% during mesh calculation.
interior_fit = fit(Z(end,:)',R(end,:)','pchipinterp');

% Spacing the mesh evenly along a 2D vector with constant length on the
% interior curve.
t = (Z(end,1):dz:Z(end,end))./Z(end,end);
meshSpacingXY = interparc(t,Z(end,:),R(end,:))';
meshSpacingXY = meshSpacingXY(:,2:end);
meshSpacingX = [0 meshSpacingXY(1,:)]; meshSpacingY = interior_fit(meshSpacingX);

% Calculate tangent vectors
dx = gradient(meshSpacingX);
dy = gradient(meshSpacingY);
% Calculate orthogonal vectors
orthogonal_dx = -dy;
orthogonal_dy = dx;

for i = 1:length(meshSpacingX)
    % Get derivative at each point
    dfdx(i) = orthogonal_dy(i)/orthogonal_dx(i);
    
    % Calculate orthogonal vector
    orthogonal_vector = [1, dfdx(i)];
    orthogonal_vector = wallThickness(i) * orthogonal_vector / norm(orthogonal_vector);
    
    if(dfdx(i) == -Inf)
        tangentPoint(:,i) = [meshSpacingX(i);interior_fit(meshSpacingX(i))+wallThickness(i)];
    elseif(dfdx(i) > 0)
        tangentPoint(:,i) = [meshSpacingX(i)+orthogonal_vector(1);interior_fit(meshSpacingX(i))+orthogonal_vector(2)];
    elseif(dfdx(i) <= 0)
        tangentPoint(:,i) = [meshSpacingX(i)-orthogonal_vector(1);interior_fit(meshSpacingX(i))-orthogonal_vector(2)];
    end
    Z(:,i) = linspace(meshSpacingX(i),tangentPoint(1,i),n_r);
    R(:,i) = linspace(meshSpacingY(i),tangentPoint(2,i),n_r);
end

%%  Export Full Contour Data to '.csv' in [cm]:

exteriorEngineContour = [Z(end,:);R(end,:)];

z_data=zeros(1,size(exteriorEngineContour,2));
output_matrix=[exteriorEngineContour(1,:);exteriorEngineContour(2,:);z_data]';

outputFolder = 'output';
fileName = 'exteriorEngineContour.txt';
filePath = fullfile(outputFolder,fileName);
writematrix(output_matrix,filePath)

% Outputs file path where the 'contour.csv' file was generated
fprintf('''%s'' generated at %s\n',fileName,filePath);

%% COMPUTATIONAL DOMAIN MESH
eta_all = 1:n_r;
zeta_all = 1:n_z;
% Zeta & Eta Transformantion Matrices
Zeta = repmat(zeta_all,n_r,1);
Eta = repmat(eta_all,n_z,1)';

end