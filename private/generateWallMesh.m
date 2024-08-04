function [Z,R,Zeta,Eta,localWallThickness] = generateWallMesh(interiorEngineContour,exteriorEngineContour,n_z,n_r)
%% INITIAL MESH
interiorFit = fit(interiorEngineContour(1,:)',interiorEngineContour(2,:)','pchipinterp');

% Import and initialize mesh points
Z = zeros(n_r,n_z); R = Z;
Z(end,:) = interiorEngineContour(1,:); R(end,:) = interiorEngineContour(2,:);

%% INTERIOR POINT SPACING
% Spacing the mesh evenly along a 2D vector with constant length on the
% interior curve.
dz = (Z(end,end)-Z(end,1))/(n_z-1);

t = (Z(end,1):dz:Z(end,end))./Z(end,end);
meshSpacingXY = interparc(t,Z(end,:),R(end,:))';
meshSpacingXY = meshSpacingXY(:,2:end);
meshSpacingX = [0 meshSpacingXY(1,:)]; meshSpacingY = interiorFit(meshSpacingX)';
interiorPoints = [meshSpacingX;meshSpacingY];


%% CREATING TANGENT LINES
% Calculate tangent vectors
dx = gradient(interiorPoints(1,:));
dy = gradient(interiorPoints(2,:));
% Calculate orthogonal vectors
ortho_dx = -dy;
ortho_dy = dx;

for i=1:length(interiorPoints)
    % Get the derivative at each point
    dfdx(i) = ortho_dy(i)/ortho_dx(i);

    % Calculate orthogonal vector
    orthogonal_vector = [1, dfdx(i)];
    orthogonal_vector = 1 * orthogonal_vector / norm(orthogonal_vector);
    
    if(dfdx(i) == -Inf)
        [intersect_x(i),intersect_y(i)] = intersections(exteriorEngineContour(1,:),exteriorEngineContour(2,:),[interiorPoints(1,i) interiorPoints(1,i)],[interiorPoints(2,i) interiorPoints(2,i)+1]);

    elseif(dfdx(i) > 0)
        
        [intersect_x(i),intersect_y(i)] = intersections(exteriorEngineContour(1,:),exteriorEngineContour(2,:),[interiorPoints(1,i) interiorPoints(1,i)+orthogonal_vector(1)],[interiorPoints(2,i) interiorFit(interiorPoints(1,i))+orthogonal_vector(2)]);

    elseif(dfdx(i) <= 0)
       [intersect_x(i),intersect_y(i)] = intersections(exteriorEngineContour(1,:),exteriorEngineContour(2,:),[interiorPoints(1,i) interiorPoints(1,i)-orthogonal_vector(1)],[interiorPoints(2,i) interiorFit(interiorPoints(1,i))-orthogonal_vector(2)]);

    end

    Z(:,i) = linspace(interiorPoints(1,i),intersect_x(i),n_r);
    R(:,i) = linspace(interiorPoints(2,i),intersect_y(i),n_r);

end

localWallThickness = abs(R(1,:)-R(end,:));

%% COMPUTATIONAL DOMAIN MESH
eta_all = 1:n_r;
zeta_all = 1:n_z;
% Zeta & Eta Transformantion Matrices
Zeta = repmat(zeta_all,n_r,1);
Eta = repmat(eta_all,n_z,1)';

%pcolor(Z,R,Zeta); axis equal; grid on; grid minor;
end