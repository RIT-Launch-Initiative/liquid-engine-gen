function [x_noz,y_noz] = methodOfCharacteristics(Rt,gam,theta_n)
% Functions to be used:
% 1) flowprandtlmeyer(gamma,mach)
%    outputs: [mach, prantdl-meyer angle, mach angle]
% 2) flowisentropic(gamma,mach)
%    outputs: [mach, temp, density, area ratio]
% startup
% % test data
% % gam_c = 1.2926, gam_e = 1.16219
% Rt = 0.2757; gam = 1.2926; Me = 2.7997;
% [~,nu_e,mu_e] = flowprandtlmeyer(gam,Me);
% theta_n = nu_e/2;

% Re = 0.58895; This is nominally what the exit radius should be
% AeAt = 3.999

hold on
nLines = 100;
%% INITIALIZE CHARACTERISTIC POINTS
nPoints = nLines + nLines*(nLines+1)/2;
isWallLocation = zeros(1,nPoints);
isCenterLineLocation = zeros(1,nPoints);
x_location = zeros(1,nPoints);
y_location = zeros(1,nPoints);
mu = zeros(1,nPoints);
nu = zeros(1,nPoints);
M = zeros(1,nPoints);
flowAngle = zeros(1,nPoints);
% Determine the points that are on the wall
j = 1 + nLines;
k = 0;
for i=1:nPoints
    if(i == j + k)
        isWallLocation(i) = 1;
        k = k + 1;
        j = j + nLines - k;
    end
end
% Determine the points that are on the centerline
for i=1:nPoints
    if(i == 1)
        isCenterLineLocation(i) = 1;
    end
    if(i > 1)
        if(isWallLocation(i-1) == 1)
            isCenterLineLocation(i) = 1;
        end
    end
end
%% INITIAL VALUES
% Divide max wall angle into divisions until minimum wall angle of 0 [deg.]
theta_div = zeros(1,nLines);
temp = theta_n - (1+mod(theta_n,1));
delta = temp/(nLines-2);
for i=1:nLines-2
    theta_div(i+1) = theta_div(i) + delta;
end
theta_div(end) = theta_n;
%% METHOD OF CHARACTERISTICS
% Point (a) at the throat
xa = 0;
flowAngle_a = 0.00000000001; % [deg.] Initial small flow angle
nu_a = flowAngle_a; % [deg.] Initial Prandtl-Meyer angle
[~,~,mu_a] = flowprandtlmeyer(gam,nu_a,'nu'); % [-] Mach number at (a)

flowAngle(1) = 0; % [deg.] Flow angle at (1) is 0 [deg.] due to being on a centerline.
nu(1) = flowAngle_a + nu_a; % [deg.]
% Mach number and mach angle at point (1)
[M(1),~,mu(1)] = flowprandtlmeyer(gam,nu(1),'nu');
slope_a = ((flowAngle_a - mu_a) + (flowAngle(1) - mu(1)))/2;

x_location(1) = Rt*tand(90+slope_a);
%plot(x_location(1),y_location(1),'ob');
% Points 2 to nLines+1
for i=2:nLines+1
    % Assign flow and P-M Angle
    if(isWallLocation(i) == 0)
        flowAngle(i) = theta_div(i);
        nu(i) = theta_div(i);
        [M(i),~,mu(i)] = flowprandtlmeyer(gam,nu(i),'nu');
        % Find theta_ax = nu_ax and the Mach angle
        nu_ax = (2*(theta_div(i)) + (flowAngle(i-1) - nu(i-1)))/2;
        [~,~,mu_ax] = flowprandtlmeyer(gam,nu_ax,'nu');
        % Find left and right running characteristic slopes
        tht_b = (flowAngle(i-1)+mu(i-1)+flowAngle(i)+mu(i))/2;
        tht_t = (nu_ax-mu_ax+flowAngle(i)-mu(i))/2;
        [x_location(i),y_location(i)] = returnXYIntersectionPoint(xa,Rt,tht_t,x_location(i-1),y_location(i-1),tht_b);

        %%plot(x_location(i),y_location(i),'ok');
    else
        flowAngle(i) = flowAngle(i-1);
        nu(i) = nu(i-1);
        M(i) = M(i-1);
        mu(i) = mu(i-1);
        % Find left and second left running (both point up), the first is
        % the interpolated wall angle (max for first reflection).
        tht_t = theta_n;
        tht_b = (flowAngle(i-1) + mu(i-1) + flowAngle(i) + mu(i))/2;
        [x_location(i),y_location(i)] = returnXYIntersectionPoint(xa,Rt,tht_t,x_location(i-1),y_location(i-1),tht_b);

        %plot(x_location(i),y_location(i),'or');
    end
end

% Remaining points
j = 0; k = 1;
for i=nLines+2:nPoints
    % If the point lies on the centerline
    if(isCenterLineLocation(i) == 1)
        flowAngle(i) = 0;
        nu(i) = flowAngle(i-(nLines-j)) + nu(i-(nLines-j)) - flowAngle(i);
        [M(i),~,mu(i)] = flowprandtlmeyer(gam,nu(i),'nu');
        tht_t = (flowAngle(i-(nLines-j)) - mu(i-(nLines-j)) + flowAngle(i) - mu(i))/2;
        tht_b = 0;
        [x_location(i),y_location(i)] = returnXYIntersectionPoint(x_location(i-(nLines-j)),y_location(i-(nLines-j)),tht_t,x_location(i-(nLines-j+1)),y_location(i-(nLines-j+1)),tht_b);

        %plot(x_location(i),y_location(i),'ob');
    end
    % If the point is not on the centerline or on the wall
    if(isCenterLineLocation(i) == 0 && isWallLocation(i) == 0)
        flowAngle(i) = theta_div(k);
        nu(i) = (flowAngle(i-(nLines-j)) + nu(i-(nLines-j)) - (flowAngle(i-1) - nu(i-1)))/2;
        [M(i),~,mu(i)] = flowprandtlmeyer(gam,nu(i),'nu');
        tht_t = (flowAngle(i-(nLines-j)) - mu(i-(nLines-j)) + flowAngle(i) - mu(i))/2;
        tht_b = (flowAngle(i-1) + mu(i-1) + flowAngle(i) + mu(i))/2;
        [x_location(i),y_location(i)] = returnXYIntersectionPoint(x_location(i-(nLines-j)),y_location(i-(nLines-j)),tht_t,x_location(i-1),y_location(i-1),tht_b);
        k = k + 1;

        %plot(x_location(i),y_location(i),'ok');
    end
    % If the point is on a wall
    if(isWallLocation(i) == 1)
        flowAngle(i) = flowAngle(i-1);
        nu(i) = nu(i-1);
        M(i) = M(i-1);
        mu(i) = mu(i-1);
        tht_t = flowAngle(i-1);
        tht_b = (flowAngle(i-1) + mu(i-1) + flowAngle(i) + mu(i))/2;
        [x_location(i),y_location(i)] = returnXYIntersectionPoint(x_location(i-(nLines-j)),y_location(i-(nLines-j)),tht_t,x_location(i-1),y_location(i-1),tht_b);
        k = 1;
        j = j + 1;

        %plot(x_location(i),y_location(i),'or');
    end
end

x_noz = []; y_noz = [];
for i=1:length(isWallLocation)
    if(isWallLocation(i))
        x_noz = [x_noz x_location(i)];
        y_noz = [y_noz y_location(i)];
    end
end
%plot(x_noz,y_noz,'--r')
% axis equal; grid on; grid minor;ylim([0 y_noz(end)*2])
%y_noz(end)
