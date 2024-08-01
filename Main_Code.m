%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                                              %%%%%%%%%%%%%
%%%%%%%%%%%%%                      // Cihan Tiryaki \\                     %%%%%%%%%%%%%
%%%%%%%%%%%%%                    // Shock Calculator  \\                   %%%%%%%%%%%%%
%%%%%%%%%%%%%                                                              %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT: 
%       - Freestream Mach Number (M)
%       - Altitude (h)
%       - Ramp Number
%       - Ramp Angle or Shock Angle
%       - Normal Shock 

% OUTPUT:
%       - Total temperature/pressure variation across the ramps
%       - Static temperature/pressure variation across the ramps
%       - Density variation across the geometry of ramps
%       - Mach Number and Velocity variation across the ramps
%       - TPR (Total Pressure Recovery) variation across the ramps

% FUNCTION USED:
% Atmosphere Model: 
%        - Calculates the properties of air at free-flow altitude
% Calculate Shock Angle: 
%        - Calculates the shock angle by entering the ramp angle using numerical methods
% Main Code : 
%        - Allows compilation of all documents and printing of graphics, geometries and notebooks

clc;
clear all;
close all;

% Initial Parameters
gam = 1.4; % Specific heat ratio
R = 287.05; % Specific gas constant for dry air [J/(kg*K)]

% Prompt the user to select the situation
sys = input(sprintf("Which situation do you want?\n 1: Ramp Angle (theta)\n 2: Shock Angle (beta)\n "));

% Prompt the user for freestream conditions
h = input("Input Freestream Condition (in meters): ");
M = input("Input Mach Number: ");
num_ramps = input("Input Number of Ramps (up to 4): ");

% Arrays for ramp and shock angles
theta_targets = zeros(1, num_ramps); % Array for ramp angles
B = zeros(1, num_ramps); % Array for shock angles

% Get ramp or shock angles from the user
if sys == 1
    for i = 1:num_ramps
        theta_targets(i) = input(sprintf("Input Ramp %d Angle (in degrees): ", i));
    end
elseif sys == 2
    for i = 1:num_ramps
        B(i) = input(sprintf("Input Shock %d Angle (in degrees): ", i));
        theta_targets(i) = calc_theta(B(i), M);
    end
else
    error("Invalid input for situation selection. Choose 1 or 2.");
end

% Prompt the user for the presence of a normal shock wave
normal_shock = input("Is there a normal shock wave? (yes=1, no=0): ");

% Calculate air properties at the given altitude using the atmosphere model
[T, P, a, rho, g] = atmosphere_model1(h);
M_values = M;
T_values = T;
P_values = P;
rho_values = rho;
a_values = a;
V_values = a * M;
TPR_values = 1;

% Calculate total temperature and pressure values
Pt_values = P * (1 + (gam - 1) / 2 * M^2)^(gam / (gam - 1));
Tt_values = T * (1 + (gam - 1) / 2 * M^2);

% Coordinates for ramp and shock waves
x_coords = [0];
y_coords = [0];

% Prepare the figure for plotting
figure;
hold on;
switch sys
    case 1
        for i = 1:num_ramps
            B = clac_beta1(theta_targets(i), M_values(end));
            B_values(i) = B;

            % Normal shock Mach number
            M_n1 = M_values(end) * sind(B);
            M_n2 = sqrt((1 + (gam - 1) / 2 * M_n1^2) / (gam * M_n1^2 - (gam - 1) / 2));
            M_next = M_n2 / sind(B - theta_targets(i));
            M_values = [M_values, M_next];

            % Update pressure and temperature values
            Pt = P_values(end) * (1 + (gam - 1) / 2 * M_values(end)^2)^(gam / (gam - 1));
            P_next = P_values(end) * (2 * gam * M_values(end)^2 - (gam - 1)) / (gam + 1);
            Pt_next = P_next * (((1 + ((gam - 1) / 2) * (M_next^2))) ^ (gam / (gam - 1)));
            T_next = T_values(end) * ((1 + 2 * gam * (M_values(end)^2 - 1) / (gam + 1)) * ((2 + (gam - 1) * (M_values(end)^2)) / ((gam + 1) * (M_values(end)^2))));
            rho_next = rho_values(end) * (((gam + 1) * M_values(end)^2) / ((gam - 1) * M_values(end)^2 + 2));
            a_next = sqrt(gam * R * T_next);
            V_next = a_next * M_next;
            TPR_next = Pt_next / Pt;

            % Update total temperature and pressure values
            Tt_next = T_next * (1 + (gam - 1) / 2 * M_next^2);
            Pt_values = [Pt_values, Pt_next];
            Tt_values = [Tt_values, Tt_next];

            P_values = [P_values, P_next];
            T_values = [T_values, T_next];
            rho_values = [rho_values, rho_next];
            a_values = [a_values, a_next];
            V_values = [V_values, V_next];
            TPR_values = [TPR_values, TPR_next];

            % Calculate coordinates for ramp and shock wave
            x_coords = [x_coords, x_coords(end) + 1];
            y_coords = [y_coords, y_coords(end) + tand(theta_targets(i))];
            plot(x_coords(end-1:end), y_coords(end-1:end), 'k-', 'LineWidth', 2);

            shock_x = [x_coords(end-1), x_coords(end)];
            shock_y = [y_coords(end-1), y_coords(end-1) + tand(B)];
            plot(shock_x, shock_y, '--', 'LineWidth', 2);

            % Text for static values after each shock
            text_offset_y = 0.05 * i; % Vertical offset for each ramp
            text(shock_x(2) + 0.1, shock_y(2) + text_offset_y, ...
                 sprintf('M%d: %.2f\nT%d: %.2f K\nP%d: %.2f Pa\nRho%d: %.2f kg/m^3\nV%d: %.2f m/s\nPt%d: %.2f Pa\nTt%d: %.2f K', ...
                 i, M_values(i+1), i, T_values(i+1), i, P_values(i+1), i, rho_values(i+1), i, V_values(i+1), i, Pt_values(i+1), i, Tt_values(i+1)), ...
                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        end
    case 2
        for i = 1:num_ramps
            M_n1 = M_values(end) * sind(B(i));
            M_n2 = sqrt((1 + (gam - 1) / 2 * M_n1^2) / (gam * M_n1^2 - (gam - 1) / 2));
            M_next = M_n2 / sind(B(i) - theta_targets(i));
            M_values = [M_values, M_next];

            % Update pressure and temperature values
            Pt = P_values(end) * (1 + (gam - 1) / 2 * M_values(end)^2)^(gam / (gam - 1));
            P_next = P_values(end) * (2 * gam * M_values(end)^2 - (gam - 1)) / (gam + 1);
            Pt_next = P_next * (((1 + ((gam - 1) / 2) * (M_next^2))) ^ (gam / (gam - 1)));
            T_next = T_values(end) * ((1 + 2 * gam * (M_values(end)^2 - 1) / (gam + 1)) * ((2 + (gam - 1) * (M_values(end)^2)) / ((gam + 1) * (M_values(end)^2))));
            rho_next = rho_values(end) * (((gam + 1) * M_values(end)^2) / ((gam - 1) * M_values(end)^2 + 2));
            a_next = sqrt(gam * R * T_next);
            V_next = a_next * M_next;
            TPR_next = Pt_next / Pt;

            % Update total temperature and pressure values
            Tt_next = T_next * (1 + (gam - 1) / 2 * M_next^2);
            Pt_values = [Pt_values, Pt_next];
            Tt_values = [Tt_values, Tt_next];

            P_values = [P_values, P_next];
            T_values = [T_values, T_next];
            rho_values = [rho_values, rho_next];
            a_values = [a_values, a_next];
            V_values = [V_values, V_next];
            TPR_values = [TPR_values, TPR_next];

            % Calculate coordinates for ramp and shock wave
            x_coords = [x_coords, x_coords(end) + 1];
            y_coords = [y_coords, y_coords(end) + tand(theta_targets(i))];
            plot(x_coords(end-1:end), y_coords(end-1:end), 'k-', 'LineWidth', 2);

            shock_x = [x_coords(end-1), x_coords(end)];
            shock_y = [y_coords(end-1), y_coords(end-1) + tand(B(i))];
            plot(shock_x, shock_y, '--', 'LineWidth', 2);

            % Text for static values after each shock
            text_offset_y = 0.05 * i; % Vertical offset for each ramp
            text(shock_x(2) + 0.1, shock_y(2) + text_offset_y, ...
                 sprintf('M%d: %.2f\nT%d: %.2f K\nP%d: %.2f Pa\nRho%d: %.2f kg/m^3\nV%d: %.2f m/s\nPt%d: %.2f Pa\nTt%d: %.2f K', ...
                 i, M_values(i+1), i, T_values(i+1), i, P_values(i+1), i, rho_values(i+1), i, V_values(i+1), i, Pt_values(i+1), i, Tt_values(i+1)), ...
                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        end
end

% Values before the first shock wave
text(x_coords(1) - 0.1, y_coords(2), ...
     sprintf('M0: %.2f\nT0: %.2f K\nP0: %.2f Pa\nRho0: %.2f kg/m^3\nV0: %.2f m/s\nPt0: %.2f Pa\nTt0: %.2f K', M, T, P, rho, V_values(1), Pt_values(1), Tt_values(1)), ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

if normal_shock == 1
    % Normal shock wave calculations
    B_normal = 90;
    M_normal = sqrt(((gam - 1) * M_values(end)^2 + 2) / (2 * gam * M_values(end)^2 - (gam - 1)));
    T_normal = T_values(end) * ((1 + 2 * gam * (M_values(end)^2 - 1) / (gam + 1)) * ((2 + (gam - 1) * (M_values(end)^2)) / ((gam + 1) * (M_values(end)^2))));
    P_normal = P_values(end) * (2 * gam * M_values(end)^2 - (gam - 1)) / (gam + 1);
    Pt_normal = P_normal * (1 + (gam - 1) / 2 * M_normal^2)^(gam / (gam - 1));
    rho_normal = rho_values(end) * ((gam + 1) * M_values(end)^2) / ((gam - 1) * M_values(end)^2 + 2);
    a_normal = sqrt(gam * R * T_normal);
    V_normal = a_normal * M_normal;
    Tt_normal = Tt_values(end);

    M_values = [M_values, M_normal];
    T_values = [T_values, T_normal];
    P_values = [P_values, P_normal];
    Pt_values = [Pt_values, Pt_normal];
    Tt_values = [Tt_values, Tt_normal];
    rho_values = [rho_values, rho_normal];
    a_values = [a_values, a_normal];
    V_values = [V_values, V_normal];

    % Plot normal shock wave
    x_coords = [x_coords, x_coords(end)]; % Normal shock at the same x position
    y_coords = [y_coords, y_coords(end) + 1]; % Normal shock going upward
    plot([x_coords(end-1), x_coords(end-1)], [y_coords(end-1), y_coords(end)], 'r-', 'LineWidth', 2);

    % Plot ramp after normal shock
    x_coords = [x_coords, x_coords(end) + 1];
    y_coords = [y_coords(end - 1), y_coords(end - 1)];
    plot(x_coords(end-1:end), y_coords(end-1:end), 'k-', 'LineWidth', 2);

    % Text for normal shock values
    text(x_coords(end) + 0.1, y_coords(end), ...
         sprintf('M%d: %.2f\nT%d: %.2f K\nP%d: %.2f Pa\nRho%d: %.2f kg/m^3\nV%d: %.2f m/s\nPt%d: %.2f Pa\nTt%d: %.2f K', ...
         num_ramps + 1, M_normal, num_ramps + 1, T_normal, num_ramps + 1, P_normal, num_ramps + 1, rho_normal, num_ramps + 1, V_normal, num_ramps + 1, Pt_normal, num_ramps + 1, Tt_normal), ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

%=====================================================================================================%
        % Draw Plot %
%=====================================================================================================%

title('Shock Wave and Ramp Geometry');
xlabel('Ramp Length (m)');
ylabel('Height (m)');
legend('Ramp', 'Shock Wave', 'Normal Shock');
grid on;
axis equal;
hold off;

% Main Figure
fig = figure('Position', [100, 100, 1000, 600]);

% Create tab group
tab_group = uitabgroup(fig, 'Position', [0.05, 0.1, 0.85, 0.85]);

% First tab (Mach Number)
tab1 = uitab(tab_group, 'Title', 'Mach Number');

% Mach Number
ax1 = axes(tab1, 'Position', [0.1, 0.2, 0.35, 0.7]);
title(ax1, 'Mach Number');
xlabel(ax1, 'Position (m)');
ylabel(ax1, 'Mach Number');
grid(ax1, 'on');
hold on;
grid on;
ax1.GridColor = [0, 0, 0]; % Grid line color (black)

% Velocity
ax2 = axes(tab1, 'Position', [0.55, 0.2, 0.35, 0.7]);
title(ax2, 'Velocity');
xlabel(ax2, 'Position (m)');
ylabel(ax2, 'Velocity (m/s)');
grid(ax2, 'on'); 
hold on;
grid on;
ax2.GridColor = [0, 0, 0]; % Grid line color (black)

% Second tab (Total Values)
tab2 = uitab(tab_group, 'Title', 'Total Values');

% Total Pressure
ax3 = axes(tab2, 'Position', [0.1, 0.2, 0.35, 0.7]);
title(ax3, 'Total Pressure');
xlabel(ax3, 'Position (m)');
ylabel(ax3, 'Total Pressure (Pa)');
grid(ax3, 'on'); 
ax3.GridColor = [0, 0, 0]; % Grid line color (black)

% Total Temperature
ax4 = axes(tab2, 'Position', [0.55, 0.2, 0.35, 0.7]);
title(ax4, 'Total Temperature');
xlabel(ax4, 'Position (m)');
ylabel(ax4, 'Total Temperature (K)');
grid(ax4, 'on'); 
ax4.GridColor = [0, 0, 0]; % Grid line color (black)

% Third tab (Static Values)
tab3 = uitab(tab_group, 'Title', 'Static Values');

% Static Pressure
ax5 = axes(tab3, 'Position', [0.1, 0.55, 0.35, 0.35]);
title(ax5, 'Static Pressure');
xlabel(ax5, 'Position (m)');
ylabel(ax5, 'Static Pressure (Pa)');
grid(ax5, 'on'); 
ax5.GridColor = [0, 0, 0]; % Grid line color (black)

% Static Temperature
ax6 = axes(tab3, 'Position', [0.55, 0.55, 0.35, 0.35]);
title(ax6, 'Static Temperature');
xlabel(ax6, 'Position (m)');
ylabel(ax6, 'Static Temperature (K)');
grid(ax6, 'on');
ax6.GridColor = [0, 0, 0]; % Grid line color (black)

% Density
ax7 = axes(tab3, 'Position', [0.1, 0.08, 0.35, 0.35]);
title(ax7, 'Density');
xlabel(ax7, 'Position (m)');
ylabel(ax7, 'Density (kg/m^3)');
grid(ax7, 'on');  
ax7.GridColor = [0, 0, 0]; % Grid line color (black)

% Define the positions for plotting
x_position = linspace(0, length(x_coords) - 1, length(Pt_values));

% Plot Mach Number
plot(ax1, x_position, M_values, 'b-');
legend(ax1, 'Mach Number');

% Plot Velocity
plot(ax2, x_position, V_values, 'b-');
legend(ax2, 'Velocity');

% Plot Total Pressure
plot(ax3, x_position, Pt_values, 'b-');
legend(ax3, 'Total Pressure');
grid(ax3, 'on');

% Plot Total Temperature
plot(ax4, x_position, Tt_values, 'b-');
legend(ax4, 'Total Temperature');
grid(ax4, 'on');

% Plot Static Pressure
plot(ax5, x_position, P_values, 'b-');
legend(ax5, 'Static Pressure');
grid(ax5, 'on');

% Plot Static Temperature
plot(ax6, x_position, T_values, 'b-');
legend(ax6, 'Static Temperature');
grid(ax6, 'on');

% Plot Density
plot(ax7, x_position, rho_values, 'b-');
legend(ax7, 'Density');
grid(ax7, 'on');
