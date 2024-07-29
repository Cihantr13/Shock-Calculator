%% MODELLING OF ATMOSPHERE
%
%
%
%  This moddeleme contains calculating to freestram flow paramaters with variable altitude


function [T , P , a , rho , g ] = atmosphere_model1(h)

    gamma       = 1.4;                  % Specific heat ratio (cp/cv)
    mol         = 0.0289644;            % molar mass of air
    IR          = 8.31446261815324;     % Universal gas constant [J/K*mol]
    radi        = 6371 * 1000;          % Radius of earth [m]

    % Calculating gravitational acceleration with input altitude    
    IGF = 9.780327 * (1+0.0053024*sin(45.5/57.2957)^2 - 0.0000058*2*sin(45.5/57.2957)^2);
    FAC = -3.086 * 10^(-6) * h;
    g = IGF + FAC;

    %%% Different Atmospehere Layouts (Depends Altitude) %%%

    
    
    
    if 0<=h<=11000;

        P0 = 101325;    % Sea level pressure    [N/m2 or Pa]
        T0 = 288.15;    % Sea level temperature [K]
        L0 = 0.0065;    % Standard lapse ratio  [K/m]

        % Variable altitude parameter [1st]
        H = (h*radi) / (h+radi);

        % Temperature [K] at 'h' altitude [1st]
        T = T0-L0*(H-0); 

        % Pressure [N/m2] at 'h' altitude [1st]
        P = P0 * (1  - ((L0 * (H-0))/T0)) ^ ((g*mol)/(IR*L0));

        % Speed of sound[m/s] at 'h' altitude [1st]
        a = sqrt((gamma*IR*T)/mol); 

        % Density[kg/m3] at 'h' altitude [1st]
        rho = (P*mol)/(IR*T);               


        
    elseif h>11000 && h<=20000;

        P0 = 22632.064;     % 11000 m altitude pressure    [N/m2 or Pa]
        T0 = 216.65;        % 11000 m altitude temperature [K]
        L0 = 0;             % Standard lapse ratio         [K/m]

        % Variable altitude parameter [2nd]
        H = (h*radi) / (h+radi);

        % Temperature [K] at 'h' altitude [2nd]
        T = T0-L0*(H-11000);

        % Pressure [N/m2] at 'h' altitude [2nd]
        P = P0*exp(-mol*g*(H-11000))/(IR*T0);

        % Speed of sound[m/s] at 'h' altitude [2nd]
        a = sqrt((gamma*IR*T)/mol);

        % Density[kg/m3] at 'h' altitude [2nd]
        rho = (P*mol)/(IR*T);


        
    elseif h>20000 && h<=32000

        P0 = 5474.88867;    % 20000 m altitude pressure    [N/m2 or Pa]
        T0 = 216.65;        % 20000 m altitude temperature [K]
        L0 = -0.001;        % Standard lapse ratio         [K/m]

        % Variable altitude parameter [3rd]   
        H = (h*radi) / (h+radi);

        % Temperature [K] at 'h' altitude [3rd]
        T = T0-L0*(H-0);

        % Pressure [N/m2] at 'h' altitude [3rd]
        P = P0 * (1  - ((L0 * (H-0))/T0)) ^ ((g*mol)/(IR*L0));

        % Speed of sound[m/s] at 'h' altitude [3rd]
        a = sqrt((gamma*IR*T)/mol);

        % Density[kg/m3] at 'h' altitude [3rd]
        rho = (P*mol)/(IR*T);


end

% Results - freestream paramaters with variable altitude

fprintf("Temperarue =");
disp(T);
fprintf("Pressure =");
disp(P);
fprintf("Speed of Sound =");
disp(a);
fprintf("Density =");
disp(rho);
fprintf("Gravity =");
disp(g);


end