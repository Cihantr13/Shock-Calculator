% Function to calculate ramp angle given the shock angle
function [theta] = calc_theta(beta_target, M)
    gam = 1.4;

    if beta_target <= 90 && beta_target >= 0
        % Calculate the deflection angle for the given shock angle and Mach number
        theta = atand(2 * cotd(beta_target) * (M.^2 * sind(beta_target).^2 - 1) / (M.^2 * (gam + cosd(2 * beta_target)) + 2));
        
        if theta < 0
            error('Calculated theta value is negative, which is not physically valid! Please check the beta value.');
        end
        
        if beta_target <= 0 || beta_target >= 90
            error('Beta value must be between 30 and 41 degrees.');
        end
        
        if beta_target >= 87
            fprintf('This value is close to a normal shock!\n');
        end
        
        fprintf('Theta value for the given beta: %f degrees\n', theta);
    else
        error('Shock angle out of valid range (0 to 90 degrees)!');
    end
end
