function [B]= clac_beta1 (theta_target, M)
gam = 1.4;
k = 1;
dt = 1e-5;  % Tolerans değeri


B_min = 0;
B_max = 90;

if theta_target <= 41 && theta_target >= 0

    while (B_max - B_min) > dt
        B(1, k) = (B_min + B_max) / 2;

        tetha = atand(2 * cotd(B(1, k)) * (M.^2 * sind(B(1, k)).^2 - 1) / (M.^2 * (gam + cosd(2 * B(1, k))) + 2));

        if tetha > theta_target
            B_max = B(1, k);
        else
            B_min = B(1, k);
        end

    end
        if B(1, k)>= 87
            fprintf('This values exceed oblique shock value. Generate normal shock!\n')
        end
    fprintf('B value for the given theta:   B(1, i)= %f degrees\n', B(1, k));
else 
    error('Not to cause oblique shock!, Because exceed supersonic flight conditions! ')
end
end