function [E_list, H_list] = calculate_conserved_quantities(X_list, orbit_params)
    %init necessary lists
    num_steps = size(X_list, 1);
    E_list = zeros(num_steps, 1);
    H_list = zeros(num_steps, 1);

    %get constants
    m_p = orbit_params.m_planet;
    m_s = orbit_params.m_sun;
    G = orbit_params.G;

    %loop through each time step
    for i = 1:num_steps
        %extract position and velocity from the i-th row of X_list
        x = X_list(i, 1);
        y = X_list(i, 2);
        vx = X_list(i, 3);
        vy = X_list(i, 4);

        r = sqrt(x^2 + y^2);

        %calc E
        E_list(i) = 0.5 * m_p * (vx^2 + vy^2) - (m_s * m_p * G / r);

        %calc H
        H_list(i) = m_p * (x * vy - y * vx);
    end
end