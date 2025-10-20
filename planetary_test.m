function planetary_test()

    %call for butcher tableu
    BT_struct = forward_euler();

    %init orbit params
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 3e-6;
    orbit_params.G = 4*pi^2;

    %init initial conditions
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0; y0; dxdt0; dydt0];
    
    %setting timespan
    t_start = 0;
    t_end = 12;
    tspan = [t_start, t_end];
    h_ref = 0.001;

    my_rate_func = @(t_in, V_in) gravity_rate_func(t_in,V_in,orbit_params);
    
    [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h_ref, BT_struct);

    V_true = compute_planetary_motion(t_list, V0, orbit_params);

    figure(); 
    hold on;
    grid on;
    plot(0, 0, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 12); %plot the sun
    plot(X_list(:,1), X_list(:,2), 'b-', 'LineWidth', 1.5); %plot the calculated orbit
    plot(V_true(:,1), V_true(:,2), 'k--'); %plot the true orbit
    title('Planetary Orbit Simulation'); xlabel('X Position'); ylabel('Y Position');
    legend('Sun', 'Simulated Orbit', 'True Orbit', 'Location', 'best'); axis equal;
    hold off;

end