function planetary_test()

    %call for butcher tableu
    BT_struct = forward_euler();

    %init orbit params
    orbit_params = struct();
    orbit_params.m_sun = 1000;
    orbit_params.m_planet = 10;
    oribt_params.G = 9.8;

    %init initial conditions
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0; y0; dxdt0; dydt0];
    
    %setting timespan
    t_start = 0;
    t_end = 30;
    tspan = [t_start, t_end];
    h_ref = 0.01;

    my_rate_func = @(t_in, V_in) gravity_rate_func(t_in,V_in,orbit_params);
    
    [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(@my_rate_func, tspan, V0, h_ref, BT_struct);

    V_true = compute_planetary_motion(t_list, V0, orbit_params);
end