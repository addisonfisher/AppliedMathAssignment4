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
    
    %Calculte Local Error Between RK and Planetary Motion
    t_ref = 1;
    h_list = logspace(-5, -1, 100);
    x_approx_RK_list = zeros(length(h_list),4); 
    x_analytical_list = zeros(length(h_list),4);

    V_list = compute_planetary_motion(t_ref,V0,orbit_params);
    for i = 1:(length(h_list))
        
        [x_approx_RK, ~] = explicit_RK_step(my_rate_func,t_ref,V_list,h_list(i),BT_struct);
        V_temp = compute_planetary_motion([t_ref, t_ref+h_list(i)],V0,orbit_params);

        x_approx_RK_list(i,:) = x_approx_RK;
        x_analytical_list(i,:) = V_temp(end,:);
    end
    analytical_difference = abs(x_analytical_list - V_list');
    RK_error_list = abs(x_approx_RK_list - x_analytical_list);

    [p1, k1] = loglog_fit(h_list, RK_error_list(:,1));
    yfit1 = k1.*(h_list.^p1);
    [p2, k2] = loglog_fit(h_list, RK_error_list(:,2));
    yfit2 = k2.*(h_list.^p2);
    [p3, k3] = loglog_fit(h_list, RK_error_list(:,3));
    yfit3 = k3.*(h_list.^p3);
    [p4, k4] = loglog_fit(h_list, RK_error_list(:,4));
    yfit4 = k4.*(h_list.^p4);


    figure;
    loglog(h_list, RK_error_list(:,1),'ro','MarkerFaceColor','r'); hold on;
    loglog(h_list,yfit1,'k','LineWidth',2);
    loglog(h_list, RK_error_list(:,2),'bo','MarkerFaceColor','b');
    loglog(h_list,yfit2,'k','LineWidth',2);
    loglog(h_list, RK_error_list(:,3),'go','MarkerFaceColor','g');
    loglog(h_list,yfit3,'k','LineWidth',2);
    loglog(h_list, RK_error_list(:,4),'mo','MarkerFaceColor','m');
    loglog(h_list,yfit4,'k','LineWidth',2);
    title('Local Truncation Error Explicit RK Step for Planetary Motion')
    legend('X','','Y','','dx','','dy')
    hold off
  
%%
    % Global Truncation Error
    % On the same axes, plot the global truncation error as a function of 
    % the average step size for both the fixed a
    % and adaptive step size integrators. This should be a loglog plot. Which performed better?
    x_approx_RK_list = zeros(length(h_list),4); 
    x_analytical_list = zeros(length(h_list),4);
    error_desired = 1e-3;

    init_vals = compute_planetary_motion(t_ref,V0,orbit_params);
    for i = 1:(length(h_list))
        
        [~, x_approx_RK, ~] = explicit_RK_variable_step_integration(my_rate_func,tspan,init_vals,h_list(i),BT_struct,1,error_desired);
        V_temp = compute_planetary_motion(tspan,V0,orbit_params);

        x_approx_RK_list(i,:) = x_approx_RK;
        x_analytical_list(i,:) = V_temp(end,:);
    end
    analytical_difference = abs(x_analytical_list - V_list');
    RK_error_list = abs(x_approx_RK_list - x_analytical_list);


    

end