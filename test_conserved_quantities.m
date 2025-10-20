function test_conserved_quantities()
    orbit_params = struct();
    orbit_params.m_sun = 1.0;
    orbit_params.m_planet = 3e-6;
    orbit_params.G = 4*pi^2;

    V0 = [1.0; 0.0; 0.0; 2*pi];
    tspan = [0, 10];
    h_ref = 0.01;

    methods = struct();

    methods(1).name = 'Forward Euler';
    methods(1).bt = forward_euler();
    
    methods(2).name = 'Explicit Midpoint';
    methods(2).bt = explicit_midpoint();
    
    methods(3).name = 'Classic RK4';
    methods(3).bt = classic_fourth_order();

    results = cell(length(methods), 1);
    my_rate_func = @(t, V) gravity_rate_func(t, V, orbit_params);

    for i = 1:length(methods)
        
        %run the fixed-step integrator
        [t_list, X_list, ~, ~] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h_ref, methods(i).bt);
        
        %calculate the conserved quantities for this run
        [E_list, H_list] = calculate_conserved_quantities(X_list, orbit_params);
        
        %store the results
        results{i}.name = methods(i).name;
        results{i}.t_list = t_list;
        results{i}.E_list = E_list;
        results{i}.H_list = H_list;
    end

    figure();
    hold on;
    grid on;
    title('Conservation of Mechanical Energy');
    xlabel('Time (Years)');
    ylabel('Relative Error in Energy');
    
    for i = 1:length(methods)
        %calc error
        E0 = results{i}.E_list(1);
        relative_E_error = (results{i}.E_list - E0) / E0;
        
        plot(results{i}.t_list, relative_E_error, 'LineWidth', 2, 'DisplayName', results{i}.name);
    end
    legend('show', 'Location', 'northwest');
    hold off;
    
    figure();
    hold on;
    grid on;
    title('Conservation of Angular Momentum');
    xlabel('Time (Years)');
    ylabel('Relative Error in Angular Momentum');
    
    for i = 1:length(methods)
        H0 = results{i}.H_list(1);
        relative_H_error = (results{i}.H_list - H0) / H0;
        
        plot(results{i}.t_list, relative_H_error, 'LineWidth', 2, 'DisplayName', results{i}.name);
    end
    legend('show', 'Location', 'northwest');

end