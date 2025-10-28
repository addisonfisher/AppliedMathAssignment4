function test_conserved_quantities()
    m_sun = 1.0;
    m_planet = 3e-6;
    G = 4*pi^2;

    V0 = [1.0; 0.0; 0.0; 2*pi];
    tspan = [0, 10];
    h_ref = 0.01;

    method1_name = 'Strong Stability Preserving Runge-Kutta';
    method1_bt = strong_stability_preserving_RK();
    
    method2_name = 'Explicit Midpoint';
    method2_bt = explicit_midpoint();
    
    method3_name = 'Classic RK4';
    method3_bt = classic_fourth_order();

    results = cell(3, 1);
    my_rate_func = @(t, V) gravity_rate_func(t, V, struct('m_sun', m_sun, 'm_planet', m_planet, 'G', G));

    for i = 1:3
        %run the fixed-step integrator
        if i == 1
            [t_list, X_list, ~, ~] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h_ref, method1_bt);
            method_name = method1_name;
        elseif i == 2
            [t_list, X_list, ~, ~] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h_ref, method2_bt);
            method_name = method2_name;
        else
            [t_list, X_list, ~, ~] = explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h_ref, method3_bt);
            method_name = method3_name;
        end
        
        %calculate the conserved quantities for this run
        [E_list, H_list] = calculate_conserved_quantities(X_list, struct('m_sun', m_sun, 'm_planet', m_planet, 'G', G));
        
        %store the results
        results{i}.name = method_name;
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
    
    for i = 1:3
        %calc error
        E0 = results{i}.E_list(1);
        relative_E_error = (results{i}.E_list - E0) / E0;
        
        plot(results{i}.t_list, abs(relative_E_error), 'LineWidth', 2, 'DisplayName', results{i}.name);
    end
    legend('show', 'Location', 'northwest');
    hold off;
    
    figure();
    hold on;
    grid on;
    title('Conservation of Angular Momentum');
    xlabel('Time (Years)');
    ylabel('Relative Error in Angular Momentum');
    
    for i = 1:3
        H0 = results{i}.H_list(1);
        relative_H_error = (results{i}.H_list - H0) / H0;
        
        plot(results{i}.t_list, abs(relative_H_error), 'LineWidth', 2, 'DisplayName', results{i}.name);
    end
    legend('show', 'Location', 'northwest');

end