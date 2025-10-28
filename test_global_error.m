function test_global_error()

    rate_func = @rate_func01;
    solution_func = @solution01;
    
    tspan = [0, 5];
    t_start = tspan(1);
    t_end = tspan(2);
    X0 = solution_func(t_start);
    
    X_true_final = solution_func(t_end);

    BT_struct = classic_fourth_order();

    n_points = 20;
    h_list = logspace(-3, 0, n_points);
    
    global_error_list = zeros(n_points, 1);
    num_evals_list = zeros(n_points, 1);

    for i = 1:length(h_list)
        current_h = h_list(i);
        
        [~, X_list, ~, num_evals] = explicit_RK_fixed_step_integration(...
            rate_func, tspan, X0, current_h, BT_struct);
        
        X_approx_final = X_list(end, :)';
        
        global_error_list(i) = norm(X_approx_final - X_true_final);
        num_evals_list(i) = num_evals;
    end

    
    valid_indices_h = find(global_error_list > 1e-15 & isfinite(global_error_list));
    
    if isempty(valid_indices_h)
        error('No valid error data for h plot.');
    end
    
    [p_h, k_h] = loglog_fit(h_list(valid_indices_h), global_error_list(valid_indices_h));
    
    y_fit_h = k_h * (h_list .^ p_h);

    fprintf('Fit vs h complete. Error scales with h^p where p = %.4f\n', p_h);

    
    valid_indices_evals = find(global_error_list > 1e-15 & isfinite(global_error_list) & num_evals_list > 0);
    
    if isempty(valid_indices_evals)
        error('No valid error data for num_evals plot.');
    end
    
    [p_evals, k_evals] = loglog_fit(num_evals_list(valid_indices_evals), global_error_list(valid_indices_evals));
    
    y_fit_evals = k_evals * (num_evals_list .^ p_evals);

    fprintf('Fit vs num_evals complete. Error scales with N^p where p = %.4f\n', p_evals);


    figure;
    hold on;
    
    loglog(h_list, global_error_list, 'bo', ...
           'MarkerFaceColor', 'b', 'DisplayName', 'Global Error');
           
    loglog(h_list, y_fit_h, 'k--', 'LineWidth', 2, ...
           'DisplayName', sprintf('Fit Line (p = %.2f)', p_h));
    
    hold off;
    
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Step Size (h)');
    ylabel('Global Truncation Error at t=5');
    title('Global Error vs. Step Size (Classic RK4)');
    legend('show', 'Location', 'northwest');
    grid on;

    
    figure;
    hold on;
    
    loglog(num_evals_list, global_error_list, 'ro', ...
           'MarkerFaceColor', 'r', 'DisplayName', 'Global Error');
           
    loglog(num_evals_list, y_fit_evals, 'k--', 'LineWidth', 2, ...
           'DisplayName', sprintf('Fit Line (p = %.2f)', p_evals));
           
    hold off;
    
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Number of Function Evaluations');
    ylabel('Global Truncation Error at t=5');
    title('Global Error vs. Function Evaluations (Classic RK4)');
    legend('show', 'Location', 'northwest');
    grid on;


    fprintf('\n--- Global Error vs. Step Size Table ---\n');
    
    table_indices = valid_indices_h(round(linspace(1, length(valid_indices_h), 5)));
    
    header_format = '%-15s | %-20s | %-15s\n';
    fprintf(header_format, 'Step Size (h)', sprintf('Global Error (O(h^%.2f))', p_h), 'Num Evals');
    fprintf(repmat('-', 1, 56 + length(sprintf('%.2f', p_h))), '\n');
    
    row_format = '%-15.4e | %-20.4e | %-15d\n';
    for k = 1:length(table_indices)
        idx = table_indices(k);
        
        h_val = h_list(idx);
        error_val = global_error_list(idx);
        evals_val = num_evals_list(idx);
        
        fprintf(row_format, h_val, error_val, evals_val);
    end
 
end