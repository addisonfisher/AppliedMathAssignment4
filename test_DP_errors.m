function test_DP_errors
    %initialize Dormand-Prince method parameters
    BT_struct = dormand_prince();
    t_ref = 0.5;

    rate_func = @rate_func01;
    solution_func = @solution01;

    XA = solution_func(t_ref); 
    h_list = logspace(-8, -1, 100);
    num_steps = length(h_list);

    XB1_results = zeros(1, num_steps);
    XB2_results = zeros(1, num_steps);
    true_solution_list = zeros(1, num_steps); 
    
    %loop through each step size to compute results
    for i = 1:num_steps
        current_h = h_list(i);
        [temp_XB1, temp_XB2, ~] = RK_step_embedded(rate_func, t_ref, XA, current_h, BT_struct);
        XB1_results(i) = temp_XB1;
        XB2_results(i) = temp_XB2;
        true_solution_list(i) = solution_func(t_ref + current_h);
    end
    
    %calculate local truncation errors and differences
    LTE_XB1 = abs(XB1_results - true_solution_list);
    LTE_XB2 = abs(XB2_results - true_solution_list);
    analytical_diff = abs(true_solution_list - XA); 
    proxy_error = abs(XB1_results - XB2_results);



    %filter for XB1
    idx1 = LTE_XB1 > 1e-15;
    [p_XB1, ~] = generate_error_fit(h_list(idx1), LTE_XB1(idx1));
    
    %filter for XB2
    idx2 = LTE_XB2 > 1e-15;
    [p_XB2, ~] = generate_error_fit(h_list(idx2), LTE_XB2(idx2));
    
    %filter for the proxy
    idx_proxy = proxy_error > 1e-15;
    [p_proxy, ~] = generate_error_fit(h_list(idx_proxy), proxy_error(idx_proxy));

    fprintf('XB1 scales with h^p where p = %.4f\n', p_XB1);
    fprintf('XB2 scales with h^p where p = %.4f\n', p_XB2);
    fprintf('Proxy Error scales with h^p where p = %.4f\n', p_proxy);

    figure();
    hold on;
    set(gca, 'XScale', 'log', 'YScale', 'log')
    loglog(h_list, LTE_XB1, 'r.-', 'DisplayName', 'LTE of XB1');
    loglog(h_list, LTE_XB2, 'b.-', 'DisplayName', 'LTE of XB2');
    loglog(h_list, analytical_diff, 'k--', 'DisplayName', '|f(t_ref+h) - f(t_ref)|');
    loglog(h_list, proxy_error, 'g.-', 'DisplayName', 'Error Proxy |XB1 - XB2|');
    title('Dormand-Prince Error Analysis');
    xlabel('Step Size (h)');
    ylabel('Error');
    legend('show', 'Location', 'northwest');
    grid on;
    hold off;

    figure();
    loglog(proxy_error, LTE_XB1, 'r.-', 'DisplayName', 'LTE of XB1');
    hold on;
    loglog(proxy_error, LTE_XB2, 'b.-', 'DisplayName', 'LTE of XB2');
    title('Error vs. Error Proxy |XB1 - XB2|');
    xlabel('Error Proxy |XB1 - XB2|');
    ylabel('Local Truncation Error (LTE)');
    legend('show', 'Location', 'northwest');
    grid on;
    axis equal;
    hold off;
end

function [p,k] = generate_error_fit(x_regression,y_regression)
    %generate Y, X1, and X2
    %note that I use the transpose operator (')
    %to convert the result from a row vector to a column
    %If you are copy-pasting, the ' character may not work correctly
    Y = log(y_regression)';
    X1 = log(x_regression)';
    X2 = ones(length(X1),1);
    %run the regression
    coeff_vec = regress(Y,[X1,X2]);
    %pull out the coefficients from the fit
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end