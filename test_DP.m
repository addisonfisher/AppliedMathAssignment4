function test_DP
    BT_struct = dormand_prince();
    
    %setting up orbit params
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 3e-6;
    orbit_params.G = 4*pi^2;
    
    %creating anon func for gravity rate func
    my_rate_func = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);
    
    %defining timestep inputs
    t = 0;
    h = logspace(1,-3,100);
    num_steps = length(h);
    XA = [8; 0; 0; 1.5];

    XB1 = zeros(4, num_steps);
    XB2 = zeros(4, num_steps);
    diff_results = zeros(4, num_steps);
    num_evals = zeros(1, num_steps);
    
    %calling function
    for i = 1:length(h)
        current_h = h(i);
        [temp_XB1, temp_XB2, temp_num_evals] = RK_step_embedded(my_rate_func, t, XA, current_h, BT_struct);

        XB1(:, i) = temp_XB1;
        XB2(:, i) = temp_XB2;
        diff_results(:, i) = temp_XB1 - temp_XB2;
        num_evals(i) = temp_num_evals;
    end
    
    disp('First estimate (XB1):');
    disp(XB1);
    
    disp('Second estimate (XB2):');
    disp(XB2);
    
    disp('Difference (XB1 - XB2):');
    disp(diff_results);
    
    disp('Number of evaluations:');
    disp(num_evals);

    error_proxy_norm = zeros(1, num_steps);
    for i = 1:num_steps
       error_proxy_norm(i) = norm(diff_results(:, i)); 
    end
    
    figure;
    loglog(h, error_proxy_norm, 'b.-');
    title('Dormand-Prince Error Proxy vs. Step Size');
    xlabel('Step Size (h)');
    ylabel('Error Proxy Norm |XB1 - XB2|');
    grid on;
end