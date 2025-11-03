function planetary_experiments()
    %Orbit params
    orbit_params = struct();
    orbit_params.m_sun = 1;                 
    orbit_params.m_planet = 3e-6;           
    orbit_params.G = 4*pi^2;                

    x0 = 0.4;
    y0 = 0;

    dxdt0 = 0;
    dydt0 = 12;

    V0 = [x0; y0; dxdt0; dydt0];

    t_start = 0;
    t_end = 3; 
    tspan = [t_start, t_end];

    BT_struct = dormand_prince();

    % rate function handle
    my_rate_func = @(t,V) gravity_rate_func(t,V,orbit_params);

    %Analytical Solution
    t_ref_dense = linspace(t_start, t_end, 20000); 
    V_sol = compute_planetary_motion(t_ref_dense, V0, orbit_params);

    % --- ADAPTIVE STEP INTEGRATION EXPERIMENTS ---
    tolerances = logspace(-8, -2, 12);  % desired local error thresholds to test
    n_tol = length(tolerances);

    adapt_avg_h   = zeros(n_tol,1);
    adapt_num_evals = zeros(n_tol,1);
    adapt_max_global_error = zeros(n_tol,1);
    adapt_frac_failed = zeros(n_tol,1);

    for i = 1:n_tol
        tol = tolerances(i);
        fprintf('Adaptive run %d/%d: tol = %.2e\n', i, n_tol, tol);

        h_init = 0.01; 
        p_method = 4;  

        [t_list, X_list, h_avg, num_evals, frac_failed] = ...
            explicit_RK_variable_step_integration(my_rate_func, tspan, V0, h_init, BT_struct, p_method, tol);

        % Interpolate reference onto solver times to compute global error
        X_ref_at_t = zeros(size(X_list));
        for k = 1:4
            X_ref_at_t(:,k) = interp1(t_ref_dense, V_sol(:,k), t_list, 'spline');
        end
        err_vec = vecnorm( X_ref_at_t - X_list, 2, 2 ); % Euclidean error per time
        global_err = max(err_vec);

        adapt_avg_h(i) = h_avg;
        adapt_num_evals(i) = num_evals;
        adapt_max_global_error(i) = global_err;
        adapt_frac_failed(i) = frac_failed;
    end

    % --- FIXED STEP INTEGRATION EXPERIMENTS ---

    h_fixed_list = logspace(-3.5, -1, 14); 
    n_h = length(h_fixed_list);

    fixed_avg_h = zeros(n_h,1);
    fixed_num_evals = zeros(n_h,1);
    fixed_max_global_error = zeros(n_h,1);

    % Using the first row of B 
    
    BT_update_weights = BT_struct;
    BT_update_weights.B = BT_struct.B(1,:);

    for j = 1:n_h
        h_try = h_fixed_list(j);
        fprintf('Fixed run %d/%d: h = %.3e\n', j, n_h, h_try);

        [t_list_f, X_list_f, h_avg_f, num_evals_f] = ...
            explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h_try, BT_update_weights);

        X_ref_at_t_f = zeros(size(X_list_f));

        for k = 1:4
            X_ref_at_t_f(:,k) = interp1(t_ref_dense, V_sol(:,k), t_list_f, 'spline');
        end
        err_vec_f = vecnorm( X_ref_at_t_f - X_list_f, 2, 2 );
        global_err_f = max(err_vec_f);

        fixed_avg_h(j) = h_avg_f;
        fixed_num_evals(j) = num_evals_f;
        fixed_max_global_error(j) = global_err_f;
    end

    %global truncation error vs average step size (loglog)
    figure;
    loglog(adapt_avg_h, adapt_max_global_error, 'ro-','MarkerFaceColor','r','DisplayName','Adaptive RK');
    hold on;
    loglog(fixed_avg_h, fixed_max_global_error, 'bs-','MarkerFaceColor','b','DisplayName','Fixed-step RK');
    xlabel('Average step size (h)'); ylabel('Global max error (max_t ||x_ref-x||)');
    title('Global truncation error vs average step size');
    legend('Location','best'); grid on; hold off;

    %global truncation error vs number of function evaluations (loglog)
    figure;
    loglog(adapt_num_evals, adapt_max_global_error, 'ro-','MarkerFaceColor','r','DisplayName','Adaptive RK');
    hold on;
    loglog(fixed_num_evals, fixed_max_global_error, 'bs-','MarkerFaceColor','b','DisplayName','Fixed-step RK');
    xlabel('Number of function evaluations'); ylabel('Global max error');
    title('Global error vs function evaluations');
    legend('Location','best'); grid on; hold off;

    %failure rate vs average step size
    figure;
    semilogx(adapt_avg_h, adapt_frac_failed, 'ko-','MarkerFaceColor','k');
    xlabel('Average step size (h)'); ylabel('Failure rate (fraction of rejected steps)');
    title('Adaptive step failure rate vs average step size');
    grid on;

    single_tol = 1e-6; % choose a tolerance that shows separate time points
    fprintf('Running single adaptive case for detailed plots, tol = %.2e\n', single_tol);
    [t_list_s, X_list_s, h_avg_s, num_evals_s, frac_failed_s] = ...
        explicit_RK_variable_step_integration(my_rate_func, tspan, V0, 0.01, BT_struct, p_method, single_tol);

    % obtain reference on t_list_s
    X_ref_at_t_s = zeros(size(X_list_s));
    for k = 1:4
        X_ref_at_t_s(:,k) = interp1(t_ref_dense, V_sol(:,k), t_list_s, 'spline');
    end
    % Position vs time plot (line + dots)
    figure;
    subplot(2,1,1);
    plot(t_list_s, X_list_s(:,1), 'ro-','markerfacecolor','k','markersize',2,'LineWidth',1); hold on;
    plot(t_list_s, X_list_s(:,2), 'bo-','markerfacecolor','k','markersize',2,'LineWidth',1);
    plot(t_list_s, X_ref_at_t_s(:,1), 'r--','DisplayName','Ref x'); % optional reference
    title(sprintf('Positions vs time (tol=%.2e) -- red:x, blue:y', single_tol));
    xlabel('t'); ylabel('position'); legend('x_{num}','y_{num}','x_{ref}'); grid on; hold off;

    % Velocity vs time
    subplot(2,1,2);
    plot(t_list_s, X_list_s(:,3), 'go-','markerfacecolor','k','markersize',2,'LineWidth',1); hold on;
    plot(t_list_s, X_list_s(:,4), 'mo-','markerfacecolor','k','markersize',2,'LineWidth',1);
    title('Velocities vs time (dx, dy)');
    xlabel('t'); ylabel('velocity'); legend('dx','dy'); grid on; hold off;

    %Step size vs Radius
    h_list_from_adaptive = diff(t_list_s);
    r_list = sqrt( X_list_s(1:end-1,1).^2 + X_list_s(1:end-1,2).^2 );
    figure;
    semilogy(r_list, h_list_from_adaptive, 'ko','MarkerFaceColor','k');
    xlabel('Distance to sun r'); ylabel('Step size h (log scale)');
    title(sprintf('Adaptive step size vs distance (tol=%.2e)', single_tol));
    grid on;

    fprintf('\n--- Adaptive summary (tol, avg_h, num_evals, max_global_err, frac_failed) ---\n');
    for i = 1:n_tol
        fprintf('tol=%.2e  avg_h=%.3e  evals=%d  max_err=%.3e  frac_failed=%.3e\n', ...
            tolerances(i), adapt_avg_h(i), adapt_num_evals(i), adapt_max_global_error(i), adapt_frac_failed(i));
    end
    fprintf('\n--- Fixed-step summary (h, avg_h, num_evals, max_global_err) ---\n');
    for j = 1:n_h
        fprintf('h=%.3e  avg_h=%.3e  evals=%d  max_err=%.3e\n', ...
            h_fixed_list(j), fixed_avg_h(j), fixed_num_evals(j), fixed_max_global_error(j));
    end
end


