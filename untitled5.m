% --- 1. Define System Parameters ---
orbit_params = struct();
orbit_params.m_sun = 1.0;
orbit_params.G = 4*pi^2; % Use units where G is convenient

% --- 2. Define Initial Conditions for an Elliptical Orbit ---
% V = [x_p; y_p; vx_p; vy_p]
X0 = [1.0; 0.0; 0.0; 1.5*pi]; % Planet starts at (1,0) with some initial velocity

% --- 3. Define Integration Settings ---
tspan = [0, 5]; % Integrate for 5 time units
h_ref = 0.01;   % Desired step size

% --- 4. Define the Butcher Tableau for RK4 ---
BT_struct_RK4 = struct();
BT_struct_RK4.C = [0; 0.5; 0.5; 1];
BT_struct_RK4.A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
BT_struct_RK4.B = [1/6, 1/3, 1/3, 1/6];

% --- 5. Run the Integration ---
% Make sure gravity_rate_func is on your MATLAB path!
% We create an anonymous function to pass orbit_params correctly
rate_func = @(t, V) gravity_rate_func(t, V, orbit_params);

[t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(...
    rate_func, tspan, X0, h_ref, BT_struct_RK4);

% --- 6. Plot the Results ---
figure;
plot(X_list(:, 1), X_list(:, 2), '-o', 'MarkerSize', 2);
hold on;
plot(0, 0, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 10); % Sun
title('Planetary Orbit Simulation (RK4)');
xlabel('x position');
ylabel('y position');
axis equal;
grid on;
legend('Planet Path', 'Star');