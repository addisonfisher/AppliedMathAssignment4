function planetary_test()

% Forward Euler Tableau
FEBT = struct; % Forward Euler Butcher Table
FEBT.name = 'Forward Euler';
FEBT.a = 0;
FEBT.b = 1;
FEBT.c = 0;

planetary_test_function

[t_list,X_list,h_avg, num_evals] = explicit_midpoint_fixed_step_integration();


end