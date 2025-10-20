% Local error for simple ODE (from assignment 3) using current RK step
% function

tref = .5;
hspan = [-5, -1, 100];
BT_struct = forward_euler();

[h_list,analytical_difference,RK_error_list] = truncation_error(tref, hspan, @rate_func01, BT_struct);

figure;
loglog(h_list,RK_error_list,'ro'); hold on;
loglog(h_list,analytical_difference,'bo');