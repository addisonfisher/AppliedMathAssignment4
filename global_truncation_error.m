% function [eror_list] = truncation_error(analytical_soln, t_ref, hspan)
% 
%     h_list = logspace(hspan(1),hspan(2),hspan(3));
%     x_approx_list = [];
% 
%     XA = analytical_soln(t_ref);
%     for i = 1:hspan(3)
% 
%         [x_approx,~] = forward_euler_step(@analytical_soln,t_ref,XA,hspan(i));
%         x_approx_list(end+1) = x_approx;
%     end
%     error_list = x_approx_list-XA;
% 
% end

function [h_list,, g_analytical_difference] = global_truncation_error(tspan, hspan, test_function)

    h_list = logspace(hspan(1), hspan(2), hspan(3));
    t0 = tspan(1);
    tf = tspan(2);

    x_approx_fel_list = []; % Forward Euler Local
    x_approx_expmid_list = []; % Explicit Midpoint
    x_analytical_list = []; % Analytical Solution
    global tot_evals_fel;
    global tot_evals_expmid;
    tot_evals_fel = [];
    tot_evals_expmid = [];

    X0 = solution01(t0);
    for i = 1:(length(h_list)) % subtracting by two to make error_list the same size as t_list in asst3. Fix later.
        
        [~,x_approx_fel, ~, fel_evals] = forward_euler_fixed_step_integration2(test_function, tspan, X0, h_list(i));
        [~,x_approx_expmid,~, expmid_evals] = explicit_midpoint_fixed_step_integration(test_function, tspan, X0, h_list(i));


        x_approx_fel_list(end+1) = x_approx_fel(end);
        x_approx_expmid_list(end+1) = x_approx_expmid(end);
        x_analytical_list(end+1) = solution01(tf);
        tot_evals_fel(end+1) = fel_evals;
        tot_evals_expmid(end+1) = expmid_evals;
    end
    g_fel_error_list = abs(x_approx_fel_list - x_analytical_list);
    g_expmid_error_list = abs(x_approx_expmid_list - x_analytical_list);
    g_analytical_difference = abs(x_analytical_list - x_analytical_list);
    % next step: plot as a function of hspan

    

end