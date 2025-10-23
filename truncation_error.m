function [h_list,analytical_difference,RK_error_list] = truncation_error(t_ref, hspan, test_function,BT_struct)

    h_list = logspace(hspan(1), hspan(2), hspan(3));
    x_approx_RK_list = zeros(length(h_list),1); 
    x_analytical_list = zeros(length(h_list),1);

    V_list = compute_planetary_motion(t_ref,V0,orbit_params);
    for i = 1:(length(h_list)) % subtracting by two to make error_list the same size as t_list in asst3. Fix later.
        
        [x_approx_RK, ~] = explicit_RK_step(test_function,t_ref,V_list,h_list(i),BT_struct);
        V_temp = compute_planetary_motion([t_ref, t_ref+h_list(i)],V0,orbit_params);

        x_approx_RK_list(i,:) = x_approx_RK;
        x_analytical_list(i,:) = V_temp(end,:);
    end
    analytical_difference = abs(x_analytical_list-V_list); % |X(t + h) − X(t)|
    RK_error_list = norm(x_approx_RK_list - x_analytical_list);

end