function [t_list,X_list,h_avg,num_evals,step_failure_rate] = explicit_RK_variable_step_integration ...
    (rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)

    failed_steps = 0;
    attempted_steps = 0;

    t_start = tspan(1);
    t_end = tspan(2);

    h_curr = h_ref;

    t_list = t_start;
    X_list = X0'; 
    
    total_num_evals = 0;

    while t_list(end) < t_end
        current_t = t_list(end);
        current_X = X_list(end, :)'; 

        % Adjust step to not overshoot t_end
        if current_t + h_curr > t_end
            h_curr = t_end - current_t;
        end

        redo = true;
        while redo
            attempted_steps = attempted_steps + 1;

            [next_X, step_evals, failed, ~, h_next, redo] = ...
                explicit_RK_variable_step(rate_func_in, current_t, current_X, h_curr, BT_struct, p, error_desired);
            
            total_num_evals = total_num_evals + step_evals;
            
            if redo
                failed_steps = failed_steps + failed;
                h_curr = h_next;  
            else
                t_list(end + 1, 1) = current_t + h_curr;
                X_list(end + 1, :) = next_X';
                h_curr = h_next;  
            end
        end
    end

    step_failure_rate = failed_steps / attempted_steps;
    num_evals = total_num_evals;
    h_avg = mean(diff(t_list));
end
