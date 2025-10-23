%local truncation error experiments for embedded methods

n_samples = 60;
h_ref_list = logspace(-6,1,n_samples);

abs_diff_list =  zeros(1,n_samples);
t_r_error_list1 = zeros(1,n_samples);
t_r_error_list2 = zeros(1,n_samples);

for n = 1:length(h_ref_list)
        h_ref = h_ref_list(n);
        V_list = compute_planetary_motion(t_space(1) + h_ref, Vo, orbit_params);
        [XB1, XB2, ~ ] = explicitRK_step_embedded(myrate, t_Spacn(1), V0, h_ref, BT_strcut1);
        abs_diff_list(n) = norm(Vlist- V0);
        tr_error_list1 = norm(XB1-Vlist);
        tr_error_list2 = norm(XB2-Vlist);
end

filter_params = struct();
filter_params.minyval = 1e-10;
filter_params.maxyval = 1e-6;


[p1, k1] = loglog_fit(h_ref_list, tr_error_list1, filter_params);
[p2, k2] = loglog_fit(h_ref_list, tr_error_list2, filter_params);


figure(2);
loglog(h_ref_list, abs_diff_list,'ro');
hold on;
loglog(h_ref_list, tr_error_list1, 'bo');

loglog(h_ref_list, tr_error_list2, 'go');


%Notes, 
% for a given method, local trucation error is proportional to the timestep h^p, p is an integer
% 
%     embedded methos spits estimates
% 
%     XB1 is not an element of XB2 for X at the next timestep
% 
%    XB1 is proportional to the timestep h^p1, p1 is an integer
%   XB2 is proportional to the timestep h^p1, p2 is an integer
% 
%   if p2>p1, and h is small, k2h^p2 << k1h^p1
% 
%       So |XB1 - XB2 | is approximately equal to k1h^p1 - which is approximately equal to the local truncation error for XB1
% 

